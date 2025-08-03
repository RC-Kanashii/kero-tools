#include <iostream>
#include <sstream>
#include <omp.h>

#include "query.hpp"
#include "logfault.h"

#include "kero-api/kero_mmap.hpp"
#include "kero-api/kero_io.hpp"
#include "kero-api/detail/util.hpp"

char const_nucleotides[4] = {'A', 'C', 'T', 'G'};

/* Decode a packed uint8_t sequence into a char array.
 * The size must be <= 4, as the uint8_t is packed with 2 bits per nucleotide.
 * The nucleotides are decoded in big-endian order.
 */
void uint8_unpacking(uint8_t packed, char *decoded, size_t size) {
    assert(size <= 4);

    size_t offset = 4 - size;
    for (size_t i = 0; i < size; i++) {
        decoded[i] = const_nucleotides[(packed >> ((3 - i - offset) * 2)) & 0b11];
    }
}

std::string decode_sequence(uint8_t *encoded, size_t size) {
    std::stringstream ss;
    char tmp_chars[4] = {0, 0, 0, 0};

    // Decode the truncated first compacted 8 bits
    size_t remnant = size % 4;
    if (remnant > 0) {
        uint8_unpacking(encoded[0], tmp_chars, remnant);
        for (size_t i = 0; i < remnant; i++) {
            ss << tmp_chars[i];
        }
        encoded += 1;
    }

    // Decode all the 8 bits packed
    size_t nb_uint_used = size / 4;
    for (size_t i = 0; i < nb_uint_used; i++) {
        uint8_unpacking(encoded[i], tmp_chars, 4);
        for (char tmp_char: tmp_chars) {
            ss << tmp_char;
        }
    }

    return ss.str();
}

/* Transform a char * sequence into a uint8_t 2-bits/nucl
 * Encoding ACTG
 * Size must be <= 4
 */
uint8_t uint8_packing(std::string sequence) {
    size_t size = sequence.length();
    assert(size <= 4);

    uint8_t val = 0;
    for (size_t i = 0; i < size; i++) {
        val <<= 2;
        val += (sequence[i] >> 1) & 0b11;
    }

    return val;
}

/* Encode the sequence into an array of uint8_t packed sequence slices.
 * The encoded sequences are organised in big endian order.
 */
void encode_sequence(const std::string &sequence, uint8_t *encoded) {
    size_t size = sequence.length();
    // Encode the truncated first 8 bits sequence
    size_t remnant = size % 4;
    if (remnant > 0) {
        encoded[0] = uint8_packing(sequence.substr(0, remnant));
        encoded += 1;
    }

    // Encode all the 8 bits packed
    size_t nb_uint_needed = size / 4;
    for (size_t i = 0; i < nb_uint_needed; i++) {
        encoded[i] = uint8_packing(sequence.substr(remnant + 4 * i, 4));
    }
}


void add_minimizer_standalone(
    uint64_t k, uint64_t m,
    const uint8_t *section_minimizer, uint8_t nb_bytes_mini,
    uint64_t nb_kmer, uint8_t *seq, uint64_t mini_pos);

Query::Query() = default;

Query::~Query() {
    delete kero_file;
    delete sh;
    delete ms;

    cleanup_kmer_infos();
}

void Query::prepare() {
    this->kero_file = new Kero_file(input_filename, "r");

    // Skip the metadata
    this->kero_file->complete_header();

    // Build nucl_map
    if (nucl_map.empty()) {
        this->nucl_map.resize(4);
        nucl_map[this->kero_file->encoding[0]] = this->kero_file->encoding[3];
        nucl_map[this->kero_file->encoding[1]] = this->kero_file->encoding[2];
        nucl_map[this->kero_file->encoding[2]] = this->kero_file->encoding[1];
        nucl_map[this->kero_file->encoding[3]] = this->kero_file->encoding[0];
    }

    assert(kero_file->indexed);

    // Read variable section first
    if (this->is_gv_read == false) {
        for (auto si: kero_file->index) {
            for (auto &it: si->index) {
                if (it.second == 'v') {
                    uint64_t cur_pos = kero_file->current_position;
                    kero_file->jump_to(kero_file->footer->beginning + it.first);
                    Section_GV sgv(kero_file);
                    sgv.close();
                    kero_file->jump_to(cur_pos);

                    this->k = kero_file->global_vars["k"];
                    this->m = kero_file->global_vars["m"];
                    this->max = kero_file->global_vars["max"];
                    this->data_size = kero_file->global_vars["data_size"];
                    this->ms = new MinimizerSearcher(k, m, kero_file->encoding);
                    break;
                }
            }
        }
        this->is_gv_read = true;
    }

    kmers.clear();
    if (!kmer_str_filename.empty()) {
        std::ifstream file(kmer_str_filename);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open kmer file: " + kmer_str_filename);
        }

        std::string line;
        while (std::getline(file, line)) {
            // Skip the annotations and empty lines
            if (line.empty() || line[0] == '#') {
                continue;
            }
            // Remove newline character at the end of the line
            if (!line.empty() && line.back() == '\n') {
                line.pop_back();
            }
            // Capitalize the line
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            if (kmer_per_line) {
                if (line.length() != k) {
                    std::cerr << "Warning: Skipping invalid kmer length: " << line
                            << " (expected " << k << " nucleotides)" << std::endl;
                    continue;
                }
                kmers.push_back(line);
            } else {
                // Split the subsequence into k-mers
                if (line.length() < k) {
                    std::cerr << "Warning: Skipping short sequence: " << line
                            << " (length " << line.length() << " < " << k << ")" << std::endl;
                    continue;
                }
                for (size_t i = 0; i <= line.length() - k; ++i) {
                    kmers.push_back(line.substr(i, k));
                }
            }
        }
        file.close();
    } else {
        // Load k-mers from the command line argument
        if (kmer_str.length() < k) {
            throw std::runtime_error("The substring length is less than k=" + std::to_string(k));
        }
        for (size_t i = 0; i <= kmer_str.length() - k; ++i) {
            kmers.push_back(kmer_str.substr(i, k));
        }
    }

    if (kmers.empty()) {
        throw std::runtime_error("No valid kmers to query");
    }

    // Find the hashtable section
    // Reverse search in the index
    if (this->is_h_proceeded == false) {
        bool find_h = false;
        for (int i = kero_file->index.size() - 1; i >= 0; i--) {
            auto si = kero_file->index[i];
            for (auto it = si->index.rbegin(); it != si->index.rend(); it++) {
                if (it->second == 'h') {
                    uint64_t cur_pos = kero_file->current_position;
                    // Note: the position in the index is a relative position to the start of the footer
                    kero_file->jump_to(kero_file->footer->beginning + it->first);
                    sh = new Section_Hashtable(kero_file);
                    kero_file->jump_to(cur_pos);
                    find_h = true;
                    break;
                }
            }
            if (find_h) {
                break;
            }
        }
        if (!find_h) {
            throw std::runtime_error("No hashtable section found in the index.");
        }
        this->is_h_proceeded = true;
    }
}

void Query::group_kmers_by_minimizer() {
    uint64_t nb_kmer_bytes = (k + 3) / 4;

    // Pre-allocate buffers for encoded and reverse complement kmers
    uint8_t encoded_kmer_buffer[32];
    uint8_t revcomp_kmer_buffer[32];

    for (size_t i = 0; i < kmers.size(); ++i) {
        const std::string &kmer_str = kmers[i];

        // Reset the buffers to zero
        memset(encoded_kmer_buffer, 0, nb_kmer_bytes);
        memset(revcomp_kmer_buffer, 0, nb_kmer_bytes);

        encode_sequence(kmer_str, encoded_kmer_buffer);
        get_revcomp_kmer(encoded_kmer_buffer, k, nb_kmer_bytes, nucl_map, revcomp_kmer_buffer);

        uint64_t minimizer_val;
        bool multi_mini = false;
        uint64_t leftmost_mini = 0;
        ms->kmer_minimizer_compute(encoded_kmer_buffer, nb_kmer_bytes, 0,
                                   minimizer_val, multi_mini, leftmost_mini);
        minimizer_val = kero::mask_mini(minimizer_val, m);

        minimizer_to_kmers[minimizer_val].emplace_back(
            kmer_str, i, encoded_kmer_buffer, revcomp_kmer_buffer, nb_kmer_bytes
        );
    }
}

void Query::cli_prepare(CLI::App *app) {
    this->subapp = app->add_subcommand("query", "Query the k-mer count in the file.");
    auto *kmer_group = subapp->add_option_group("kmer_input", "K-mer input source");
    kmer_group->add_option("-s, --kmer", kmer_str, "The k-mer (or substring) to query.");
    kmer_group->add_option("-S, --kmer-file", kmer_str_filename, "File containing k-mers (or substring) to query.")
            ->check(CLI::ExistingFile);
    subapp->add_flag("--kmer-per-line", kmer_per_line, "Each line in kmer file is a single k-mer.");
    subapp->add_option("-i, --infile", input_filename, "File in which to query the k-mer count.")
            ->required()
            ->check(CLI::ExistingFile);
    subapp->add_flag("--no-index", is_no_index, "Do not use the index to query the k-mer.");
    subapp->add_option("-o, --output", output_filename, "Output file for query results")
            ->required();
}

void Query::exec() {
    // 1. Prepare Phase
    prepare();
    group_kmers_by_minimizer();

    // 2. Parallel Execution Phase
    std::vector<std::pair<std::string, uint64_t> > results(kmers.size());
    process_minimizer_sections_batch(results);

    // 3. Output results
    output_results(results);
}

uint64_t Query::get_without_index(uint8_t *kmer, uint64_t nb_kmer_bytes) const {
    Kero_reader reader = Kero_reader(input_filename);
    auto *tmp_kmer_bytes = new uint8_t[nb_kmer_bytes];
    auto *cnt_bytes = new uint8_t[data_size];
    uint64_t cnt = 0;
    while (reader.has_next()) {
        reader.next_kmer(tmp_kmer_bytes, cnt_bytes);
        // Calculate the reverse complement
        uint8_t revcomp_kmer[nb_kmer_bytes];
        get_revcomp_kmer(kmer, this->k, nb_kmer_bytes, this->nucl_map, revcomp_kmer);
        if (memcmp(tmp_kmer_bytes, kmer, nb_kmer_bytes) == 0 || memcmp(tmp_kmer_bytes, revcomp_kmer, nb_kmer_bytes) ==
            0) {
            for (uint64_t i = 0; i < data_size; i++) {
                cnt = (cnt << 8) | cnt_bytes[i];
            }
            break;
        }
    }

    delete[] tmp_kmer_bytes;
    delete[] cnt_bytes;

    return cnt;
}

void get_revcomp_kmer(uint8_t *kmer, uint64_t k, uint64_t nb_kmer_bytes, std::vector<uint8_t> nucl_map,
                      uint8_t *result_kmer) {
    uint8_t revcomp_kmer[nb_kmer_bytes];
    uint8_t kmer_copy[nb_kmer_bytes];
    memset(revcomp_kmer, 0, nb_kmer_bytes);
    memcpy(kmer_copy, kmer, nb_kmer_bytes);

    // Construct the reverse complete of the k-mer
    for (int i = 0; i < k; i++) {
        uint8_t nucl = kmer_copy[nb_kmer_bytes - 1] & 0b11;
        rightshift8(kmer_copy, nb_kmer_bytes, 2);
        revcomp_kmer[nb_kmer_bytes - 1] |= nucl_map[nucl];
        // std::cout << "revcomp_kmer[0]: " << (short)revcomp_kmer[0] << std::endl;
        if (i < k - 1) {
            leftshift8(revcomp_kmer, nb_kmer_bytes, 2);
        }
    }

    // Mask the first byte of revcomp_kmer
    uint8_t mask = (k % 4) == 0 ? 0xFF : (1 << (k % 4) * 2) - 1;
    revcomp_kmer[0] &= mask;

    std::string kmer_str = decode_sequence(kmer, k);
    std::string revcomp_str = decode_sequence(revcomp_kmer, k);
    // std::cout << "kmer: " << kmer_str << std::endl;
    // std::cout << "revcomp: " << revcomp_str << std::endl;
    memcpy(result_kmer, revcomp_kmer, nb_kmer_bytes);
}

/* This function is a standalone version of the original add_minimizer function
 * adapted to work without relying on the Section_Minimizer class.
 * It takes the minimizer as an argument instead of using 'this->minimizer'.
 */
void add_minimizer_standalone(
    uint64_t k, uint64_t m,
    const uint8_t *section_minimizer, uint8_t nb_bytes_mini,
    uint64_t nb_kmer, uint8_t *seq, uint64_t mini_pos) {
    // This logic is a direct adaptation of Section_Minimizer::add_minimizer
    uint64_t seq_size = nb_kmer + k - 1;
    uint64_t seq_bytes = bytes_from_bit_array(2, seq_size);
    uint64_t seq_left_offset = (4 - (seq_size % 4)) % 4;
    uint64_t no_mini_size = seq_size - m;
    uint64_t no_mini_bytes = bytes_from_bit_array(2, no_mini_size);
    uint64_t no_mini_left_offset = (4 - (no_mini_size % 4)) % 4;
    leftshift8(seq, no_mini_bytes, no_mini_left_offset * 2);

    // Prepare the suffix
    uint8_t *suffix = new uint8_t[seq_bytes];
    memset(suffix, 0, seq_bytes);
    uint suff_nucl = seq_size - m - mini_pos;
    // Values inside seq before any change
    uint no_mini_suff_start_nucl = mini_pos;
    uint no_mini_suff_start_byte = no_mini_suff_start_nucl / 4;
    uint no_mini_suff_bytes = no_mini_bytes - no_mini_suff_start_byte;
    memcpy(suffix, seq + no_mini_suff_start_byte, no_mini_suff_bytes);
    // Shift to the left
    uint no_mini_suff_offset = no_mini_suff_start_nucl % 4;
    leftshift8(suffix, no_mini_suff_bytes, no_mini_suff_offset * 2);


    // Prepare the minimizer
    uint8_t *mini = new uint8_t[seq_bytes];
    memset(mini, 0, seq_bytes);
    memcpy(mini, section_minimizer, nb_bytes_mini);
    // Shift to the left
    uint mini_offset = (4 - (m % 4)) % 4;
    leftshift8(mini, nb_bytes_mini, mini_offset * 2);


    // Align the minimizer
    uint final_mini_start_nucl = mini_pos;
    uint final_mini_start_byte = mini_pos / 4;
    uint final_mini_offset = final_mini_start_nucl % 4;
    uint final_mini_byte_size = (m + final_mini_offset + 3) / 4;
    rightshift8(mini, seq_bytes, final_mini_offset * 2);

    // Merge minimizer
    seq[final_mini_start_byte] = fusion8(
        seq[final_mini_start_byte],
        mini[0],
        final_mini_offset * 2
    );
    for (uint idx = 1; idx < final_mini_byte_size; idx++) {
        seq[final_mini_start_byte + idx] = mini[idx];
    }


    // Align the suffix with the end of the minimizer
    uint final_suff_start_nucl = final_mini_start_nucl + m;
    uint final_suff_start_byte = final_suff_start_nucl / 4;
    uint final_suff_offset = final_suff_start_nucl % 4;
    uint final_suff_byte_size = (suff_nucl + final_suff_offset + 3) / 4;
    if (final_suff_byte_size > 0) {
        rightshift8(suffix, seq_bytes, final_suff_offset * 2);

        // Merge the suffix
        seq[final_suff_start_byte] = fusion8(
            seq[final_suff_start_byte],
            suffix[0],
            final_suff_offset * 2
        );
        for (uint64_t idx = 1; idx < final_suff_byte_size; idx++) {
            seq[final_suff_start_byte + idx] = suffix[idx];
        }
    }

    // Align everything to the right
    rightshift8(seq, seq_bytes, seq_left_offset * 2);

    delete[] suffix;
    delete[] mini;
}

void Query::process_minimizer_sections_batch(std::vector<std::pair<std::string, uint64_t> > &results) {
    mmap_accessor = std::make_unique<kero::Kero_Mmap_Accessor>(input_filename);
    // Group k-mers by minimizer and process them in parallel
    std::vector<uint64_t> minimizer_values;
    minimizer_values.reserve(minimizer_to_kmers.size());
    for (const auto &pair: minimizer_to_kmers) {
        minimizer_values.push_back(pair.first);
    }

#pragma omp parallel for schedule(dynamic) default(none) shared(minimizer_values, results)
    for (size_t i = 0; i < minimizer_values.size(); ++i) {
        uint64_t minimizer_val = minimizer_values[i];
        const auto &kmer_infos = minimizer_to_kmers[minimizer_val];
        process_single_minimizer_section(minimizer_val, kmer_infos, results);
    }
}

void Query::process_single_minimizer_section(
    uint64_t minimizer_val,
    const std::vector<KmerInfo> &kmer_infos,
    std::vector<std::pair<std::string, uint64_t> > &results) {
    uint64_t ms_idx = sh->mpht[minimizer_val];

    const uint8_t *section_ptr = mmap_accessor->get_ptr() + ms_idx;
    uint8_t nb_bytes_mini = (m + 3) / 4;
    uint8_t section_minimizer[nb_bytes_mini];
    memcpy(section_minimizer, section_ptr + 1, nb_bytes_mini);

    if (kero::mask_mini(section_minimizer, m) != minimizer_val) {
        // If the minimizer does not match, set all k-mers to 0
        for (const auto &kmer_info: kmer_infos) {
            results[kmer_info.original_index] = {kmer_info.kmer, 0};
        }
        return;
    }

    // Read and decompress all data columns from this section
    SectionData section_data;
    decompress_section_data(section_ptr, ms_idx, section_data);

    std::unordered_map<size_t, uint64_t> batch_results =
            batch_find_kmers_in_decompressed_data(kmer_infos, section_data);

    // Update results with counts from batch_results
    for (const auto &kmer_info: kmer_infos) {
        uint64_t count = 0;
        // Check if this k-mer is found in the batch results
        if (auto it = batch_results.find(kmer_info.original_index); it != batch_results.end()) {
            count = it->second;
        }
        results[kmer_info.original_index] = {kmer_info.kmer, count};
    }
}

/**
 * @brief Decompress all data columns from a minimizer section at once
 * @param section_ptr Pointer to the start of the section in mmap
 * @param ms_idx Section index in the file
 * @param section_data Output structure to store decompressed data
 */
void Query::decompress_section_data(const uint8_t *section_ptr, uint64_t ms_idx, SectionData &section_data) {
    uint64_t current_offset = 1; // Skip 'M' type marker

    // Copy basic parameters
    section_data.k = this->k;
    section_data.m = this->m;
    section_data.max_value = this->max;
    section_data.data_size = this->data_size;
    section_data.nb_bytes_mini = (m + 3) / 4;

    // Read section minimizer
    memcpy(section_data.section_minimizer, section_ptr + current_offset, section_data.nb_bytes_mini);
    current_offset += section_data.nb_bytes_mini;

    // Read block count and column offsets
    uint64_t n_col_offset, m_idx_col_offset, data_col_offset;
    kero::load_big_endian(section_ptr + current_offset, 8, section_data.nb_blocks);
    current_offset += 8;
    kero::load_big_endian(section_ptr + current_offset, 8, n_col_offset);
    current_offset += 8;
    kero::load_big_endian(section_ptr + current_offset, 8, m_idx_col_offset);
    current_offset += 8;
    kero::load_big_endian(section_ptr + current_offset, 8, data_col_offset);
    current_offset += 8;
    kero::load_big_endian(section_ptr + current_offset, 8, section_data.seq_col_offset);

    // Convert relative offsets to absolute positions
    n_col_offset += ms_idx;
    m_idx_col_offset += ms_idx;
    data_col_offset += ms_idx;
    section_data.seq_col_offset += ms_idx;

    // Initialize vectors with proper sizes
    section_data.n_values.resize(section_data.nb_blocks);
    section_data.m_idx_values.resize(section_data.nb_blocks);

    // Decompress N column (number of k-mers per block)
    uint64_t compressed_size;
    kero::load_big_endian(mmap_accessor->get_ptr() + n_col_offset, 8, compressed_size);
    if (compressed_size > 0) {
        std::vector<uint8_t> compressed_buf(compressed_size);
        memcpy(compressed_buf.data(), mmap_accessor->get_ptr() + n_col_offset + 8, compressed_size);
        // Ensure buffer is aligned to 8 bytes for p4ndec64
        while (compressed_buf.size() % 8 != 0) {
            compressed_buf.push_back(0);
        }
        p4ndec64(compressed_buf.data(), section_data.nb_blocks, section_data.n_values.data());
    }

    // Decompress M_idx column (minimizer positions)
    kero::load_big_endian(mmap_accessor->get_ptr() + m_idx_col_offset, 8, compressed_size);
    if (compressed_size > 0) {
        std::vector<uint8_t> compressed_buf(compressed_size);
        memcpy(compressed_buf.data(), mmap_accessor->get_ptr() + m_idx_col_offset + 8, compressed_size);
        // Ensure buffer is aligned to 8 bytes for p4ndec64
        while (compressed_buf.size() % 8 != 0) {
            compressed_buf.push_back(0);
        }
        p4ndec64(compressed_buf.data(), section_data.nb_blocks, section_data.m_idx_values.data());
    }

    // Decompress Data column (k-mer counts) if present
    if (section_data.data_size > 0) {
        uint64_t data_buf_size;
        kero::load_big_endian(mmap_accessor->get_ptr() + data_col_offset, 8, data_buf_size);
        kero::load_big_endian(mmap_accessor->get_ptr() + data_col_offset + 8, 8, compressed_size);

        if (compressed_size > 0) {
            section_data.data_values.resize(data_buf_size);
            std::vector<uint8_t> compressed_buf(compressed_size);
            memcpy(compressed_buf.data(), mmap_accessor->get_ptr() + data_col_offset + 16, compressed_size);
            // Ensure buffer is aligned to 8 bytes for p4ndec8
            while (compressed_buf.size() % 8 != 0) {
                compressed_buf.push_back(0);
            }
            p4ndec8(compressed_buf.data(), data_buf_size, section_data.data_values.data());
        }
    }
}

/**
 * @brief Extract the first k-mer from a sequence with proper alignment
 */
void Query::extract_first_kmer(const uint8_t *seq_buffer, uint8_t *extracted_kmer, uint64_t nb_kmer_bytes,
                               uint8_t nb_seq_padding, uint8_t nb_kmer_padding, uint8_t kmer_padding_mask) {
    if (nb_seq_padding == nb_kmer_padding) {
        // Direct copy when alignments match
        memcpy(extracted_kmer, seq_buffer, nb_kmer_bytes);
    } else if (nb_seq_padding > nb_kmer_padding) {
        // Sequence has more padding, need to shift left
        uint8_t *tmp_kmer = new uint8_t[nb_kmer_bytes + 1];
        memcpy(tmp_kmer, seq_buffer, nb_kmer_bytes + 1);
        leftshift8(tmp_kmer, nb_kmer_bytes + 1, nb_seq_padding - nb_kmer_padding);
        memcpy(extracted_kmer, tmp_kmer, nb_kmer_bytes);
        delete[] tmp_kmer;
    } else {
        // K-mer has more padding, need to shift right
        memcpy(extracted_kmer, seq_buffer, nb_kmer_bytes);
        rightshift8(extracted_kmer, nb_kmer_bytes, nb_kmer_padding - nb_seq_padding);
    }

    // Apply padding mask
    extracted_kmer[0] &= kmer_padding_mask;
}

/**
 * @brief Slide the k-mer window to extract the next k-mer
 */
void Query::slide_to_next_kmer(const uint8_t *seq_buffer, uint8_t *extracted_kmer, uint64_t nb_kmer_bytes,
                               int next_nucl_pos, uint8_t kmer_padding_mask) {
    // Extract next nucleotide
    uint8_t next_nucl = seq_buffer[next_nucl_pos / 8];
    next_nucl = (next_nucl >> (6 - next_nucl_pos % 8)) & 0b11;

    // Shift k-mer left by 2 bits and add new nucleotide
    leftshift8(extracted_kmer, nb_kmer_bytes, 2);
    extracted_kmer[0] &= kmer_padding_mask;
    extracted_kmer[nb_kmer_bytes - 1] |= next_nucl;
}

// Don't forget to add the output function implementation:
void Query::output_results(const std::vector<std::pair<std::string, uint64_t> > &results) {
    std::ofstream output_file(output_filename);
    if (!output_file.is_open()) {
        throw std::runtime_error("Could not open output file: " + output_filename);
    }

    for (const auto &result: results) {
        output_file << result.first << " " << result.second << std::endl;
    }
    output_file.close();
}

// Also add cleanup function to prevent memory leaks:
void Query::cleanup_kmer_infos() {
    minimizer_to_kmers.clear(); // std::vector<uint8_t> 会自动析构
}


/**
 * @brief Query all target k-mers in the decompressed SectionData.
 * This function traverses the section data only once.
 * @param kmers_to_find_batch All k-mer information corresponding to the current minimizer.
 * @param section_data Decompressed minimizer section data.
 * @return A map from original k-mer index to count.
 */
std::unordered_map<size_t, uint64_t> Query::batch_find_kmers_in_decompressed_data(
    const std::vector<KmerInfo> &kmers_to_find_batch,
    const SectionData &section_data) {
    std::unordered_map<size_t, uint64_t> found_kmer_counts_by_original_index;
    if (kmers_to_find_batch.empty()) {
        return found_kmer_counts_by_original_index;
    }

    const uint64_t k = section_data.k;
    const uint64_t m = section_data.m;
    const uint64_t data_size = section_data.data_size;
    const uint64_t nb_kmer_bytes = (k + 3) / 4;

    std::unordered_map<std::pair<const uint8_t *, uint8_t>, size_t, KmerHasher, KmerEqual> target_kmers_lookup;

    // Prepare a lookup table for target k-mers
    for (const auto &kmer_info: kmers_to_find_batch) {
        target_kmers_lookup[{kmer_info.encoded_kmer, nb_kmer_bytes}] = kmer_info.original_index;
        target_kmers_lookup[{kmer_info.revcomp_kmer, nb_kmer_bytes}] = kmer_info.original_index;
    }

    // Prepare buffers for full sequence and extracted k-mer
    uint64_t max_seq_size = section_data.max_value + k - 1;
    uint64_t seq_buffer_size = bytes_from_bit_array(2, max_seq_size);
    std::unique_ptr<uint8_t[]> full_seq_buffer_uptr(new uint8_t[seq_buffer_size]);
    std::unique_ptr<uint8_t[]> extracted_kmer_uptr(new uint8_t[nb_kmer_bytes]);
    uint8_t *full_seq_buffer = full_seq_buffer_uptr.get();
    uint8_t *extracted_kmer = extracted_kmer_uptr.get();

    uint64_t current_seq_ptr_offset = section_data.seq_col_offset;
    uint64_t current_data_ptr_offset = 0;

    // Traverse all blocks in the section
    for (uint64_t block_idx = 0; block_idx < section_data.nb_blocks; ++block_idx) {
        uint64_t nb_kmers_in_block = section_data.n_values[block_idx];
        if (nb_kmers_in_block == 0) {
            continue;
        }

        uint64_t seq_size_no_mini = nb_kmers_in_block + k - m - 1;
        uint64_t seq_bytes_no_mini = bytes_from_bit_array(2, seq_size_no_mini);

        memset(full_seq_buffer, 0, seq_buffer_size);
        memcpy(full_seq_buffer, mmap_accessor->get_ptr() + current_seq_ptr_offset, seq_bytes_no_mini);

        add_minimizer_standalone(k, m, section_data.section_minimizer, section_data.nb_bytes_mini,
                                 nb_kmers_in_block, full_seq_buffer, section_data.m_idx_values[block_idx]);

        uint8_t nb_seq_padding = ((k + nb_kmers_in_block - 1) % 4) * 2;
        uint8_t nb_kmer_padding = (k % 4) * 2;
        uint8_t kmer_padding_mask = nb_kmer_padding == 0 ? 0xFF : (1 << nb_kmer_padding) - 1;
        int next_nucl_pos = nb_seq_padding + k * 2;

        for (uint64_t kmer_idx = 0; kmer_idx < nb_kmers_in_block; ++kmer_idx) {
            if (kmer_idx == 0) {
                extract_first_kmer(full_seq_buffer, extracted_kmer, nb_kmer_bytes,
                                   nb_seq_padding, nb_kmer_padding, kmer_padding_mask);
            } else {
                slide_to_next_kmer(full_seq_buffer, extracted_kmer, nb_kmer_bytes,
                                   next_nucl_pos, kmer_padding_mask);
                next_nucl_pos += 2;
            }

            auto it = target_kmers_lookup.find({extracted_kmer, nb_kmer_bytes});
            if (it != target_kmers_lookup.end()) {
                uint64_t original_idx = it->second;
                uint64_t count = 0;
                if (data_size > 0 && section_data.data_values.data() != nullptr) {
                    const uint8_t *count_ptr = section_data.data_values.data() + current_data_ptr_offset + (
                                                   kmer_idx * data_size);
                    for (uint64_t c = 0; c < data_size; ++c) {
                        count = (count << 8) | count_ptr[c];
                    }
                }
                found_kmer_counts_by_original_index[original_idx] = count;
            }
        }
        current_seq_ptr_offset += seq_bytes_no_mini;
        current_data_ptr_offset += nb_kmers_in_block * data_size;
    }

    return found_kmer_counts_by_original_index;
}
