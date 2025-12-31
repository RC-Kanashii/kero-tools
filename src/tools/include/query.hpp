#pragma once

#include <string>
#include <iostream>
#include <vector>
#include <kero-api/kero_mmap.hpp>

#include "kero-api/kero_io.hpp"
#include "sequences.hpp"
#include "CLI/CLI.hpp"
#include "kerotools.hpp"

// declare the function
void uint8_unpacking(uint8_t packed, char * decoded, size_t size);
std::string decode_sequence(uint8_t * encoded, size_t size);
uint8_t uint8_packing(std::string sequence);
void encode_sequence(const std::string& sequence, uint8_t * encoded);


// KmerInfo is optimized to use fixed-size arrays for k-mer and reverse complement k-mer,
struct KmerInfo {
    // Use fixed-size arrays to avoid dynamic allocation overhead of vector
    uint8_t encoded_kmer[8]{};  // Assume maximum k-mer length is 128 (128+3)/4=32 bytes
    uint8_t revcomp_kmer[8]{};
    std::string kmer;
    size_t original_index;
    uint8_t nb_bytes;  // Actually used bytes in the k-mer

    KmerInfo(std::string  k, size_t idx, const uint8_t* enc, const uint8_t* rev, uint8_t bytes)
        : kmer(std::move(k)), original_index(idx), nb_bytes(bytes) {
        memcpy(encoded_kmer, enc, bytes);
        memcpy(revcomp_kmer, rev, bytes);
    }
};

// An optimized hash function for k-mers using a pair of pointers
struct KmerHasher {
    std::size_t operator()(const std::pair<const uint8_t*, uint8_t>& kmer_data) const {
        const uint8_t* data = kmer_data.first;
        uint8_t len = kmer_data.second;
        std::size_t hash = 0;
        for (uint8_t i = 0; i < len; ++i) {
            hash ^= std::hash<uint8_t>{}(data[i]) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

struct KmerEqual {
    bool operator()(const std::pair<const uint8_t*, uint8_t>& lhs,
                   const std::pair<const uint8_t*, uint8_t>& rhs) const {
        return lhs.second == rhs.second && memcmp(lhs.first, rhs.first, lhs.second) == 0;
    }
};

struct SectionData {
    // Section header info
    uint8_t section_minimizer[10]; // Assuming max possible minimizer bytes
    uint8_t nb_bytes_mini;
    uint64_t nb_blocks;
    uint64_t seq_col_offset;

    // Decompressed column data
    std::vector<uint64_t> n_values;
    std::vector<uint64_t> m_idx_values;
    std::vector<uint8_t> data_values;
    std::vector<uint8_t> seq_values;  // For ROW mode: all seq data concatenated

    // Kero parameters
    uint64_t k;
    uint64_t m;
    uint64_t max_value;
    uint64_t data_size;
};

class Query: public KeroTool {
private:
    std::string input_filename;
    std::string output_filename;
    std::string kmer_str;
    std::string kmer_str_filename;
    bool is_no_index{false};
    bool is_gv_read{false};
    bool is_h_proceeded{false};
    bool kmer_per_line{false};
    std::vector<std::string> kmers;
    void prepare();

    std::unique_ptr<kero::Kero_Mmap_Accessor> mmap_accessor; // Use a smart pointer
    // The new query function for parallel execution
    uint64_t query_kmer_parallel(const std::string& kmer_str) const;

    void group_kmers_by_minimizer();
    void process_minimizer_sections_batch(std::vector<std::pair<std::string, uint64_t>>& results);
    void process_single_minimizer_section(
        uint64_t minimizer_val,
        const std::vector<KmerInfo>& kmer_infos,
        std::vector<std::pair<std::string, uint64_t>>& results,
        const kero::Kero_Mmap_Accessor* mmap_accessor);

    void decompress_section_data(
        const uint8_t* section_ptr,
        uint64_t ms_idx,
        SectionData& section_data,
        const kero::Kero_Mmap_Accessor* mmap_accessor);

    // Find k-mers in ROW mode (direct mmap access)
    std::unordered_map<size_t, uint64_t> batch_find_kmers_in_row_mode(
        const std::vector<KmerInfo>& kmers_to_find_batch,
        const uint8_t* section_ptr,
        uint64_t ms_idx,
        const uint8_t* section_minimizer,
        const kero::Kero_Mmap_Accessor* mmap_accessor);

    // Find k-mers in decompressed data (COLUMNAR modes)
    std::unordered_map<size_t, uint64_t> batch_find_kmers_in_decompressed_data(
        const std::vector<KmerInfo>& kmers_to_find_batch,
        const SectionData& section_data);

    void extract_first_kmer(const uint8_t* seq_buffer, uint8_t* extracted_kmer, uint64_t nb_kmer_bytes,
                               uint8_t nb_seq_padding, uint8_t nb_kmer_padding, uint8_t kmer_padding_mask);
    void slide_to_next_kmer(const uint8_t* seq_buffer, uint8_t* extracted_kmer, uint64_t nb_kmer_bytes,
                               int next_nucl_pos, uint8_t kmer_padding_mask);
    void output_results(const std::vector<std::pair<std::string, uint64_t>>& results);
    void cleanup_kmer_infos();

    std::unordered_map<uint64_t, std::vector<KmerInfo>> minimizer_to_kmers;

public:
    void cli_prepare(CLI::App * subapp) override;
    void exec() override;

    void prepareForTest(std::string input_filename, std::string output_filename, std::string kmer_str, bool is_no_index);

    Kero_file* kero_file{};
    Section_Hashtable* sh{};
    MinimizerSearcher* ms{};
    std::vector<uint8_t> nucl_map;

    uint64_t k{};
    uint64_t m{};
    uint64_t max{};
    uint64_t data_size{};

    Query();
    ~Query() override;
    uint64_t get(uint8_t* kmer, uint64_t nb_kmer_bytes) const;
    uint64_t get_with_local_file(uint8_t *kmer, uint64_t nb_kmer_bytes, Kero_file *local_file) const;
    uint64_t get_without_index(uint8_t* kmer, uint64_t nb_kmer_bytes) const;
    uint64_t get_without_index_with_local_file(uint8_t *kmer, uint64_t nb_kmer_bytes, Kero_file *local_file) const;
};

void get_revcomp_kmer(uint8_t* kmer, uint64_t k, uint64_t nb_kmer_bytes, std::vector<uint8_t> nucl_map, uint8_t* result_kmer);