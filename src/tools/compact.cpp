#include <utility>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <queue>
#include <cassert>
#include <unordered_set>
#include <memory> // 包含 <memory>

#include "compact.hpp"
#include "encoding.hpp"
#include "sequences.hpp"
#include "merge.hpp"
#include "colinear.hpp"
#include "skmers.hpp"
#include "kero-api/detail/util.hpp"

using namespace std;


Compact::Compact() {
	input_filename = "";
	output_filename = "";
	sorted = false;

	this->k = 0;
	this->m = 0;
	this->bytes_compacted = 0;
	this->offset_idx = 0;
}

Compact::Compact(std::string input_filename, std::string output_filename, bool sorted) {
    this->input_filename = std::move(input_filename);
    this->output_filename = std::move(output_filename);
    this->sorted = sorted;

    this->k = 0;
    this->m = 0;
    this->bytes_compacted = 0;
    this->offset_idx = 0;
}

Compact::~Compact() = default;


void Compact::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("compact", "Read a kero file and try to compact the kmers from minimizer sections. The available ram must be sufficent to load a complete minimizer section into memory.");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "Input kero file to compact.");
	input_option->required();
	input_option->check(CLI::ExistingFile);
	CLI::Option * out_option = subapp->add_option("-o, --outfile", output_filename, "Kero to write (must be different from the input)");
	out_option->required();
	subapp->add_flag("-s, --sorted", sorted, "The output compacted superkmers will be sorted to allow binary search. Sorted superkmer have a lower compaction ratio (ie will be less compacted).");
}


void Compact::exec() {
	Kero_file infile(input_filename, "r");
	Kero_file outfile(output_filename, "w");

	outfile.write_encoding(infile.encoding);

	outfile.set_uniqueness(infile.uniqueness);
	outfile.set_canonicity(infile.canonicity);

	// Metadata transfer
	std::vector<uint8_t> metadata(infile.metadata_size);
	infile.read_metadata(metadata.data());
	outfile.write_metadata(infile.metadata_size, metadata.data());

	bool first_warning = true;

	while (infile.tellp() != infile.end_position) {
		char section_type = infile.read_section_type();

		if (section_type == 'v') {
			Section_GV isgv(&infile);
			isgv.close();

			unordered_map<string, uint64_t> to_copy;
			for (auto & p : isgv.vars) {
				if (p.first != "first_index" and p.first != "footer_size") {
					to_copy[p.first] = p.second;
				}
			}

			if (!to_copy.empty()) {
				Section_GV osgv(&outfile);
				for (auto & p : isgv.vars)
					osgv.write_var(p.first, p.second);
				osgv.close();
			}
		}
		else if (section_type == 'i') {
			Section_Index si(&infile);
			si.close();
		}
		else if (section_type == 'r') {
			if (first_warning) {
				first_warning = false;
				cerr << "WARNING: kero-tools has detected R sections inside of the file. The compact tool is only compacting kmers inside of M sections. The R sections are omitted." << endl;
			}

			Section_Raw sr(&infile);
			sr.close();
		}
		else if (section_type == 'm' or section_type == 'M') {
			uint k = outfile.global_vars["k"];
			uint m = outfile.global_vars["m"];

			if (outfile.global_vars["max"] < k - m + 1) {
				unordered_map<string, uint64_t> values(outfile.global_vars);
				Section_GV sgv(&outfile);

				for (auto & p : values)
					if (p.first != "max")
						sgv.write_var(p.first, p.second);
				sgv.write_var("max", k - m + 1);

				sgv.close();
			}

			Section_Minimizer sm(&infile);
			this->compact_section(sm, outfile);
			sm.close();
		} else if (section_type == 'h') {
			Section_Hashtable sh(&infile);
			sh.close();
		} else {
			cerr << "Unknown section type " << section_type << " in file " << infile.filename << endl;
			exit(2);
		}
	}

	infile.close();
	outfile.close();
}

/* Compact the kmers from a Section_Minimizer and write the result to the output file.
 * The compacted kmers are stored in paths, which are then written to the output file.
 */
void Compact::compact_section(Section_Minimizer & ism, Kero_file & outfile) {
	// General variables
	uint k = outfile.global_vars["k"];
	uint m = outfile.global_vars["m"];
	uint data_size = outfile.global_vars["data_size"];

	this->k = k;
	this->m = m;
	this->data_size = data_size;
	this->bytes_compacted = (k - m + 3) / 4;
    this->mini_pos_size = (static_cast<uint>(ceil(log2(ism.k + ism.max - ism.m))) + 7) / 8;
	this->offset_idx = (4 - ((k - m) % 4)) % 4;

    // 1 - Load the input section into a memory-safe structure
	KmerMatrix kmers_per_index = this->prepare_kmer_matrix(ism);

	// 2 - Compact kmers
	std::vector<std::vector<const uint8_t*>> paths;
	if (this->sorted) {
        this->mini_pos_size = (static_cast<uint>(ceil(log2(ism.k + ism.max - ism.m))) + 7) / 8;
		paths = this->sorted_assembly(kmers_per_index);
	} else {
		std::vector<std::pair<const uint8_t*, const uint8_t*>> to_compact = this->greedy_assembly(kmers_per_index);
		paths = this->pairs_to_paths(to_compact);
	}

	Section_Minimizer osm(&outfile);
	osm.write_minimizer(ism.minimizer);
	this->write_paths(paths, osm, data_size);
	osm.close();
}

/* Prepare the kmer matrix from the Section_Minimizer.
 * The kmer matrix is a vector of vectors, where each inner vector contains
 * unique_ptrs to KmerData objects representing the kmers.
 */
KmerMatrix Compact::prepare_kmer_matrix(Section_Minimizer & sm) {
	KmerMatrix kmer_matrix(sm.k - sm.m + 1);

	uint64_t max_nucl = sm.k + sm.max - 1 - sm.m;
	uint64_t max_seq_bytes = (max_nucl + 3) / 4;
    uint64_t kmer_bytes = (sm.k - sm.m + 3) / 4;
    uint64_t mini_pos_size = (static_cast<uint>(ceil(log2(sm.k + sm.max - sm.m))) + 7) / 8;

	std::vector<uint8_t> seq_buffer(max_seq_bytes);
	std::vector<uint8_t> data_buffer(sm.data_size * sm.max);

    const size_t kmer_storage_size = kmer_bytes + sm.data_size + mini_pos_size;

	// 1 - Load the input section
	for (uint n = 0; n < sm.nb_blocks; n++) {
		uint64_t mini_pos = 0xFFFFFFFFFFFFFFFF;
		// Read sequence
		uint nb_kmers = sm.read_compacted_sequence_without_mini(
			seq_buffer.data(), data_buffer.data(), mini_pos);

		// Add kmer by index
		for (uint kmer_idx = 0; kmer_idx < nb_kmers; kmer_idx++) {
			uint kmer_pos_in_matrix = sm.k - (uint)sm.m - mini_pos + kmer_idx;

			// Create a unique_ptr to manage the kmer storage
            KmerDataPtr kmer_storage = std::make_unique<uint8_t[]>(kmer_storage_size);
            uint8_t* buffer_ptr = kmer_storage.get();

			// Copy kmer sequence
			subsequence(seq_buffer.data(), sm.k - sm.m + nb_kmers - 1, buffer_ptr, kmer_idx, kmer_idx + sm.k - sm.m - 1);

            // Copy data array
			memcpy(buffer_ptr + kmer_bytes, data_buffer.data() + kmer_idx * sm.data_size, sm.data_size);

            // Write mini position
			uint64_t kmer_mini_pos = mini_pos - kmer_idx;
            uint8_t* mini_pos_ptr = buffer_ptr + kmer_bytes + sm.data_size;
            for (int b = mini_pos_size - 1; b >= 0; b--) {
				*(mini_pos_ptr + b) = kmer_mini_pos & 0xFF;
				kmer_mini_pos >>= 8;
			}

			kmer_matrix[kmer_pos_in_matrix].push_back(std::move(kmer_storage));
		}
	}
	return kmer_matrix;
}


uint Compact::mini_pos_from_buffer(const uint8_t * kmer) const {
	uint mini_pos = 0;
	const uint8_t* mini_pos_ptr = kmer + this->bytes_compacted + this->data_size;
	for (uint i = 0; i < this->mini_pos_size; i++) {
		mini_pos <<= 8;
		mini_pos += mini_pos_ptr[i];
	}
	return mini_pos;
}


int Compact::interleaved_compare_kmers(const uint8_t * kmer1, const uint8_t * kmer2) const {
	const uint mini_pos1 = this->mini_pos_from_buffer(kmer1);
	const uint mini_pos2 = this->mini_pos_from_buffer(kmer2);

	assert(mini_pos1 == mini_pos2);

	const uint used_nucl = this->k - this->m;
	const uint offset_nucl = (4 - (used_nucl % 4)) % 4;
	const uint pref_nucl = mini_pos1;
	const uint pref_bytes = (offset_nucl + pref_nucl + 3) / 4;
	const uint suff_nucl = used_nucl - pref_nucl;
	const uint suff_bytes = (suff_nucl + 3) / 4;
	const uint total_bytes = (used_nucl + 3) / 4;

	// --- Prefix ---
	int last_prefix_divergence = -1;
	const uint pref_start_mask = (1u << (8 - 2 * offset_nucl)) - 1;
	const uint pref_stop_mask = ~((1u << (2 * (suff_nucl % 4))) - 1);
	for (uint pref_byte=0 ; pref_byte<pref_bytes ; pref_byte++) {
		uint8_t byte1 = kmer1[pref_byte];
		uint8_t byte2 = kmer2[pref_byte];
		bool end_byte = pref_byte == pref_bytes-1;
		if (pref_byte == 0) {
			byte1 &= pref_start_mask;
			byte2 &= pref_start_mask;
		}
		if (end_byte) {
			byte1 &= pref_stop_mask;
			byte2 &= pref_stop_mask;
		}
		uint8_t result = byte1 xor byte2;
		for (uint8_t i=0 ; i<4 ; i++) {
			if ((not end_byte) or (end_byte and (i >= (suff_nucl % 4)))) {
				if ((result & (0b11 << (2 * i))) != 0) {
					last_prefix_divergence = pref_byte * 4 + 3 - i - offset_nucl;
					break;
				}
			}
		}
	}

	// --- Suffix ---
	int first_suffix_divergence = suff_nucl;
	int current_divergence_idx = 0;
	const uint suff_first_byte = total_bytes - suff_bytes;
	for (uint suff_byte=suff_first_byte ; suff_byte<total_bytes and first_suffix_divergence==(int)suff_nucl ; suff_byte++) {
		uint8_t byte1 = kmer1[suff_byte];
		uint8_t byte2 = kmer2[suff_byte];
		uint8_t result = byte1 xor byte2;
        for (uint8_t i=4 ; i>0 ; i--) {
            if (suff_byte == suff_first_byte and i == 4) {
                i = ((suff_nucl - 1) % 4) + 1;
			}
			if ((result & (0b11 << (2 * (i-1)))) != 0) {
				first_suffix_divergence = current_divergence_idx;
                break;
			} else {
				current_divergence_idx += 1;
			}
		}
	}

	if (last_prefix_divergence == -1 and first_suffix_divergence == (int)suff_nucl) {
		return 0;
	}
	else {
		uint pref_div_distance = pref_nucl - last_prefix_divergence - 1;
		uint nucl_pos = offset_nucl;
		if (pref_div_distance == pref_nucl) {
			nucl_pos += pref_nucl + first_suffix_divergence;
		} else if (first_suffix_divergence == (int)suff_nucl) {
			nucl_pos += last_prefix_divergence;
		}
		else if ((int)pref_div_distance < first_suffix_divergence) {
			nucl_pos += last_prefix_divergence;
		}
		else {
			nucl_pos += pref_nucl + first_suffix_divergence;
		}
		uint byte_pos = nucl_pos / 4;
		uint nucl_shift = 2 * (3 - (nucl_pos % 4));
		uint nucl1 = (kmer1[byte_pos] >> nucl_shift) & 0b11;
		uint nucl2 = (kmer2[byte_pos] >> nucl_shift) & 0b11;
		if (nucl1 < nucl2) return -1;
		else return +1;
	}
}

/* Sort the kmer matrix.
 * Each column of the kmer matrix is sorted using a custom comparator that compares
 * the interleaved representation of the kmers.
 */
void Compact::sort_matrix(KmerMatrix & kmer_matrix) {
	for (uint i = 0; i < kmer_matrix.size(); i++) {
		auto comp_function = [this](const KmerDataPtr & kmer1, const KmerDataPtr & kmer2) {
			return this->interleaved_compare_kmers(kmer1.get(), kmer2.get()) < 0;
		};
		sort(kmer_matrix[i].begin(), kmer_matrix[i].end(), comp_function);
	}
}

/* Convert pairs of kmers into a vector of pairs of indices.
 * The pairs are formed by comparing the kmer sequences in two columns.
 * The function returns a vector of pairs, where each pair contains the indices
 * of the matching kmers from the two columns.
 */
std::vector<std::pair<uint64_t, uint64_t>> Compact::pair_kmers(const std::vector<KmerDataPtr>& column1, const std::vector<KmerDataPtr>& column2) const {
	const uint nb_nucl = k - m;
    vector<pair<uint64_t, uint64_t>> pairs;
	pairs.reserve(max(column1.size(), column2.size()));

	unordered_map<uint64_t, vector<uint64_t>> index;
	for (uint64_t idx = 0; idx < column2.size(); idx++) {
		const uint8_t* kmer = column2[idx].get();
		uint64_t hash = subseq_to_uint(kmer, nb_nucl, 0, nb_nucl - 2);
		index[hash].push_back(idx);
	}

	for (uint64_t idx = 0; idx < column1.size(); idx++) {
		const uint8_t* kmer = column1[idx].get();
		uint64_t hash = subseq_to_uint(kmer, nb_nucl, 1, nb_nucl - 1);
		if (index.count(hash)) {
			for (uint64_t candidate_vect_idx : index.at(hash)) {
				const uint8_t* candidate = column2[candidate_vect_idx].get();
                if (sequence_compare(
							candidate, nb_nucl, 0, nb_nucl - 2,
							kmer, nb_nucl, 1, nb_nucl - 1
						) == 0) {
					pairs.emplace_back(idx, candidate_vect_idx);
				}
			}
		}
	}
	return pairs;
}


/* * A hash function for pairs of integers.
 * This is used to create a hash for pairs of kmers.
 */
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);

        // Mainly for demonstration purposes, i.e. works but is overly simple
        // In the real world, use sth. like boost.hash_combine
        return h1 ^ h2;
    }
};


template<>
bool std::operator<(const PairInt& l, const PairInt& r)
{
	if (l.second == r.second)
		return l.first < r.first;
	else
		return l.second < r.second;
}

ostream& operator<<(ostream& os, const PairInt& p)
{
    os << p.first << ", " << p.second;
    return os;
}


/* This function takes a vector of candidate pairs and computes the longest chain of colinear pairs.
 * It uses the Colinear class to compute scores and find the longest chain.
 */
vector<pair<uint64_t, uint64_t> > Compact::colinear_chaining(vector<pair<uint64_t, uint64_t> > & candidates) const {

    if (candidates.empty()) {
        return {};
    }

    Colinear colinear(candidates);
    colinear.compute_scores(candidates);

	return colinear.longest_chain();
}



/**
 * @brief Finalizes the sorting of k-mers by assembling them into super-kmers and ordering them lexicographically.
 * @param matrix The k-mer matrix, sorted within columns.
 * @param pairs A vector of links between k-mers in adjacent columns.
 * @return A vector where each inner vector is a "super-kmer" (an ordered sequence of k-mers).
 */
std::vector<std::vector<const uint8_t*>> Compact::polish_sort(const KmerMatrix & matrix, const std::vector<std::vector<std::pair<uint64_t, uint64_t>>> & pairs) const {
    uint8_t encoding[4] = {0, 1, 3, 2};
    Stringifyer strif(encoding);

    // Map each k-mer pointer to its column index.
    std::unordered_map<const uint8_t *, size_t> columns;
    for (size_t col = 0; col < matrix.size(); col++) {
        for (const auto& kmer_ptr : matrix[col]) {
            columns[kmer_ptr.get()] = col;
        }
    }

    using SuperKmerPtr = std::unique_ptr<std::vector<const uint8_t*>>;
    std::unordered_map<const uint8_t *, SuperKmerPtr*> rev_skmers_index;
    std::vector<SuperKmerPtr> skmers;

	// Build super kmers from pairs
    for (size_t col_idx = 0; col_idx < pairs.size(); col_idx++) {
        const std::vector<std::pair<uint64_t, uint64_t>>& column = pairs[col_idx];

        for (const auto& link : column) {
            const uint8_t* lkmer = matrix[col_idx][link.first].get();
            const uint8_t* rkmer = matrix[col_idx + 1][link.second].get();

        	// Check if the left kmer is already in a super kmer
            auto it = rev_skmers_index.find(lkmer);
            if (it == rev_skmers_index.end()) {
            	// Create a new super kmer with the left kmer
                auto new_sk = std::make_unique<std::vector<const uint8_t*>>();
                new_sk->push_back(lkmer);
                new_sk->push_back(rkmer);

                SuperKmerPtr* sk_ptr = &skmers.emplace_back(std::move(new_sk));
                rev_skmers_index[rkmer] = sk_ptr;
            } else {
            	// Extend the existing super kmer with the right kmer
                SuperKmerPtr* sk_ptr = it->second;
                (*sk_ptr)->push_back(rkmer);

                rev_skmers_index.erase(it);
                rev_skmers_index[rkmer] = sk_ptr;
            }
        }
    }

	// Create a complete index from any k-mer to its containing super-kmer.
    std::unordered_map<const uint8_t *, SuperKmerPtr*> sk_index;
    for (const auto& kv : rev_skmers_index) {
        SuperKmerPtr* sk_ptr = kv.second;
        for (const uint8_t* kmer : **sk_ptr) {
            sk_index[kmer] = sk_ptr;
        }
    }

	// Build super kmers for individual kmers that are not part of any pair
    for (size_t col = 0; col < matrix.size(); col++) {
        for (const auto& kmer_ptr : matrix[col]) {
            const uint8_t* kmer = kmer_ptr.get();
            if (sk_index.find(kmer) == sk_index.end()) {
                auto lonely_sk = std::make_unique<std::vector<const uint8_t*>>();
                lonely_sk->push_back(kmer);

                SuperKmerPtr* sk_ptr = &skmers.emplace_back(std::move(lonely_sk));
                sk_index[kmer] = sk_ptr;
            }
        }
    }

	// Initialize super kmer counts
    std::unordered_map<SuperKmerPtr*, uint64_t> sk_counts;
    for (const auto& col : matrix) {
        if (col.empty()) continue;

        auto it = sk_index.find(col[0].get());
        if (it != sk_index.end()) {
            SuperKmerPtr* sk_ptr = it->second;
            sk_counts[sk_ptr]++;
        }
    }

    std::vector<size_t> current_kmers(matrix.size(), 0);
    std::vector<std::unique_ptr<uint8_t[]>> memory;
    std::vector<interleved_t> current_interleaves;

    size_t num_kmer_nucl = this->k - this->m;

	// Initialize memory and interleaved structures
    for (size_t i = 0; i < matrix.size(); i++) {
        memory.push_back(std::make_unique<uint8_t[]>((num_kmer_nucl + 3) / 4));
        current_interleaves.push_back({nullptr, 0, 0});

        if (!matrix[i].empty()) {
            current_interleaves[i] = interleaved(matrix[i][0].get(), memory[i].get(),
                                               num_kmer_nucl, num_kmer_nucl - i);
        }
    }

    std::vector<std::vector<const uint8_t*>> result;
    bool matrix_consumed = false;

	// Add super kmers to the result until all kmers are consumed
    while (!matrix_consumed) {
        std::vector<interleved_t> candidates;

        for (size_t i = 0; i < matrix.size(); i++) {
            size_t position = matrix.size() / 2;
            if (i % 2 == 0)
                position += i / 2;
            else
                position -= 1 + i / 2;

            if (current_kmers[position] < matrix[position].size()) {
                const uint8_t* kmer = matrix[position][current_kmers[position]].get();
                auto it = sk_index.find(kmer);
                if (it != sk_index.end()) {
                    SuperKmerPtr* sk_ptr = it->second;

                	// A super-kmer is available only if all its k-mers are at the front of their respective columns.
                    auto count_it = sk_counts.find(sk_ptr);
                    if (count_it != sk_counts.end() && count_it->second == (*sk_ptr)->size()) {
                        candidates.push_back(current_interleaves[position]);
                    }
                }
            }
        }

        if (candidates.empty()) {
            matrix_consumed = true;
            break;
        }

    	// Select the lexicographically smallest candidate k-mer.
        interleved_t selected = min_interleaved(candidates.begin(), candidates.end());
        size_t col = selected.suf_size;
        const uint8_t* selected_kmer = matrix[col][current_kmers[col]].get();

        auto it = sk_index.find(selected_kmer);
        if (it != sk_index.end()) {
            SuperKmerPtr* selected_sk_ptr = it->second;
            result.push_back(**selected_sk_ptr);

            sk_counts.erase(selected_sk_ptr);

            // Advance pointers for k-mers in the chosen super-kmer and update counts.
            for (const uint8_t* kmer : **selected_sk_ptr) {
                auto col_it = columns.find(kmer);
                if (col_it != columns.end()) {
                    size_t kmer_col = col_it->second;
                    current_kmers[kmer_col]++;

                    if (current_kmers[kmer_col] < matrix[kmer_col].size()) {
                        current_interleaves[kmer_col] = interleaved(
                            matrix[kmer_col][current_kmers[kmer_col]].get(),
                            memory[kmer_col].get(),
                            num_kmer_nucl, num_kmer_nucl - kmer_col);

                        const uint8_t* next_kmer = matrix[kmer_col][current_kmers[kmer_col]].get();
                        auto next_it = sk_index.find(next_kmer);
                        if (next_it != sk_index.end()) {
                            SuperKmerPtr* next_sk_ptr = next_it->second;
                            sk_counts[next_sk_ptr]++;
                        }
                    }
                }
            }
        }

    	// Check if all columns are consumed
        matrix_consumed = true;
        for (size_t i = 0; i < matrix.size(); i++) {
            if (current_kmers[i] < matrix[i].size()) {
                matrix_consumed = false;
                break;
            }
        }
    }

    return result;
}

/**
 * @brief Performs a full sorted-assembly pipeline: sorts columns, finds co-linear pairs, and polishes the final order.
 * @param kmers The KmerMatrix.
 * @return The assembled super-kmer paths.
 */
std::vector<std::vector<const uint8_t*>> Compact::sorted_assembly(KmerMatrix & kmers) {
	// 1 - Sort the matrix column by column.
	this->sort_matrix(kmers);
	std::vector<std::vector<std::pair<uint64_t, uint64_t>>> kmer_pairs;
	for (uint i = 0; i < this->k - this->m; i++) {
		// 2 - Find all possible k-mer overlaps.
		std::vector<std::pair<uint64_t, uint64_t>> candidate_links = this->pair_kmers(kmers[i], kmers[i + 1]);
		// 3 - Filter out k-mer pairs that are not in the optimal co-linear chain.
		std::vector<std::pair<uint64_t, uint64_t>> colinear_links = this->colinear_chaining(candidate_links);
		kmer_pairs.push_back(std::move(colinear_links));
	}
	// 4 - Finalize the sort by ordering swappable super-kmers.
	return this->polish_sort(kmers, kmer_pairs);
}

/**
 * @brief Assembles k-mers using a greedy strategy, linking adjacent columns based on suffix-prefix overlaps.
 * @param kmers The KmerMatrix.
 * @return A list of (predecessor, successor) k-mer pairs.
 */
std::vector<std::pair<const uint8_t*, const uint8_t*>> Compact::greedy_assembly(KmerMatrix & kmers) {
    uint nb_nucl = k - m;
    std::vector<std::pair<const uint8_t*, const uint8_t*>> assembly;

    for (const auto& kmer_ptr : kmers[0]) {
        assembly.emplace_back(nullptr, kmer_ptr.get());
    }

    // Link column i with i+1 in each iteration.
    for (uint i = 0; i < nb_nucl; i++) {
        // Index k-mers in the current column by their suffix for quick lookups.
        std::unordered_map<uint64_t, std::vector<const uint8_t*>> index;

        for (const auto& kmer_ptr : kmers[i]) {
            const uint8_t* kmer = kmer_ptr.get();
            uint64_t val = subseq_to_uint(kmer, nb_nucl, 1, nb_nucl - 1);
            index[val].push_back(kmer);
        }

        // Try to link k-mers from the next column by their prefix.
        for (const auto& kmer_ptr : kmers[i + 1]) {
            const uint8_t* kmer = kmer_ptr.get();
            uint64_t val = subseq_to_uint(kmer, nb_nucl, 0, nb_nucl - 2);

            auto it = index.find(val);
            if (it == index.end()) {
                // If no match, start a new path.
                assembly.emplace_back(nullptr, kmer);
            } else {
                bool chaining_found = false;
                auto& candidates = it->second;

                for (auto candidate_it = candidates.begin(); candidate_it != candidates.end(); ++candidate_it) {
                    const uint8_t* candidate = *candidate_it;

                    if (sequence_compare(
                            kmer, nb_nucl, 0, nb_nucl - 2,
                            candidate, nb_nucl, 1, nb_nucl - 1
                        ) == 0) {
                        // If a match is found, add the pair to the assembly.
                        chaining_found = true;
                        assembly.emplace_back(candidate, kmer);

                        candidates.erase(candidate_it);
                        break;
                    }
                }

                if (!chaining_found) {
                    assembly.emplace_back(nullptr, kmer);
                }
            }
        }
    }

    int assembly_idx = static_cast<int>(assembly.size()) - 1;
    for (auto it = kmers[nb_nucl].rbegin(); it != kmers[nb_nucl].rend(); ++it) {
        const uint8_t* kmer = (*it).get();

        if (assembly_idx < 0 || kmer != assembly[assembly_idx].second) {
            assembly.emplace_back(nullptr, kmer);
        } else {
            assembly_idx--;
        }
    }

    return assembly;
}

/**
 * @brief Converts a list of (predecessor, successor) pairs into contiguous paths.
 * @param to_compact A vector of (predecessor, successor) pairs.
 * @return A vector of k-mer paths.
 */
std::vector<std::vector<const uint8_t*>> Compact::pairs_to_paths(const std::vector<std::pair<const uint8_t*, const uint8_t*>> & to_compact) {
	vector<vector<const uint8_t *>> paths;
	// `path_registry` maps the last k-mer of a path to the path itself.
	unordered_map<const uint8_t *, uint> path_registry;

	for (const pair<const uint8_t *, const uint8_t *> & p : to_compact) {
		// A null predecessor starts a new path.
		if (p.first == nullptr) {
			path_registry[p.second] = paths.size();
			paths.emplace_back(vector<const uint8_t *>());
			paths.back().push_back(p.second);
		}
		else {
			// A non-null predecessor extends an existing path.
			uint vec_idx = path_registry.at(p.first);
			paths[vec_idx].push_back(p.second);
			path_registry.erase(p.first);
			path_registry[p.second] = vec_idx;
		}
	}
	return paths;
}

/**
 * @brief Writes the compacted paths to the output file.
 * Each path is written as a compacted k-mer sequence with its associated minimizer position and data.
 * @param paths The vector of k-mer paths to write.
 * @param sm The Section_Minimizer to write to.
 * @param data_size The size of the data associated with each k-mer.
 */
void Compact::write_paths(const std::vector<std::vector<const uint8_t*>> & paths, Section_Minimizer & sm, const uint data_size) {
    // This function also works with raw pointers and should be fine.
	uint kmer_bytes = (k - m + 3) / 4;
	uint mini_pos_size = (static_cast<uint>(ceil(log2(sm.max + k - m))) + 7) / 8;

	uint max_skmer_bytes = (2 * (k - m) + 3) / 4;
    std::vector<uint8_t> skmer_buffer(max_skmer_bytes + 1, 0);
	uint data_bytes = (k - m + 1) * data_size;
    std::vector<uint8_t> data_buffer(data_bytes, 0);

	for (const vector<const uint8_t *> & path : paths) {
        if (path.empty()) continue;

		// Get the skmer minimizer position
		uint64_t mini_pos = 0;
		const uint8_t * mini_pos_pointer = path[0] + kmer_bytes + data_size;
		for (uint b=0 ; b<mini_pos_size ; b++) {
			mini_pos = (mini_pos << 8) + mini_pos_pointer[b];
		}

		uint skmer_size = k - m - 1 + path.size();
		uint skmer_offset = (4 - (skmer_size % 4)) % 4;

		// Save the first kmer
		memcpy(skmer_buffer.data(), path[0], kmer_bytes);
		leftshift8(skmer_buffer.data(), kmer_bytes, 2 * skmer_offset);
		rightshift8(skmer_buffer.data(), kmer_bytes + 1, 2 * skmer_offset);

        // Save the first data
		memcpy(data_buffer.data(), path[0] + kmer_bytes, data_size);

		// Compact kmer+data one by one
		for (uint kmer_idx = 1 ; kmer_idx<path.size() ; kmer_idx++) {
			const uint8_t * kmer = path[kmer_idx];
			uint compact_nucl_pos = skmer_offset + k - m - 1 + kmer_idx;
			uint compact_byte = compact_nucl_pos / 4;
			uint compact_shift = 3 - (compact_nucl_pos % 4);
			uint8_t last_nucl = kmer[kmer_bytes - 1] & 0b11;
			skmer_buffer[compact_byte] |= last_nucl << (2 * compact_shift);
			memcpy(data_buffer.data() + kmer_idx * data_size, kmer + kmer_bytes, data_size);
		}
		sm.write_compacted_sequence_without_mini(skmer_buffer.data(), skmer_size, mini_pos, data_buffer.data());
	}
}
