#include <vector>
#include <string>
#include <cstring>
#include <fstream>

#include "outstr.hpp"
#include "encoding.hpp"


using namespace std;


Outstr::Outstr() {
	input_filename = "";
	output_filename = "";
	revcomp = false;
}

void Outstr::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("outstr", "Output kmer sequences as string. One kmer per line is printed");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "The file to print");
	input_option->required();
	input_option->check(CLI::ExistingFile);
	subapp->add_option("-o, --outfile", output_filename, "Output text file (if not specified, output to stdout)");
	subapp->add_flag("-c, --reverse-complement", revcomp, "Print the minimal value between a kmer and its reverse complement");
}


string format_data(uint8_t * data, size_t data_size) {
	if (data_size == 0)
		return "";
	else if (data_size < 8) {
		uint val = data[0];
		for (uint i=1 ; i<data_size ; i++) {
			val <<= 8;
			val += data[i];
		}
		return std::to_string(val);
	} else {
		string val = "";
		for (uint i=0 ; i<data_size ; i++) {
			val += "[" + std::to_string((uint)data[i]) + "]";
		}
		return val;
	}
}



bool inf_eq(uint8_t * seq1, uint8_t * seq2, uint64_t size) {
	uint64_t nb_bytes = (size + 3) / 4;

	// Test first byte
	uint16_t mask = (1 << ((size % 4) * 2)) - 1;
	uint8_t byte_seq1 = seq1[0] & mask;
	uint8_t byte_seq2 = seq2[0] & mask;
	if (byte_seq1 != byte_seq2)
		return byte_seq1 < byte_seq2;

	// Test all remaining bytes
	for (uint idx=1 ; idx<nb_bytes ; idx++)
		if (seq1[idx] != seq2[idx])
			return seq1[idx] < seq2[idx];

	// Equals
	return true;
}


void Outstr::exec() {
	// Read the encoding and prepare the translator
	Kero_reader reader = Kero_reader(input_filename);
	Stringifyer strif(reader.get_encoding());

	// Prepare revcomp
	RevComp rc(reader.get_encoding());
	uint8_t * rc_copy = new uint8_t[1];
	uint64_t k = 0;

	// Prepare sequence and data buffers
	uint8_t * nucleotides = nullptr;
	uint8_t * data = nullptr;

	// Open output file if specified, otherwise use stdout
	ofstream outfile;
	ostream* out = &cout;
	if (!output_filename.empty()) {
		outfile.open(output_filename);
		if (!outfile.is_open()) {
			cerr << "Error: Could not open output file: " << output_filename << endl;
			return;
		}
		out = &outfile;
	}

	while (reader.next_kmer(nucleotides, data)) {

		if (not revcomp) {
			*out << strif.translate(nucleotides, reader.k) << " ";
			*out << format_data(data, reader.data_size) << '\n';
		} else {
			// Change the size of rev comp datastruct if k changes
			if (reader.k != k) {
				k = reader.k;
				delete[] rc_copy;
				rc_copy = new uint8_t[(k+3) / 4];
			}

			// Get the reverse complement
			memcpy(rc_copy, nucleotides, (k+3)/4);
			rc.rev_comp(rc_copy, k);

			if (inf_eq(nucleotides, rc_copy, k)) {
				*out << strif.translate(nucleotides, k) << " ";
				*out << format_data(data, reader.data_size) << '\n';
			} else {
				*out << strif.translate(rc_copy, k) << " ";
				*out << format_data(data, reader.data_size) << '\n';
			}
		}
	}

	// Close output file if it was opened
	if (outfile.is_open()) {
		outfile.close();
	}
}
