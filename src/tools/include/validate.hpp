#include <string>

#include "CLI/CLI.hpp"
#include "kerotools.hpp"
#include "encoding.hpp"


#ifndef VALID_H
#define VALID_H

class Validate: public KeroTool {
private:
	std::string input_filename;
	std::string output_filename;
	bool index_only;
	bool verbose;
	Stringifyer strif;

public:
	Validate();
	void cli_prepare(CLI::App * subapp);
	void exec();

	bool is_valid_r_section(Kero_file & infile);
	bool is_valid_M_section(Kero_file & infile);
	bool is_valid_i_section(Kero_file & infile);
	bool is_valid_h_section(Kero_file & infile);
};

#endif