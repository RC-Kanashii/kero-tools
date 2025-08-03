#include <string>

#include "CLI/CLI.hpp"
#include "kerotools.hpp"


#ifndef DATARM_H
#define DATARM_H

class DataRm: public KeroTool {
private:
	std::string input_filename;
	std::string output_filename;

public:
	DataRm();
	void cli_prepare(CLI::App * subapp);
	void exec();
};

#endif