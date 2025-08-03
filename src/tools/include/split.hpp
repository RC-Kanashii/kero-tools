// #include <experimental/filesystem>
#include <string>
#include <iostream>

#include "CLI/CLI.hpp"
#include "kerotools.hpp"


#ifndef SPLIT_H
#define SPLIT_H

class Split: public KeroTool {
private:
	std::string input_filename;
	std::string output_dirname;

public:
	Split();
	void cli_prepare(CLI::App * subapp);
	void exec();
};


#endif