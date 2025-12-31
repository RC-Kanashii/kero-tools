#include <iostream>
#include <vector>
#include <sys/resource.h>

#include "kerotools.hpp"

#include "bucket.hpp"
#include "compact.hpp"
#include "datarm.hpp"
#include "disjoin.hpp"
#include "index.hpp"
#include "instr.hpp"
#include "merge.hpp"
#include "outstr.hpp"
#include "shuffle.hpp"
#include "sort.hpp"
#include "split.hpp"
#include "translate.hpp"
#include "validate.hpp"
#include "query.hpp"


using namespace std;


KeroTool * parse_args(int argc, char** argv, vector<KeroTool *> tools) {
	// Main command
	CLI::App app{"kero-tools is a software for kero file manipulations. For more details on kero format, please refer to https://github.com/Kmer-File-Format/kero-reference"};
	app.require_subcommand(1);
	CLI::Option * help =	app.get_help_ptr();

	// Subcommands prepare
	for (KeroTool * tool : tools) {
		tool->cli_prepare(&app);
	}

	// Parsing and return status if wrong
	try {
    app.parse(argc, argv);
	} catch (const CLI::ParseError &e) {
    auto val = app.exit(e);
		if (val != 0) {
			exit(val);
		}
	}

	// Help detection
	if (!help->empty()) {
		return nullptr;
	}

	// Read the command line return
	for (KeroTool * tool : tools) {
		if (tool->subapp->parsed()) {
			// Help on tool triggered
			if (not tool->subapp->get_help_ptr()->empty()) {
				return nullptr;
			} else {
				return tool;
			}
		}
	}

	return nullptr;
}

int main(int argc, char** argv) {
	// --- System calls for optimization ---
	// Remove interactive synchronization for speedup I/O
	// ios_base::sync_with_stdio(false);
	// Raise the number of simultaneous file descriptors to maximum
	struct rlimit nb_file_descriptors;
	getrlimit(RLIMIT_NOFILE, &nb_file_descriptors);
	nb_file_descriptors.rlim_cur = nb_file_descriptors.rlim_max;
	setrlimit(RLIMIT_NOFILE, &nb_file_descriptors);


	// --- Prepare tools ---
	vector<KeroTool *> tools;
	tools.push_back(new Bucket());
	tools.push_back(new Compact());
	tools.push_back(new DataRm());
	tools.push_back(new Disjoin());
	tools.push_back(new Index());
	tools.push_back(new Instr());
	tools.push_back(new Merge());
	tools.push_back(new Outstr());
	tools.push_back(new Shuffle());
	tools.push_back(new Sort());
	tools.push_back(new Split());
	tools.push_back(new Translate());
	tools.push_back(new Validate());
	tools.push_back(new Query());

	try {
		// Get the one selected
		KeroTool * tool = parse_args(argc, argv, tools);

		if (tool != nullptr)
			tool->exec();
	} catch (const std::runtime_error& e) {
		std::cerr << "Runtime error: " << e.what() << std::endl;
		return 1;
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	} catch (...) {
		std::cerr << "Unknown error occurred" << std::endl;
		return 1;
	}

	for (KeroTool * tool : tools)
		delete tool;

	return 0;
}