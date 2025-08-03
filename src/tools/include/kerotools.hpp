#pragma once

#include "CLI/CLI.hpp"
#include "kero-api/kero_io.hpp"

#ifndef KERO_TOOLS
#define KERO_TOOLS


class KeroTool {
public:
	CLI::App * subapp = nullptr;
	virtual void cli_prepare(CLI::App * subapp) = 0;
	virtual void exec() = 0;
	virtual ~KeroTool() {};
};


#endif
