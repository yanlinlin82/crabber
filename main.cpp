#include <iostream>
#include <string>
#include "DepthStat.h"
#include "Annotate.h"
#include "RegionGet.h"
#include "RegionCount.h"
#include "version.h"

static void PrintUsage(const char* progname)
{
	std::cout << "\n"
		"Program: Crabber catches crabs :P\n"
		"Version: " << VERSION << "\n"
		"\n"
		"Usage: " << progname << " <command> [options]\n"
		"\n"
		"Commands:\n"
		"    depth-stat     stat coverage depth\n"
		"    region-get     extract sequences in regions\n"
		"    region-count   count bases in regions\n"
		"    annotate       annotate genetic mutations\n"
		<< std::endl;
}

int main(int argc, char* const argv[])
{
	if (argc < 2) {
		PrintUsage(argv[0]);
		return 1;
	}

	std::string cmd(argv[1]);
	if (cmd == "depth-stat") {
		return DepthStat_main(argc - 1, argv + 1);
	} else if (cmd == "region-get") {
		return RegionGet_main(argc - 1, argv + 1);
	} else if (cmd == "region-count") {
		return RegionCount_main(argc - 1, argv + 1);
	} else if (cmd == "annotate") {
		return Annotate_main(argc - 1, argv + 1);
	} else {
		std::cerr << "Error: Unknown command '" << argv[1] << "'!\n" << std::endl;
		PrintUsage(argv[0]);
		return 1;
	}
}
