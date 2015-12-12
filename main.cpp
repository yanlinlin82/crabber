#include <iostream>
#include <string>
#include "annotate.h"
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
		"    annotate     annotate genetic mutations\n"
		<< std::endl;
}

int main(int argc, char* const argv[])
{
	if (argc < 2) {
		PrintUsage(argv[0]);
		return 1;
	}

	std::string cmd(argv[1]);
	if (cmd == "annotate") {
		return main_annotate(argc - 1, argv + 1);
	} else {
		std::cerr << "Error: Unknown command '" << argv[1] << "'!\n" << std::endl;
		PrintUsage(argv[0]);
		return 1;
	}
}
