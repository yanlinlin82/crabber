#include <string>
#include <iostream>
#include <fstream>
#include "Fasta.h"
#include "LineSplit.h"
#include "RegionGet.h"

const int LINE_WIDTH = 60;

static bool Process(const std::string& filename, const Fasta& fa)
{
	std::ifstream file(filename, std::ios::in);
	if (!file.is_open()) {
		std::cerr << "Error: Can not open file '" << filename << "'!" << std::endl;
		return false;
	}

	size_t lineNo = 0;
	std::string line;
	while (std::getline(file, line)) {
		++lineNo;
		if (line.empty() || line[0] == '#') continue;

		LineSplit sp;
		sp.Split(line, '\t');

		try {
			std::string chrom = sp.GetField(0);
			int start = stoi(sp.GetField(1));
			int end = stoi(sp.GetField(2));

			if (start < 0) {
				start = 0;
			}
			if (start >= end) {
				std::cerr << "Skip invalid region: " << chrom << ":" << start + 1 << "-" << end << std::endl;
				continue;
			}
			size_t len = fa.GetLength(chrom);
			if (len == 0) {
				std::cerr << "Skip non-existed sequence: '" << chrom << "'!" << std::endl;
				continue;
			}
			if (start >= static_cast<int>(len)) {
				std::cerr << "Skip non-existed region: " << chrom << ":" << start + 1 << "-" << end << std::endl;
				continue;
			}

			std::cout << ">" << chrom << ":" << start + 1 << "-" << end << std::endl;

			if (end > static_cast<int>(len)) {
				end = len;
			}
			for (int i = start; i < end; i += LINE_WIDTH) {
				int size = LINE_WIDTH;
				if (i + size > end) {
					size = end - i;
				}
				std::cout << fa.GetSeq(chrom, i + 1, size) << std::endl;
			}
		} catch (const std::exception& e) {
			std::cerr << "Unexpected error in line " << lineNo << " of file '" << filename << "'! " << e.what() << std::endl;
			file.close();
			return false;
		}
	}

	file.close();
	return true;
}
static void PrintUsage()
{
	std::cout << "\n"
		"Usage:  crabber region-get <ref.fa> <region.bed>\n"
		"\n"
		"Input:\n"
		"   <ref.fa>        reference genome in FASTA format\n"
		"   <region.bed>    target region to extract sequences\n"
		<< std::endl;
}

int RegionGet_main(int argc, char* const argv[])
{
	std::vector<std::string> args(argv, argv + argc);
	if (args.size() < 3) {
		PrintUsage();
		return 1;
	}

	Fasta fa;
	if (!fa.Load(args[1])) {
		return 1;
	}

	if (!Process(args[2], fa)) {
		return 1;
	}
	return 0;
}
