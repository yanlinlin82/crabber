#include <iostream>
#include <fstream>
#include <map>
#include "LineSplit.h"
#include "DepthStat.h"

static bool Process(const std::string& filename)
{
	std::ifstream file(filename, std::ios::in);
	if (!file.is_open()) {
		std::cerr << "Error: Can not open file '" << filename << "'!" << std::endl;
		return false;
	}

	std::map<int, int> depthStat;

	std::string line;
	int lineNo = 0;
	long long totalBases = 0;
	while (std::getline(file, line)) {
		++lineNo;

		LineSplit sp;
		sp.Split(line, '\t', 5);

		try {
			int depth = stoi(sp.GetField(3));
			++depthStat[depth];
			totalBases += depth;

		} catch (const std::exception& e) {
			std::cerr << "Unexpected error in line " << lineNo << " of file '" << filename << "'! " << e.what() << std::endl;
			file.close();
			return false;
		}
	}
	file.close();

	std::cout << "depth\tcount\tratio" << std::endl;
	long long bases = 0;
	for (auto it = depthStat.begin(); it != depthStat.end(); ++it) {
		bases += static_cast<long long>(it->first) * it->second;
		double ratio = static_cast<double>(bases) / totalBases;
		std::cout << it->first << '\t' << it->second << '\t' << ratio << std::endl;
	}
	return true;
}

int DepthStat_main(int argc, char* const argv[])
{
	if (argc < 2) {
		std::cout << "Usage: crabber depth-stat <x.mpileup>" << std::endl;
		return 1;
	}

	if (!Process(argv[1])) {
		return 1;
	}
	return 0;
}
