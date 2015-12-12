#include <iostream>
#include <fstream>
#include "String.h"
#include "Fasta.h"

bool Fasta::Load(const std::string& filename, bool verbose)
{
	std::ifstream file(filename, std::ios::in);
	if (!file.is_open()) {
		std::cerr << "Error: Can not open file '" << filename << "'!" << std::endl;
		return false;
	}

	if (verbose) {
		std::cerr << "Loading fasta '" << filename << "'" << std::endl;
	}

	std::string chrom;
	std::string line;
	while (std::getline(file, line)) {
		if (line.empty()) continue;
		if (line[0] == '>') {
			chrom = TrimLeft(line.substr(1));
			std::string::size_type pos = chrom.find_first_of(" \t");
			if (pos != std::string::npos) {
				chrom = chrom.substr(0, pos);
			}
			if (verbose) {
				std::cerr << "  loading '" << chrom << "'\r" << std::flush;
			}
		} else {
			seq_[chrom] += Trim(line);
		}
	}
	file.close();

	if (verbose) {
		std::cerr << "Total " << seq_.size() << " sequence(s) loaded" << std::endl;
	}
	return true;
}

bool Fasta::Has(const std::string& chrom) const
{
	return (seq_.find(chrom) != seq_.end());
}

std::string Fasta::GetSeq(const std::string& chrom, size_t pos, size_t size) const
{
	std::string res;
	auto it = seq_.find(chrom);
	if (it != seq_.end()) {
		std::string s = it->second.substr(pos, size);
		for (size_t i = 0; i < s.size(); ++i) {
			res += std::toupper(s[i]);
		}
	}
	return res;
}
