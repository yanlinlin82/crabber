#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cassert>
#include "annotate.h"

std::string trim_left(const std::string& str, const std::string& drop = " \t")
{
	std::string::size_type pos = str.find_first_not_of(drop);
	if (pos == std::string::npos) {
		return "";
	}
	return str.substr(pos);
}

std::string trim_right(const std::string& str, const std::string& drop = " \t")
{
	std::string::size_type pos = str.find_last_not_of(drop);
	if (pos == std::string::npos) {
		return "";
	}
	return str.substr(0, pos + 1);
}

std::string trim(const std::string& str, const std::string& drop = " \t")
{
	return trim_left(trim_right(str, drop), drop);
}

static bool LoadFasta(const std::string& filename, std::map<std::string, std::string>& ref)
{
	std::ifstream file(filename, std::ios::in);
	if (!file.is_open()) {
		std::cerr << "Error: Can not open file '" << filename << "'!" << std::endl;
		return false;
	}

	//std::cerr << "Loading fasta '" << filename << "'" << std::endl;
	std::string chrom;
	std::string line;
	while (std::getline(file, line)) {
		if (line.empty()) continue;
		if (line[0] == '>') {
			chrom = trim_left(line.substr(1));
			std::string::size_type pos = chrom.find_first_of(" \t");
			if (pos != std::string::npos) {
				chrom = chrom.substr(0, pos);
			}
			//std::cerr << "  loading '" << chrom << "'\r" << std::flush;
		} else {
			ref[chrom] += trim(line);
		}
	}
	//std::cerr << "Total " << ref.size() << " sequence(s) loaded" << std::endl;
	file.close();
#if 0
	for (auto it = ref.begin(); it != ref.end(); ++it) {
		std::cout << it->first << '\t' << it->second.size() << std::endl;
	}
#endif
	return true;
}

void Split(const char* s, size_t* pos, size_t count, char sep)
{
	*pos++ = 0;
	--count;

	size_t i;
	for (i = 0; s[i] && count > 0; ++i) {
		if (s[i] == sep) {
			*pos++ = i + 1;
			--count;
		}
	}
	if (count > 0) {
		*pos = i;
	}
}

std::string GetField(const std::string& line, size_t index, size_t* pos)
{
	return line.substr(pos[index], pos[index + 1] - pos[index] - 1);
}

class Transcript
{
public:
	Transcript(): txStart_(0), txEnd_(0) { }

	Transcript(const std::string& name, const std::string& name2, int txStart, int txEnd):
		name_(name), name2_(name2), txStart_(txStart), txEnd_(txEnd)
	{
	}

	bool SetExons(int exonCount, const std::string& exonStarts, const std::string& exonEnds);
public:
	std::string name_;
	std::string name2_;
	std::string strand_;
	int txStart_;
	int txEnd_;
	int cdsStart_;
	int cdsEnd_;
	int exonCount_;
	std::string exonStarts_;
	std::string exonEnds_;
	std::vector<std::pair<int, int>> exons_;
};

bool Transcript::SetExons(int exonCount, const std::string& exonStarts, const std::string& exonEnds)
{
	std::vector<size_t> pos, pos2;
	pos.resize(exonCount + 1);
	pos2.resize(exonCount + 1);

	Split(exonStarts.c_str(), &pos[0], exonCount + 1, ',');
	Split(exonEnds.c_str(), &pos2[0], exonCount + 1, ',');
	for (int i = 0; i < exonCount; ++i) {
		exons_.push_back(std::make_pair(
					stoi(GetField(exonStarts.c_str(), i, &pos[0])) - txStart_,
					stoi(GetField(exonEnds.c_str(), i, &pos2[0])) - txStart_));
	}
	return true;
}

bool LoadRefGene(const std::string& filename, std::map<std::string, std::vector<Transcript>>& data)
{
	std::fstream file(filename.c_str(), std::ios::in);
	if (!file.is_open()) {
		std::cerr << "Can not open file '" << filename << "'!" << std::endl;
		return false;
	}

	//std::cerr << "Loading refGene '" << filename << "'" << std::endl;
	time_t t0 = time(NULL);
	size_t lineNo = 0;
	int count = 0;
	std::string line;
	while (std::getline(file, line)) {
		++lineNo;

		time_t t = time(NULL);
		if (t != t0) {
			//std::cerr << "  " << count << " record(s) loaded\r" << std::flush;
			t0 = t;
		}

		if (line.empty() || line[0] == '#') continue;

		size_t pos[20] = { };
		Split(line.c_str(), pos, sizeof(pos) / sizeof(pos[0]), '\t');

		try {
			std::string cdsStartStat = GetField(line, 13, pos);
			std::string cdsEndStat = GetField(line, 14, pos);
			if (cdsStartStat != "cmpl" || cdsEndStat != "cmpl") continue;

			std::string name = GetField(line, 1, pos);
			std::string name2 = GetField(line, 12, pos);
			std::string chrom = GetField(line, 2, pos);
			std::string strand = GetField(line, 3, pos);
			int txStart = stoi(GetField(line, 4, pos));
			int txEnd = stoi(GetField(line, 5, pos));
			int cdsStart = stoi(GetField(line, 6, pos));
			int cdsEnd = stoi(GetField(line, 7, pos));
			int exonCount = stoi(GetField(line, 8, pos));
			std::string exonStarts = GetField(line, 9, pos);
			std::string exonEnds = GetField(line, 10, pos);

			Transcript item(name, name2, txStart, txEnd);
			item.strand_ = strand;
			item.cdsStart_ = cdsStart;
			item.cdsEnd_ = cdsEnd;
			if (!item.SetExons(exonCount, exonStarts, exonEnds)) {
				std::cerr << "Invalid exon info at line " << lineNo << std::endl;
			}
			item.exonCount_ = exonCount;
			item.exonStarts_ = exonStarts;
			item.exonEnds_ = exonEnds;
			data[chrom].push_back(item);

			++count;
		} catch (const std::exception& e) {
			std::cerr << "Unexpected error in line " << lineNo << " of file '" << filename << "'! " << e.what() << std::endl;
			file.close();
			return false;
		}
	}
	//std::cerr << "Total " << count << " record(s) loaded" << std::endl;
	file.close();
	return true;
}

int GetTxPos(const std::vector<std::pair<int, int>>& exons, int pos)
{
	int count = 0;
	for (size_t i = 0; i < exons.size(); ++i) {
		int start = exons[i].first;
		int end = exons[i].second;
		if (pos < start) {
			throw std::runtime_error("bad pos (" + std::to_string(pos) + " < " + std::to_string(start) + ")");
		}
		if (pos < end) {
			return count + (pos - start);
		}
		count += end - start;
	}
	throw std::runtime_error("bad pos (" + std::to_string(pos) + " >= " + std::to_string(exons.back().second) + ")");
}

int GetTxPosRev(const std::vector<std::pair<int, int>>& exons, int pos)
{
	int count = 0;
	for (size_t i = exons.size(); i > 0; --i) {
		int start = exons[i - 1].first;
		int end = exons[i - 1].second;
		if (pos >= end) {
			throw std::runtime_error("bad pos (" + std::to_string(pos) + " >= " + std::to_string(end) + ")");
		}
		if (pos >= start) {
			return count + (end - pos);
		}
		count += end - start;
	}
	throw std::runtime_error("bad pos (" + std::to_string(pos) + " < " + std::to_string(exons.front().first) + ")");
}

bool Convert(const std::string& chrom, const Transcript& trans, int pos, const std::string& ref, const std::string alt, const std::string& seq)
{
#if 0
	static bool firstTime = true;
	if (firstTime) {
		//firstTime = false;
		std::cout << "Convert: " << trans.name2_ << '\t' << pos << '\t' << trans.strand_ << '\t' << seq.size() << std::endl;
		std::cout << "   cds: " << trans.cdsStart_ << ", " << trans.cdsEnd_ << std::endl;
		std::cout << "   cds: " << trans.cdsStart_ - trans.txStart_ << ", " << trans.cdsEnd_ - trans.txStart_ << std::endl;
		for (size_t i = 0; i < trans.exons_.size(); ++i) {
			std::cout << "   (" << trans.exons_[i].first + trans.txStart_ << '\t' << trans.exons_[i].second + trans.txStart_ << ")  "
				<< trans.exons_[i].first << '\t' << trans.exons_[i].second << std::endl;
		}
	}
#endif
	std::string res = ".";
	std::string type = ".";
	std::string type2 = ".";

	const auto& exons = trans.exons_;
	int cdsStart = trans.cdsStart_ - trans.txStart_;
	int cdsEnd = trans.cdsEnd_ - trans.txStart_;

	if (trans.strand_ == "+") {
		int cdsStartTxPos = GetTxPos(exons, cdsStart);
		int cdsEndTxPos = GetTxPos(exons, cdsEnd - 1) + 1;

		for (size_t i = 0; i < exons.size(); ++i) {
			int start = exons[i].first;
			int end = exons[i].second;
			if (pos < cdsStart) { // 5'-UTR
				if (pos < start) { // intron
					if (i == 0 || (pos - exons[i - 1].second) > (start - pos)) {
						res = "c.-" + std::to_string(cdsStartTxPos - GetTxPos(exons, start)) + "-" + std::to_string(start - pos);
					} else {
						res = "c.-" + std::to_string(cdsStartTxPos - GetTxPos(exons, exons[i - 1].second - 1)) + "+" + std::to_string(pos - (exons[i - 1].second - 1));
					}
					type = "Intron(" + std::to_string(i) + "/" + std::to_string(exons.size() - 1) + ")";
					break;
				} else if (pos < end) { // exon
					res = "c.-" + std::to_string(cdsStartTxPos - GetTxPos(exons, pos));
					type = "Exon(" + std::to_string(i + 1) + "/" + std::to_string(exons.size()) + ")";
					type2 = "5'-UTR";
					break;
				}
			} else if (pos >= cdsEnd) { // 3'-UTR
				if (pos < start) { // intro
					if (i == 0 || (pos - exons[i - 1].second) > (start - pos)) {
						res = "c.*" + std::to_string(GetTxPos(exons, start) - cdsEndTxPos + 1) + "-" + std::to_string(start - pos);
					} else {
						res = "c.*" + std::to_string(GetTxPos(exons, exons[i - 1].second - 1) - cdsEndTxPos + 1) + "+" + std::to_string(pos - (exons[i - 1].second - 1));
					}
					type = "Intron(" + std::to_string(i) + "/" + std::to_string(exons.size() - 1) + ")";
					break;
				} else if (pos < end) { // exon
					res = "c.*" + std::to_string(GetTxPos(exons, pos) - cdsEndTxPos + 1);
					type = "Exon(" + std::to_string(i + 1) + "/" + std::to_string(exons.size()) + ")";
					type2 = "3'-UTR";
					break;
				}
			} else { // CDS
				if (pos < start) { // intron
					if (i == 0 || (pos - exons[i - 1].second) > (start - pos)) {
						res = "c." + std::to_string(GetTxPos(exons, start) - cdsStartTxPos + 1) + "-" + std::to_string(start - pos);
					} else {
						res = "c." + std::to_string(GetTxPos(exons, exons[i - 1].second - 1) - cdsStartTxPos + 1) + "+" + std::to_string(pos - (exons[i - 1].second - 1));
					}
					type = "Intron(" + std::to_string(i) + "/" + std::to_string(exons.size() - 1) + ")";
					break;
				} else if (pos < end) { // exon
					res = "c." + std::to_string(GetTxPos(exons, pos) - cdsStartTxPos + 1);
					type = "Exon(" + std::to_string(i + 1) + "/" + std::to_string(exons.size()) + ")";
					type2 = "CDS";
					break;
				}
			}
		}
	} else {
		assert(trans.strand_ == "-");
		int cdsStartTxPos = GetTxPosRev(exons, cdsStart);
		int cdsEndTxPos = GetTxPosRev(exons, cdsEnd - 1) + 1;

		for (size_t i = exons.size(); i > 0; --i) {
			int start = exons[i - 1].first;
			int end = exons[i - 1].second;
			if (pos >= cdsEnd) { // 5'-UTR
				if (pos >= end) { // intron
					if (i == exons.size() || (pos - end) < (exons[i].first - pos)) {
						res = "c.-" + std::to_string(cdsEndTxPos - GetTxPosRev(exons, end - 1)) + "-" + std::to_string(pos - end);
					} else {
						res = "c.-" + std::to_string(cdsEndTxPos - GetTxPosRev(exons, exons[i].first)) + "+" + std::to_string(exons[i].first - pos);
					}
					type = "Intron(" + std::to_string(exons.size() - i) + "/" + std::to_string(exons.size() - 1) + ")";
					break;
				} else if (pos >= start) { // exon
					res = "c.-" + std::to_string(cdsEndTxPos - GetTxPosRev(exons, pos) + 1);
					type = "Exon(" + std::to_string(exons.size() - i + 1) + "/" + std::to_string(exons.size()) + ")";
					type2 = "5'-UTR";
					break;
				}
			} else if (pos < cdsStart) { // 3'-UTR
				if (pos >= end) { // intron
					if (i == exons.size() || (pos - end) < (exons[i].first - pos)) {
						res = "c.*" + std::to_string(cdsStartTxPos - GetTxPosRev(exons, end - 1)) + "-" + std::to_string(pos - end);
					} else {
						res = "c.*" + std::to_string(cdsStartTxPos - GetTxPosRev(exons, exons[i].first)) + "+" + std::to_string(exons[i].first - pos);
					}
					type = "Intron(" + std::to_string(exons.size() - i) + "/" + std::to_string(exons.size() - 1) + ")";
					break;
				} else if (pos >= start) { // exon
					res = "c.*" + std::to_string(cdsStartTxPos - GetTxPosRev(exons, pos) + 1);
					type = "Exon(" + std::to_string(exons.size() - i + 1) + "/" + std::to_string(exons.size()) + ")";
					type2 = "3'-UTR";
					break;
				}
			} else { // CDS
				if (pos >= end) { // intron
					if (i == exons.size() || (pos - end) < (exons[i].first - pos)) {
						res = "c." + std::to_string(GetTxPosRev(exons, end - 1) - cdsEndTxPos) + "-" + std::to_string(pos - end + 1);
					} else {
						res = "c." + std::to_string(GetTxPosRev(exons, exons[i].first) - cdsEndTxPos) + "+" + std::to_string(exons[i].first - pos);
					}
					type = "Intron(" + std::to_string(exons.size() - i) + "/" + std::to_string(exons.size() - 1) + ")";
					break;
				} else if (pos >= start) { // exon
					res = "c." + std::to_string(GetTxPosRev(exons, pos - 1) - cdsEndTxPos + 1);
					type = "Exon(" + std::to_string(exons.size() - i + 1) + "/" + std::to_string(exons.size()) + ")";
					type2 = "CDS";
					break;
				}
			}
		}
	}

	//assert(ref == seq.substr(pos, 1));
	res += std::toupper(seq.substr(pos, 1)[0]);
	res += ">" + alt;

	std::cout << chrom << "\t" << trans.txStart_ + pos + 1 << "\t" << trans.strand_ << "\t" << trans.name_ << "\t" << trans.name2_ << "\t" << res << "\t" << type << "\t" << type2 << std::endl;
	return true;
}

bool ProcessItem(const std::string& chrom, int pos, const std::string& alleleRef, const std::string& alleleAlt,
		const std::map<std::string, std::vector<Transcript>>& data,
		const std::map<std::string, std::string>& ref)
{
	//std::cout << chrom << '\t' << pos << '\t' << alleleRef << '\t' << alleleAlt << std::endl;
	//std::cout << "Process: " << chrom << ':' << pos << ',' << alleleRef << '>' << alleleAlt << std::endl;

	auto it = data.find(chrom);
	if (it != data.end()) {
		auto trans = it->second;
		bool found = false;
		for (size_t i = 0; i < trans.size(); ++i) {
			if (pos >= trans[i].txStart_ && pos < trans[i].txEnd_) {
				//std::cerr << "Found in " << trans[i].name_ << " (" << chrom << ':' << trans[i].txStart_ << "-" << trans[i].txEnd_ << ")" << std::endl;

				auto it2 = ref.find(chrom);
				if (it2 == ref.end()) {
					std::cerr << "Can not found sequence '" << chrom << "' in ref fasta" << std::endl;
					return false;
				}

				Convert(chrom, trans[i], pos - trans[i].txStart_, alleleRef, alleleAlt, it2->second.substr(trans[i].txStart_, trans[i].txEnd_));
				found = true;
			}
		}
		if (!found) {
			std::cout << chrom << "\t" << pos + 1 << "\t.\t.\t.\t.\tIntergenic\t." << std::endl;
		}
	}
	return true;
}

bool Process(const std::string& filename,
		const std::map<std::string, std::vector<Transcript>>& data,
		const std::map<std::string, std::string>& ref)
{
	std::ifstream file(filename, std::ios::in);
	if (!file.is_open()) {
		std::cerr << "Error: Can not open file '" << filename << "'!" << std::endl;
		return false;
	}

	std::cout << "chrom\tpos\tstrand\tname\tname2\tmutate\tsegment\ttype" << std::endl;

	size_t lineNo = 0;
	std::string line;
	while (std::getline(file, line)) {
		++lineNo;
		if (line.empty() || line[0] == '#') continue;

		size_t pos[20] = { };
		Split(line.c_str(), pos, sizeof(pos) / sizeof(pos[0]), '\t');

		try {
			std::string chrom = GetField(line, 0, pos);
			int genomePos = stoi(GetField(line, 1, pos));
			std::string alleleRef = GetField(line, 3, pos);
			std::string alleleAlt = GetField(line, 4, pos);

			ProcessItem(chrom, genomePos - 1, alleleRef, alleleAlt, data, ref);
		} catch (const std::exception& e) {
			std::cerr << "Unexpected error in line " << lineNo << " of file '" << filename << "'! " << e.what() << std::endl;
			file.close();
			return false;
		}
	}

	file.close();
	return true;
}

int main_annotate(int argc, char* const argv[])
{
	std::vector<std::string> args(argv, argv + argc);
	if (args.size() < 4) {
		std::cout << "\n"
			"Usage:  crabber annotate <x.vcf> <refGene.tsv> <ref.fa>\n"
			"\n"
			"Input:\n"
			"   <x.vcf>         input SNV list in VCF format\n"
			"   <refGene.tsv>   track data downloaded from UCSC table browser\n"
			"   <ref.fa>        reference genome in FASTA format\n"
			<< std::endl;
		return 1;
	}

	//std::cout << "fasta: " << args[3] << std::endl;
	//std::cout << "refGene: " << args[2] << std::endl;
	//std::cout << "vcf: " << args[1] << std::endl;

	std::map<std::string, std::string> ref;
	if (!LoadFasta(args[3], ref)) {
		return 1;
	}

	std::map<std::string, std::vector<Transcript>> data;
	if (!LoadRefGene(args[2], data)) {
		return 1;
	}

	if (!Process(args[1], data, ref)) {
		return 1;
	}
	return 0;
}
