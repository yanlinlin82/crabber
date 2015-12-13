#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cassert>
#include "annotate.h"
#include "String.h"
#include "Fasta.h"
#include "LineSplit.h"
#include "Transcript.h"

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

		LineSplit sp;
		sp.Split(line, '\t');

		try {
			std::string cdsStartStat = sp.GetField(13);
			std::string cdsEndStat = sp.GetField(14);
			if (cdsStartStat != "cmpl" || cdsEndStat != "cmpl") continue;

			std::string name = sp.GetField(1);
			std::string name2 = sp.GetField(12);
			std::string chrom = sp.GetField(2);
			std::string strand = sp.GetField(3);
			int txStart = stoi(sp.GetField(4));
			int txEnd = stoi(sp.GetField(5));
			int cdsStart = stoi(sp.GetField(6));
			int cdsEnd = stoi(sp.GetField(7));
			int exonCount = stoi(sp.GetField(8));
			std::string exonStarts = sp.GetField(9);
			std::string exonEnds = sp.GetField(10);

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
			throw std::runtime_error("GetTxPos - bad pos (" + std::to_string(pos) + " < " + std::to_string(start) + ")");
		}
		if (pos < end) {
			return count + (pos - start);
		}
		count += end - start;
	}
	throw std::runtime_error("GetTxPos - bad pos (" + std::to_string(pos) + " >= " + std::to_string(exons.back().second) + ")");
}

int GetTxPosRev(const std::vector<std::pair<int, int>>& exons, int pos)
{
	int count = 0;
	for (size_t i = exons.size(); i > 0; --i) {
		int start = exons[i - 1].first;
		int end = exons[i - 1].second;
		if (pos >= end) {
			throw std::runtime_error("GetTxPosRev - bad pos (" + std::to_string(pos) + " >= " + std::to_string(end) + ")");
		}
		if (pos >= start) {
			return count + (end - pos);
		}
		count += end - start;
	}
	throw std::runtime_error("GetTxPosRev - bad pos (" + std::to_string(pos) + " < " + std::to_string(exons.front().first) + ")");
}

const char* CODON_TABLE[4][4][4] = {
	{
		{ "F", "F", "L", "L" }, { "S", "S", "S", "S" },
		{ "Y", "Y", "*", "*" }, { "C", "C", "*", "W" },
	},
	{
		{ "L", "L", "L", "L" }, { "P", "P", "P", "P" },
		{ "H", "H", "Q", "Q" }, { "R", "R", "R", "R" },
	},
	{
		{ "I", "I", "I", "M" }, { "T", "T", "T", "T" },
		{ "N", "N", "K", "K" }, { "S", "S", "R", "R" },
	},
	{
		{ "V", "V", "V", "V" }, { "A", "A", "A", "A" },
		{ "D", "D", "E", "E" }, { "G", "G", "G", "G" },
	},
};

const char* CODON_TABLE_3[4][4][4] = {
	{
		{ "Phe", "Phe", "Leu", "Leu" }, { "Ser", "Ser", "Ser", "Ser" },
		{ "Tyr", "Tyr", "*", "*" }, { "Cys", "Cys", "*", "Trp" }
	},
	{
		{ "Leu", "Leu", "Leu", "Leu" }, { "Pro", "Pro", "Pro", "Pro" },
		{ "His", "His", "Gln", "Gln" }, { "Arg", "Arg", "Arg", "Arg" },
	},
	{
		{ "Ile", "Ile", "Ile", "Met" }, { "Thr", "Thr", "Thr", "Thr" },
		{ "Asn", "Asn", "Lys", "Lys" }, { "Ser", "Ser", "Arg", "Arg" },
	},
	{
		{ "Val", "Val", "Val", "Val" }, { "Ala", "Ala", "Ala", "Ala" },
		{ "Asp", "Asp", "Asp", "Asp" }, { "Gly", "Gly", "Gly", "Gly" },
	},
};

int GetBaseIndex(char c)
{
	if (c == 'T' || c == 't') return 0;
	if (c == 'C' || c == 'c') return 1;
	if (c == 'A' || c == 'a') return 2;
	if (c == 'G' || c == 'g') return 3;
	throw std::runtime_error("Invalid base '" + std::string(1, c) + "'");
}

std::string CompBase(const std::string& base)
{
	if (base == "A" || base == "a") {
		return "T";
	} else if (base == "C" || base == "c") {
		return "G";
	} else if (base == "G" || base == "g") {
		return "C";
	} else if (base == "T" || base == "t") {
		return "t";
	} else {
		return base;
	}
}

std::string BaseToAA(const std::string& codon)
{
	return CODON_TABLE[GetBaseIndex(codon[0])][GetBaseIndex(codon[1])][GetBaseIndex(codon[2])];
}

std::string BaseToAA3(const std::string& codon)
{
	return CODON_TABLE_3[GetBaseIndex(codon[0])][GetBaseIndex(codon[1])][GetBaseIndex(codon[2])];
}

std::string GetMutType(const std::string& aa1, const std::string& aa2)
{
	if (aa1 == aa2) {
		return "Synonymous";
	} else if (aa1 == "*") {
		return "Stop-codon-loss";
	} else if (aa2 == "*") {
		return "Stop-codon-gain";
	} else {
		return "Non-synonymous";
	}
}

bool Convert(const std::string& chrom, const Transcript& trans, int pos, const std::string& ref, const std::string alt, const Fasta& fa)
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
	std::string codon1 = "";
	std::string codon2 = "";
	std::string mutAA = ".";
	std::string mutAA3 = ".";
	std::string mutType = ".";

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
					int mutPos = GetTxPos(exons, pos) - cdsStartTxPos;
					res = "c." + std::to_string(mutPos + 1);
					type = "Exon(" + std::to_string(i + 1) + "/" + std::to_string(exons.size()) + ")";
					type2 = "CDS";

					if (mutPos % 3 == 0) {
						std::string base1 = fa.GetSeq(chrom, trans.txStart_ + pos);

						int pos2 = pos + 1;
						int index = i;
						if (pos2 >= exons[i].second) {
							++index;
							assert(static_cast<size_t>(index) < exons.size());
							pos2 = exons[index].first;
						}
						std::string base2 = fa.GetSeq(chrom, trans.txStart_ + pos2);

						int pos3 = pos2 + 1;
						if (pos3 >= exons[index].second) {
							++index;
							assert(static_cast<size_t>(index) < exons.size());
							pos3 = exons[index].first;
						}
						std::string base3 = fa.GetSeq(chrom, trans.txStart_ + pos3);

						codon1 = base1 + base2 + base3;
						codon2 = alt + base2 + base3;
						std::string aa1 = BaseToAA(codon1);
						std::string aa2 = BaseToAA(codon2);
						mutAA = "p." + aa1 + std::to_string(mutPos / 3) + aa2;
						mutAA3 = "p." + BaseToAA3(codon1) + std::to_string(mutPos / 3) + BaseToAA3(codon2);
						mutType = GetMutType(aa1, aa2);
					} else if (mutPos % 3 == 1) {
						std::string base2 = fa.GetSeq(chrom, trans.txStart_ + pos);

						int pos2 = pos - 1;
						if (pos2 < exons[i].first) {
							assert(i > 0);
							pos2 = exons[i - 1].second - 1;
						}
						std::string base1 = fa.GetSeq(chrom, trans.txStart_ + pos2);

						int pos3 = pos + 1;
						if (pos3 >= exons[i].second) {
							assert(i + 1 < exons.size());
							pos3 = exons[i + 1].first;
						}
						std::string base3 = fa.GetSeq(chrom, trans.txStart_ + pos3);

						codon1 = base1 + base2 + base3;
						codon2 = base1 + alt + base3;
						std::string aa1 = BaseToAA(codon1);
						std::string aa2 = BaseToAA(codon2);
						mutAA = "p." + aa1 + std::to_string(mutPos / 3) + aa2;
						mutAA3 = "p." + BaseToAA3(codon1) + std::to_string(mutPos / 3) + BaseToAA3(codon2);
						mutType = GetMutType(aa1, aa2);
					} else {
						std::string base3 = fa.GetSeq(chrom, trans.txStart_ + pos);

						int pos2 = pos - 1;
						int index = i;
						if (pos2 < exons[index].first) {
							--index;
							assert(index >= 0);
							pos2 = exons[index].second - 1;
						}
						std::string base2 = fa.GetSeq(chrom, trans.txStart_ + pos2);

						int pos3 = pos2 - 1;
						if (pos3 < exons[index].first) {
							--index;
							assert(index >= 0);
							pos3 = exons[index].second - 1;
						}
						std::string base1 = fa.GetSeq(chrom, trans.txStart_ + pos3);

						codon1 = base1 + base2 + base3;
						codon2 = base1 + base2 + alt;
						std::string aa1 = BaseToAA(codon1);
						std::string aa2 = BaseToAA(codon2);
						mutAA = "p." + aa1 + std::to_string(mutPos / 3) + aa2;
						mutAA3 = "p." + BaseToAA3(codon1) + std::to_string(mutPos / 3) + BaseToAA3(codon2);
						mutType = GetMutType(aa1, aa2);
					}
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
					int mutPos = GetTxPosRev(exons, pos - 1) - cdsEndTxPos;
					res = "c." + std::to_string(mutPos + 1);
					type = "Exon(" + std::to_string(exons.size() - i + 1) + "/" + std::to_string(exons.size()) + ")";
					type2 = "CDS";

					if (mutPos % 3 == 2) {
						std::string base3 = CompBase(fa.GetSeq(chrom, trans.txStart_ + pos));

						int pos2 = pos + 1;
						int index = i - 1;
						if (pos2 >= exons[index].second) {
							++index;
							assert(static_cast<size_t>(index) < exons.size());
							pos2 = exons[index].first;
						}
						std::string base2 = CompBase(fa.GetSeq(chrom, trans.txStart_ + pos2));

						int pos3 = pos2 + 1;
						if (pos3 >= exons[index].second) {
							++index;
							assert(static_cast<size_t>(index) < exons.size());
							pos3 = exons[index].first;
						}
						std::string base1 = CompBase(fa.GetSeq(chrom, trans.txStart_ + pos3));

						codon1 = base1 + base2 + base3;
						codon2 = base1 + base2 + alt;
						std::string aa1 = BaseToAA(codon1);
						std::string aa2 = BaseToAA(codon2);
						mutAA = "p." + aa1 + std::to_string(mutPos / 3) + aa2;
						mutAA3 = "p." + BaseToAA3(codon1) + std::to_string(mutPos / 3) + BaseToAA3(codon2);
						mutType = GetMutType(aa1, aa2);
					} else if (mutPos % 3 == 1) {
						std::string base2 = CompBase(fa.GetSeq(chrom, trans.txStart_ + pos));
						int pos2 = pos - 1;
						int index = i - 1;
						if (pos2 < exons[index].first) {
							assert(index > 0);
							pos2 = exons[index - 1].second - 1;
						}
						std::string base3 = CompBase(fa.GetSeq(chrom, trans.txStart_ + pos2));

						int pos3 = pos + 1;
						index = i - 1;
						if (pos3 >= exons[index].second) {
							assert(index + 1 < static_cast<int>(exons.size()));
							pos3 = exons[index + 1].first;
						}
						std::string base1 = CompBase(fa.GetSeq(chrom, trans.txStart_ + pos3));

						codon1 = base1 + base2 + base3;
						codon2 = base1 + alt + base3;
						std::string aa1 = BaseToAA(codon1);
						std::string aa2 = BaseToAA(codon2);
						mutAA = "p." + aa1 + std::to_string(mutPos / 3) + aa2;
						mutAA3 = "p." + BaseToAA3(codon1) + std::to_string(mutPos / 3) + BaseToAA3(codon2);
						mutType = GetMutType(aa1, aa2);
					} else {
						std::string base1 = CompBase(fa.GetSeq(chrom, trans.txStart_ + pos));

						int pos2 = pos - 1;
						int index = i - 1;
						if (pos2 < exons[index].first) {
							--index;
							assert(index >= 0);
							pos2 = exons[index].second - 1;
						}
						std::string base2 = CompBase(fa.GetSeq(chrom, trans.txStart_ + pos2));

						int pos3 = pos2 - 1;
						if (pos3 < exons[index].first) {
							--index;
							assert(index >= 0);
							pos3 = exons[index].second - 1;
						}
						std::string base3 = CompBase(fa.GetSeq(chrom, trans.txStart_ + pos3));

						codon1 = base1 + base2 + base3;
						codon2 = alt + base2 + base3;
						std::string aa1 = BaseToAA(codon1);
						std::string aa2 = BaseToAA(codon2);
						mutAA = "p." + aa1 + std::to_string(mutPos / 3) + aa2;
						mutAA3 = "p." + BaseToAA3(codon1) + std::to_string(mutPos / 3) + BaseToAA3(codon2);
						mutType = GetMutType(aa1, aa2);
					}
					break;
				}
			}
		}
	}

	//assert(ref == seq.substr(pos, 1));
	res += fa.GetSeq(chrom, trans.txStart_ + pos, 1) + ">" + alt;

	std::cout << chrom << "\t" << trans.txStart_ + pos + 1 << "\t" << trans.strand_ << "\t" << trans.name_ << "\t" << trans.name2_
		<< "\t" << res << "\t" << type << "\t" << type2 << "\t" << codon1 << "\t" << codon2 << "\t" << mutAA << "\t" << mutAA3 << "\t" << mutType << std::endl;
	return true;
}

bool ProcessItem(const std::string& chrom, int pos, const std::string& alleleRef, const std::string& alleleAlt,
		const std::map<std::string, std::vector<Transcript>>& data, const Fasta& fa)
{
	auto it = data.find(chrom);
	if (it != data.end()) {
		auto trans = it->second;
		bool found = false;
		for (size_t i = 0; i < trans.size(); ++i) {
			if (pos >= trans[i].txStart_ && pos < trans[i].txEnd_) {
				//std::cerr << "Found in " << trans[i].name_ << " (" << chrom << ':' << trans[i].txStart_ << "-" << trans[i].txEnd_ << ")" << std::endl;

				if (!fa.Has(chrom)) {
					std::cerr << "Can not found sequence '" << chrom << "' in ref fasta" << std::endl;
					return false;
				}

				Convert(chrom, trans[i], pos - trans[i].txStart_, alleleRef, alleleAlt, fa);
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
		const std::map<std::string, std::vector<Transcript>>& data, const Fasta& fa)
{
	std::ifstream file(filename, std::ios::in);
	if (!file.is_open()) {
		std::cerr << "Error: Can not open file '" << filename << "'!" << std::endl;
		return false;
	}

	std::cout << "chrom\tpos\tstrand\tname\tname2\tmutate\tsegment\ttype\tcodon1\tcodon2\tmutAA\tmutAA3\tmutType" << std::endl;

	size_t lineNo = 0;
	std::string line;
	while (std::getline(file, line)) {
		++lineNo;
		if (line.empty() || line[0] == '#') continue;

		LineSplit sp;
		sp.Split(line, '\t');

		try {
			std::string chrom = sp.GetField(0);
			int genomePos = stoi(sp.GetField(1));
			std::string alleleRef = sp.GetField(3);
			std::string alleleAlt = sp.GetField(4);

			ProcessItem(chrom, genomePos - 1, alleleRef, alleleAlt, data, fa);
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
		"Usage:  crabber annotate <x.vcf> <refGene.tsv> <ref.fa>\n"
		"\n"
		"Input:\n"
		"   <x.vcf>         input SNV list in VCF format\n"
		"   <refGene.tsv>   track data downloaded from UCSC table browser\n"
		"   <ref.fa>        reference genome in FASTA format\n"
		<< std::endl;
}

int main_annotate(int argc, char* const argv[])
{
	std::vector<std::string> args(argv, argv + argc);
	if (args.size() < 4) {
		PrintUsage();
		return 1;
	}

	Fasta fa;
	if (!fa.Load(args[3])) {
		return 1;
	}

	std::map<std::string, std::vector<Transcript>> data;
	if (!LoadRefGene(args[2], data)) {
		return 1;
	}

	if (!Process(args[1], data, fa)) {
		return 1;
	}
	return 0;
}
