#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cassert>
#include "Annotate.h"
#include "String.h"
#include "Fasta.h"
#include "LineSplit.h"
#include "Transcript.h"

static std::vector<std::string> Split(const std::string& s, const std::string& sep = "\t ", size_t count = 0)
{
	std::vector<std::string> a;
	size_t lastPos = s.find_first_not_of(sep);
	while (lastPos != std::string::npos && (count == 0 || a.size() < count)) {
		if (a.size() + 1 == count) {
			a.push_back(s.substr(lastPos));
		} else {
			size_t pos = s.find_first_of(sep, lastPos);
			if (pos == std::string::npos) {
				a.push_back(s.substr(lastPos));
			} else {
				a.push_back(s.substr(lastPos, pos - lastPos));
			}
			lastPos = s.find_first_not_of(sep, pos);
		}
	}
	return a;
}

bool LoadRefGene(const std::string& filename, std::map<std::string, std::vector<Transcript>>& data)
{
	std::ifstream file(filename, std::ios::in);
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
		{ "Asp", "Asp", "Glu", "Glu" }, { "Gly", "Gly", "Gly", "Gly" },
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
		return "A";
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

bool Convert(const std::string& chrom, const Transcript& trans, int pos, const std::string& ref, const std::string alt,
		const Fasta& fa, const std::vector<std::string>& fields)
{
	std::string res = ".";
	std::string type = ".";
	std::string type2 = ".";
	std::string codon1 = ".";
	std::string codon2 = ".";
	std::string mutAA = ".";
	std::string mutAA3 = ".";
	std::string mutType = ".";

	const auto& exons = trans.exons_;
	int cdsStart = trans.cdsStart_ - trans.txStart_;
	int cdsEnd = trans.cdsEnd_ - trans.txStart_;

	if (cdsStart == cdsEnd) {
		mutType = "Unknown";
	} else if (trans.strand_ == "+") {
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

					if (alt == "*") {
						mutType = "Frameshift";
					} else if (alt == ".") {
						mutType = "Unknown";
					} else if (alt.size() > 1) {
						mutType = ((alt.size() - 2) % 3 == 0) ? "Inframe" : "Frameshift";
					} else if (mutPos % 3 == 0) {
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
						mutAA = "p." + aa1 + std::to_string(mutPos / 3 + 1) + aa2;
						mutAA3 = "p." + BaseToAA3(codon1) + std::to_string(mutPos / 3 + 1) + BaseToAA3(codon2);
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
						mutAA = "p." + aa1 + std::to_string(mutPos / 3 + 1) + aa2;
						mutAA3 = "p." + BaseToAA3(codon1) + std::to_string(mutPos / 3 + 1) + BaseToAA3(codon2);
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
						mutAA = "p." + aa1 + std::to_string(mutPos / 3 + 1) + aa2;
						mutAA3 = "p." + BaseToAA3(codon1) + std::to_string(mutPos / 3 + 1) + BaseToAA3(codon2);
						mutType = GetMutType(aa1, aa2);
					}
					break;
				}
			}
		}
	} else {
		assert(trans.strand_ == "-");
		int cdsStartTxPos = GetTxPos(exons, cdsStart);
		int cdsEndTxPos = GetTxPos(exons, cdsEnd - 1) + 1;

		for (size_t i = exons.size(); i > 0; --i) {
			int start = exons[i - 1].first;
			int end = exons[i - 1].second;
			if (pos >= cdsEnd) { // 5'-UTR
				if (pos >= end) { // intron
					if (i == exons.size() || (pos - end) < (exons[i].first - pos)) {
						res = "c.-" + std::to_string(GetTxPos(exons, end - 1) - cdsEndTxPos + 1) + "-" + std::to_string(pos - end + 1);
					} else {
						res = "c.-" + std::to_string(GetTxPos(exons, exons[i].first) - cdsEndTxPos + 1) + "+" + std::to_string(exons[i].first - pos);
					}
					type = "Intron(" + std::to_string(exons.size() - i) + "/" + std::to_string(exons.size() - 1) + ")";
					break;
				} else if (pos >= start) { // exon
					res = "c.-" + std::to_string(GetTxPos(exons, pos) - cdsEndTxPos + 1);
					type = "Exon(" + std::to_string(exons.size() - i + 1) + "/" + std::to_string(exons.size()) + ")";
					type2 = "5'-UTR";
					break;
				}
			} else if (pos < cdsStart) { // 3'-UTR
				if (pos >= end) { // intron
					if (i == exons.size() || (pos - end) < (exons[i].first - pos)) {
						res = "c.*" + std::to_string(cdsStartTxPos - GetTxPos(exons, end - 1)) + "-" + std::to_string(pos - end + 1);
					} else {
						res = "c.*" + std::to_string(cdsStartTxPos - GetTxPos(exons, exons[i].first)) + "+" + std::to_string(exons[i].first - pos);
					}
					type = "Intron(" + std::to_string(exons.size() - i) + "/" + std::to_string(exons.size() - 1) + ")";
					break;
				} else if (pos >= start) { // exon
					res = "c.*" + std::to_string(cdsStartTxPos - GetTxPos(exons, pos));
					type = "Exon(" + std::to_string(exons.size() - i + 1) + "/" + std::to_string(exons.size()) + ")";
					type2 = "3'-UTR";
					break;
				}
			} else { // CDS
				if (pos >= end) { // intron
					if (i == exons.size() || (pos - end) < (exons[i].first - pos)) {
						res = "c." + std::to_string(cdsEndTxPos - GetTxPos(exons, end - 1)) + "-" + std::to_string(pos - end + 1);
					} else {
						res = "c." + std::to_string(cdsEndTxPos - GetTxPos(exons, exons[i].first)) + "+" + std::to_string(exons[i].first - pos);
					}
					type = "Intron(" + std::to_string(exons.size() - i) + "/" + std::to_string(exons.size() - 1) + ")";
					break;
				} else if (pos >= start) { // exon
					int mutPos = cdsEndTxPos - GetTxPos(exons, pos) - 1;
					res = "c." + std::to_string(mutPos + 1);
					type = "Exon(" + std::to_string(exons.size() - i + 1) + "/" + std::to_string(exons.size()) + ")";
					type2 = "CDS";

					if (alt == "*") {
						mutType = "Frameshift";
					} else if (alt == ".") {
						mutType = "Unknown";
					} else if (alt.size() > 1) {
						mutType = ((alt.size() - 2) % 3 == 0) ? "Inframe" : "Frameshift";
					} else if (mutPos % 3 == 2) {
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
						codon2 = base1 + base2 + CompBase(alt);
						std::string aa1 = BaseToAA(codon1);
						std::string aa2 = BaseToAA(codon2);
						mutAA = "p." + aa1 + std::to_string(mutPos / 3 + 1) + aa2;
						mutAA3 = "p." + BaseToAA3(codon1) + std::to_string(mutPos / 3 + 1) + BaseToAA3(codon2);
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
						codon2 = base1 + CompBase(alt) + base3;
						std::string aa1 = BaseToAA(codon1);
						std::string aa2 = BaseToAA(codon2);
						mutAA = "p." + aa1 + std::to_string(mutPos / 3 + 1) + aa2;
						mutAA3 = "p." + BaseToAA3(codon1) + std::to_string(mutPos / 3 + 1) + BaseToAA3(codon2);
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
						codon2 = CompBase(alt) + base2 + base3;
						std::string aa1 = BaseToAA(codon1);
						std::string aa2 = BaseToAA(codon2);
						mutAA = "p." + aa1 + std::to_string(mutPos / 3 + 1) + aa2;
						mutAA3 = "p." + BaseToAA3(codon1) + std::to_string(mutPos / 3 + 1) + BaseToAA3(codon2);
						mutType = GetMutType(aa1, aa2);
					}
					break;
				}
			}
		}
	}

	//assert(ref == seq.substr(pos, 1));
	res += fa.GetSeq(chrom, trans.txStart_ + pos, ref.size()) + ">" + alt;

	for (size_t i = 0; i + 1 < fields.size(); ++i) {
		std::cout << fields[i] << "\t";
	}
	std::cout << trans.strand_ << "\t" << trans.name_ << "\t" << trans.name2_
		<< "\t" << res << "\t" << type << "\t" << type2 << "\t" << codon1 << "\t" << codon2 << "\t" << mutAA << "\t" << mutAA3 << "\t" << mutType
		<< "\t" << fields.back() << std::endl;
	return true;
}

bool ProcessItem(const std::string& chrom, int pos, const std::string& alleleRef, const std::string& alleleAlt,
		const std::map<std::string, std::vector<Transcript>>& data, const Fasta& fa,
		const std::vector<std::string>& fields, bool outputFirstOnly)
{
	auto it = data.find(chrom);
	if (it != data.end()) {
		auto trans = it->second;
		bool found = false;
		for (size_t i = 0; i < trans.size(); ++i) {
			if (pos >= trans[i].txStart_ && pos < trans[i].txEnd_) {
				if (!fa.Has(chrom)) {
					std::cerr << "Can not found sequence '" << chrom << "' in ref fasta" << std::endl;
					return false;
				}

				Convert(chrom, trans[i], pos - trans[i].txStart_, alleleRef, alleleAlt, fa, fields);
				found = true;
				if (outputFirstOnly) {
					break;
				}
			}
		}
		if (!found) {
			for (size_t i = 0; i + 1 < fields.size(); ++i) {
				std::cout << fields[i] << "\t";
			}
			std::cout << ".\t.\t.\t.\tIntergenic\t.\t.\t.\t.\t.\t.\t" << fields.back() << std::endl;
		}
	}
	return true;
}

static void OutputHeader(const std::vector<std::string>& fields, size_t insertPos)
{
	for (size_t i = 0; i < insertPos && i < fields.size(); ++i) {
		std::cout << fields[i] << '\t';
	}
	std::cout << "strand\tname\tname2\tmutate\tsegment\ttype\tcodon1\tcodon2\tmutAA\tmutAA3\tmutType";
	for (size_t i = insertPos; i < fields.size(); ++i) {
		std::cout << '\t' << fields[i];
	}
	std::cout << std::endl;
}

static bool Process(const std::string& filename, bool tsvFile, bool hasHeader,
		const std::map<std::string, std::vector<Transcript>>& data, const Fasta& fa, bool outputFirstOnly)
{
	std::ifstream file(filename, std::ios::in);
	if (!file.is_open()) {
		std::cerr << "Error: Can not open file '" << filename << "'!" << std::endl;
		return false;
	}

	std::vector<std::string> fields;
	size_t lineNo = 0;
	std::string line;
	while (std::getline(file, line)) {
		++lineNo;
		if (tsvFile) {
			if (line.empty() || line[0] == '#') continue;
			if (lineNo == 1 && hasHeader) {
				fields = Split(line, "\t");
				OutputHeader(fields, 5);
				continue;
			}
		} else {
			if (line.empty()) continue;
			if (line[0] == '#') {
				if (line[1] == '#') continue;
				if (hasHeader) {
					fields = Split(line.substr(1), "\t");
					OutputHeader(fields, 5);
				}
				continue;
			}
		}

		try {
			fields = Split(line, "\t", 6);

			std::string chrom = fields[0];
			int genomePos = std::stoi(fields[1]);
			std::string alleleRef = fields[3];
			std::string alleleAlt = fields[4];

			ProcessItem(chrom, genomePos - 1, alleleRef, alleleAlt, data, fa, fields, outputFirstOnly);
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
		"Usage:  crabber annotate [options] <x.vcf> <refGene.tsv> <ref.fa>\n"
		"\n"
		"Input:\n"
		"   <x.vcf>         input SNV list in VCF format\n"
		"   <refGene.tsv>   track data downloaded from UCSC table browser\n"
		"   <ref.fa>        reference genome in FASTA format\n"
		"\n"
		"Options:\n"
		"   -1              output only first matched script, default to output all\n"
		"   -T              TSV input, with columns: chrom, start, end, ref, alt...\n"
		"   -H              input file has header, output with header\n"
		<< std::endl;
}

int Annotate_main(int argc, char* const argv[])
{
	bool outputFirstOnly = false;
	std::string inputFile;
	std::string refGeneFile;
	std::string refFastaFile;
	bool tsvInput = false;
	bool hasHeader = false;

	std::vector<std::string> args(argv, argv + argc);
	std::vector<std::string> restArgs;
	for (size_t i = 1; i < args.size(); ++i) {
		if (args[i] == "-1") {
			outputFirstOnly = true;
		} else if (args[i] == "-T") {
			tsvInput = true;
		} else if (args[i] == "-H") {
			hasHeader = true;
		} else {
			restArgs.push_back(args[i]);
		}
	}
	if (restArgs.size() < 3) {
		PrintUsage();
		return 1;
	}
	inputFile = restArgs[0];
	refGeneFile = restArgs[1];
	refFastaFile = restArgs[2];

	Fasta fa;
	if (!fa.Load(refFastaFile)) {
		return 1;
	}

	std::map<std::string, std::vector<Transcript>> data;
	if (!LoadRefGene(refGeneFile, data)) {
		return 1;
	}

	if (!Process(inputFile, tsvInput, hasHeader, data, fa, outputFirstOnly)) {
		return 1;
	}
	return 0;
}
