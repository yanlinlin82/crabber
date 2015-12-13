#ifndef __TRANSCRIPT_H__
#define __TRANSCRIPT_H__

#include <string>
#include <vector>
#include <utility>

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

#endif
