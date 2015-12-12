#ifndef __FASTA_H__
#define __FASTA_H__

#include <map>
#include <string>

class Fasta
{
public:
	bool Load(const std::string& filename, bool verbose = false);

	bool Has(const std::string& chrom) const;

	std::string GetSeq(const std::string& chrom, size_t pos, size_t size = 1) const;
private:
	std::map<std::string, std::string> seq_;
};

#endif
