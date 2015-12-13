#include "Transcript.h"
#include "LineSplit.h"

bool Transcript::SetExons(int exonCount, const std::string& exonStarts, const std::string& exonEnds)
{
	std::vector<size_t> pos, pos2;
	pos.resize(exonCount + 1);
	pos2.resize(exonCount + 1);

	LineSplit sp, sp2;
	sp.Split(exonStarts, ',');
	sp2.Split(exonEnds, ',');
	for (int i = 0; i < exonCount; ++i) {
		exons_.push_back(std::make_pair(
					stoi(sp.GetField(i)) - txStart_,
					stoi(sp2.GetField(i)) - txStart_));
	}
	return true;
}
