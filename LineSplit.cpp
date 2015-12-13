#include "LineSplit.h"

std::string LineSplit::GetField(size_t index) const
{
	if (index >= 0 && index + 1 < pos_.size() && (maxCount_ == 0 || index + 1 < maxCount_)) {
		return std::string(pos_[index], pos_[index + 1] - 1);
	} else if (index < pos_.size()) {
		return std::string(pos_[index]);
	} else {
		return "";
	}
}

size_t LineSplit::Split(const std::string& s, char sep, size_t maxCount)
{
	pos_.clear();
	maxCount_ = maxCount;

	pos_.push_back(s.c_str());
	for (const char *p = s.c_str(); maxCount == 0 || pos_.size() <= maxCount; ++p) {
		if (*p == sep) {
			pos_.push_back(p + 1);
		} else if (*p == '\r' || *p == '\n' || *p == '\0') {
			pos_.push_back(p + 1);
			break;
		}
	}
	return pos_.size() - 1;
}
