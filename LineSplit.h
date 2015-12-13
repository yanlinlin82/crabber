#ifndef __LINE_SPLIT_H__
#define __LINE_SPLIT_H__

#include <vector>
#include <string>

class LineSplit
{
public:
	LineSplit(): maxCount_(0) { }
public:
	size_t Split(const std::string& s, char sep = '\t', size_t maxCount = 0);
	std::string GetField(size_t index) const;
private:
	std::vector<const char*> pos_;
	size_t maxCount_;
};

#endif
