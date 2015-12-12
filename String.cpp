#include "String.h"

std::string TrimLeft(const std::string& str, const std::string& drop)
{
	std::string::size_type pos = str.find_first_not_of(drop);
	if (pos == std::string::npos) {
		return "";
	}
	return str.substr(pos);
}

std::string TrimRight(const std::string& str, const std::string& drop)
{
	std::string::size_type pos = str.find_last_not_of(drop);
	if (pos == std::string::npos) {
		return "";
	}
	return str.substr(0, pos + 1);
}

std::string Trim(const std::string& str, const std::string& drop)
{
	return TrimLeft(TrimRight(str, drop), drop);
}
