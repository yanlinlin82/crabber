#ifndef __STR_FUNC_H__
#define __STR_FUNC_H__

#include <string>

std::string TrimLeft(const std::string& str, const std::string& drop = " \t");
std::string Trim_right(const std::string& str, const std::string& drop = " \t");
std::string Trim(const std::string& str, const std::string& drop = " \t");

#endif
