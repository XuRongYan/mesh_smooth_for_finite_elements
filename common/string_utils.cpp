//
// Created by 徐溶延 on 2020/11/18.
//

#include "string_utils.h"
#include <sstream>
#include <algorithm>

std::string common::findSuffix(const std::string &s) {
	size_t suffix_id = s.rfind('.');
	return s.substr(suffix_id);
}

std::string common::findPrefix(const std::string &s) {
	size_t suffix_id = s.rfind('.');
    int divide = s.rfind('/');
    if (divide == -1) divide = 0;
	return s.substr(divide + 1, suffix_id);
}

void common::trimLeftTrailingSpaces(std::string &input) {
    input.erase(input.begin(), std::find_if(input.begin(), input.end(), [](int ch) {
		return !isspace(ch);
	}));
}

void common::trimRightTrailingSpaces(std::string &input) {
    input.erase(std::find_if(input.rbegin(), input.rend(), [](int ch) {
		return !isspace(ch);
	}).base(), input.end());
}


std::vector<int> common::string2IntegerVector(std::string &str) {
	std::vector<int> output;
	trimLeftTrailingSpaces(str);
	trimRightTrailingSpaces(str);
	std::stringstream ss;
	ss.str(str);
	std::string item;
	char delim = ',';
	while (getline(ss, item, delim)) {
		output.emplace_back(std::stoi(item));
	}
	return output;
}




