//
// Created by 徐溶延 on 2020/11/18.
//

#ifndef MESH_CUTTING_STRING_UTILS_H
#define MESH_CUTTING_STRING_UTILS_H
#include <string>
#include <vector>

namespace common {
std::string findSuffix(const std::string &s);

std::string findPrefix(const std::string &s);

void trimLeftTrailingSpaces(std::string &input);

void trimRightTrailingSpaces(std::string &input);

std::vector<int> string2IntegerVector(std::string &str);
}


#endif //MESH_CUTTING_STRING_UTILS_H
