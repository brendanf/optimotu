#ifndef OPTIMOTU_PADSTRINGS_H_INCLUDED
#define OPTIMOTU_PADSTRINGS_H_INCLUDED

#include <utility>
#include <vector>
#include <string>
#include <cstring>

std::pair<char*, size_t> pad_strings(const std::vector<std::string> &seq);

#endif //OPTIMOTU_PADSTRINGS_H_INCLUDED
