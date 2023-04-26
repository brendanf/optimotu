#ifndef OPTIMOTU_PADSTRINGS_H_INCLUDED
#define OPTIMOTU_PADSTRINGS_H_INCLUDED

#include <vector>
#include <string>
#include <memory>

std::shared_ptr<char[]> pad_strings(
    const std::vector<std::string> &seq,
    std::size_t & seq_width
);

#endif //OPTIMOTU_PADSTRINGS_H_INCLUDED
