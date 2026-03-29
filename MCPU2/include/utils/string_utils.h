#ifndef UTILS_STRING_H
#define UTILS_STRING_H

#include <algorithm>
#include <cctype>
#include <string>

// Trims whitespace from both ends of the string
inline std::string trim(const std::string &s) {
    auto start = s.find_first_not_of(" \t\n\r\f\v");
    auto end   = s.find_last_not_of(" \t\n\r\f\v");
    return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
}

#endif  // UTILS_STRING_H