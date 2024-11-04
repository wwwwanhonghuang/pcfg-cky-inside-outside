#include "utils/string_helper.hpp"
std::string remove_quotes(const std::string& str) {
    if (!str.empty() && str.front() == '\'' && str.back() == '\'') {
        return str.substr(1, str.size() - 2);
    }
    return str;
}
