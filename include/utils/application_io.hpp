#ifndef HPP_APPLICATION_IO
#define HPP_APPLICATION_IO
#include <vector>
#include <fstream>
#include <string>
#include <bits/stdc++.h>
#include <map>
#include "grammar/grammar.hpp"
std::vector<std::vector<uint32_t>> parse_input_file(const std::string& file_path, pcfg* grammar, int limit = -1, int max_sentence_length = 500);
void save_data_set_to_file(std::string file_path, std::vector<std::vector<uint32_t>> sentences, pcfg* grammar);
#endif
