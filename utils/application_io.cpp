#include "application_io.hpp"

std::vector<std::vector<uint32_t>> parse_input_file(const std::string& file_path, pcfg* grammar){
    std::vector<std::vector<uint32_t>> sentences;
    std::string line;
	std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the input file at path: " << file_path << std::endl;
        return sentences;
    }
    int N = grammar->N();

    while (std::getline(file, line)) {
        if(line == "")
            continue;
        std::vector<uint32_t> input_words;

        std::string word;
        std::stringstream line_string_stream(line);
        while (getline(line_string_stream, word, ' ')) {
            input_words.push_back(grammar->terminate_map.find(std::string("\'") + word + std::string("\'"))->second + N);
        }
        sentences.push_back(input_words);
    }
    return sentences;
}