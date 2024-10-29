#include "application_io.hpp"
#include "utils/string_helper.hpp"

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


void save_data_set_to_file(std::string file_path, std::vector<std::vector<uint32_t>> sentences, pcfg* grammar){
    std::ofstream output_file(file_path);
    if (!output_file.is_open()) {
        std::cerr << "Error opening file: " << file_path << std::endl;
        return;
    }
    int N = grammar->N();
    for(size_t i = 0; i < sentences.size(); i++){
        const auto& sentence = sentences[i];
        size_t sequence_length = sentence.size();
        for(size_t j = 0; j < sequence_length; j++){
            output_file << remove_quotes(grammar->reversed_terminate_map[sentence[j] - N]);
            if(j != sequence_length - 1){
                output_file << " ";
            }
        }
        output_file << std::endl;
    }
}