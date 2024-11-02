#include<iostream>
#include <vector>
#include <fstream>
#include <string>
#include <bits/stdc++.h>
#include <map>
#include "utils/application_io.hpp"
#include "grammar/grammar_parser.hpp"


std::vector<std::vector<uint32_t>> read_input_file(const std::string& file_path){
    std::vector<std::vector<uint32_t>> sentences;
    std::string line;
	std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the input file at path: " << file_path << std::endl;
        return sentences;
    }

    while (std::getline(file, line)) {
        if(line == "")
            continue;
        std::vector<uint32_t> input_words;

        std::string word;
        std::stringstream line_string_stream(line);
        while (getline(line_string_stream, word, ' ')) {
            uint32_t word_id = std::stoi(word);            
            input_words.push_back(word_id);
        }
        sentences.push_back(input_words);
    }
    return sentences;
}

int main(int argc, char* argv[]){
    std::string grammar_filename = argc > 1 ? std::string(argv[1]) : "grammar.pcfg";
    std::string input_filename = argc > 2 ? std::string(argv[2]) : "eeg_sentences.txt";
    std::string output_filename = argc > 3 ? std::string(argv[3]) : "sentences_converted.txt";
    int delay = 2;    
    pcfg* grammar = prepare_grammar(grammar_filename);
    std::vector<std::vector<uint32_t>> sentences = read_input_file(input_filename);
    int N = grammar->N();
    for(auto& sentence : sentences){
        size_t sentence_length = sentence.size();
        std::vector<uint32_t> phase_array;
        for(int begin = delay - 1; begin < sentence_length; begin ++){
            uint32_t s = 0;
            for(int j = 0; j < delay; j++){
                s *= 4;
                s += (sentence[begin - j] - 1);
            }
            phase_array.emplace_back(s + 1);
        }
        assert(phase_array.size() > 1 && (phase_array[0] == phase_array[phase_array.size() - 1]));
        sentence.clear();
        sentence.assign(phase_array.begin(), phase_array.end());
    }
    std::ofstream output_file(output_filename);
    for(int sentence_id = 0; sentence_id < sentences.size(); sentence_id++){
        std::vector<uint32_t>& sentence = sentences[sentence_id];
        size_t sentence_length = sentence.size();
        std::vector<uint32_t> phase_array;
        for(int begin = 0; begin < sentence_length; begin ++){
            uint32_t s = sentence[begin];
            output_file << s;
            if(begin != sentence_length - 1){
                output_file << " ";
            }
        }
        if(sentence_id != sentences.size() - 1)
            output_file << std::endl;
        assert(sentence.size() > 1 && (sentence[0] == sentence[sentence.size() - 1]));
        
    }
    return 0;
}
