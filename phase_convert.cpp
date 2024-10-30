#include<iostream>
#include <vector>
#include <fstream>
#include <string>
#include <bits/stdc++.h>
#include <map>
#include "utils/application_io.hpp"
#include "grammar/grammar_parser.hpp"


int main(int argc, char* argv[]){
    std::string grammar_filename = argc > 1 ? std::string(argv[1]) : "grammar.pcfg";
    std::string input_filename = argc > 2 ? std::string(argv[2]) : "eeg_sentences_test.txt";
    std::string output_filename = argc > 3 ? std::string(argv[3]) : "sentences_converted.txt";
    int delay = 2;    
    pcfg* grammar = prepare_grammar(grammar_filename);
    std::vector<std::vector<uint32_t>> sentences = parse_input_file(input_filename, grammar);
    int N = grammar->N();
    for(auto& sentence : sentences){
        size_t sentence_length = sentence.size();
        std::vector<uint32_t> phase_array;
        for(int begin = delay - 1; begin < sentence_length; begin ++){
            uint32_t s = 0;
            for(int j = 0; j < delay; j++){
                s *= 4;
                s += sentence[begin - j] - N;
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