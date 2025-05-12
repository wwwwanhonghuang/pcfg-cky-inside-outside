#ifndef H_GRAMMAR_PARSER
#define H_GRAMMAR_PARSER
#include <string>
#include "grammar.hpp"
pcfg* prepare_grammar(const std::string& path, bool original_in_log_form);
#endif