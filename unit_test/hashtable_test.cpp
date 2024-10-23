#define BOOST_TEST_MODULE GRAMMAR_TABLE_UNITEST_Module
#include <boost/test/included/unit_test.hpp>
#include "grammar.hpp"
float test_hashtable(pcfg* grammar) {
    return 0.0f;
}

BOOST_AUTO_TEST_CASE(test_hashtable) {
    pcfg* grammar = prepare_grammar("grammar.pcfg");
    test_hashtable(grammar);
    BOOST_CHECK_EQUAL(result, 0.0f);
}