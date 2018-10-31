#include <iostream>
 
#include "antlr4-runtime.h"
#include "circuit_specificationLexer.h"
#include "circuit_specificationParser.h"
 
using namespace std;
using namespace antlr4;
 
int main(int argc, const char* argv[]) {
    std::ifstream stream;
    stream.open("circuit_specification.txt");
    
    ANTLRInputStream input(stream);
    circuit_specificationLexer lexer(&input);
    CommonTokenStream tokens(&lexer);
    circuit_specificationParser parser(&tokens);    
 	
 	cout << "Hello world" << endl;
    return 0;
}