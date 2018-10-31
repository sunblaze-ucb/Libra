JAVA=java

GRAMMAR=circuit_specification.g4
GENERATED=generated

CC=g++
LIBS=../../lib/libantlr4-runtime.a
JAR=../../lib/antlr-4.7.1-complete.jar
CCARGS=-c -I ../../include -I $(GENERATED) -std=c++11
LDARGS=-g

OUTPUT=output

all: front_end
front_end: dirs antlr4 front_end.cpp
	$(CC) $(CCARGS) $(GENERATED)/circuit_specificationLexer.cpp -o $(OUTPUT)/circuit_specificationLexer.o
	$(CC) $(CCARGS) $(GENERATED)/circuit_specificationParser.cpp -o $(OUTPUT)/circuit_specificationParser.o
	$(CC) $(CCARGS) front_end.cpp -o $(OUTPUT)/front_end.o
	$(CC) $(LDARGS) $(OUTPUT)/front_end.o $(OUTPUT)/circuit_specificationLexer.o $(OUTPUT)/circuit_specificationParser.o $(LIBS) -o front_end

antlr4:$(GRAMMAR)
	$(JAVA) -jar $(JAR) -Dlanguage=Cpp -o $(GENERATED) $(GRAMMAR)

dirs:
	mkdir -p $(OUTPUT)
	mkdir -p $(GENERATED)

clean:
	rm -rf $(OUTPUT)
	rm -rf $(GENERATED)