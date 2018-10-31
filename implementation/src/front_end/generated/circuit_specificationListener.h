
// Generated from circuit_specification.g4 by ANTLR 4.7.1

#pragma once


#include "antlr4-runtime.h"
#include "circuit_specificationParser.h"


/**
 * This interface defines an abstract listener for a parse tree produced by circuit_specificationParser.
 */
class  circuit_specificationListener : public antlr4::tree::ParseTreeListener {
public:

  virtual void enterProg(circuit_specificationParser::ProgContext *ctx) = 0;
  virtual void exitProg(circuit_specificationParser::ProgContext *ctx) = 0;

  virtual void enterExpr(circuit_specificationParser::ExprContext *ctx) = 0;
  virtual void exitExpr(circuit_specificationParser::ExprContext *ctx) = 0;


};

