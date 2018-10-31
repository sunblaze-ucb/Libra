
// Generated from circuit_specification.g4 by ANTLR 4.7.1

#pragma once


#include "antlr4-runtime.h"
#include "circuit_specificationListener.h"


/**
 * This class provides an empty implementation of circuit_specificationListener,
 * which can be extended to create a listener which only needs to handle a subset
 * of the available methods.
 */
class  circuit_specificationBaseListener : public circuit_specificationListener {
public:

  virtual void enterProg(circuit_specificationParser::ProgContext * /*ctx*/) override { }
  virtual void exitProg(circuit_specificationParser::ProgContext * /*ctx*/) override { }

  virtual void enterExpr(circuit_specificationParser::ExprContext * /*ctx*/) override { }
  virtual void exitExpr(circuit_specificationParser::ExprContext * /*ctx*/) override { }


  virtual void enterEveryRule(antlr4::ParserRuleContext * /*ctx*/) override { }
  virtual void exitEveryRule(antlr4::ParserRuleContext * /*ctx*/) override { }
  virtual void visitTerminal(antlr4::tree::TerminalNode * /*node*/) override { }
  virtual void visitErrorNode(antlr4::tree::ErrorNode * /*node*/) override { }

};

