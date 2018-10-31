
// Generated from circuit_specification.g4 by ANTLR 4.7.1


#include "circuit_specificationListener.h"

#include "circuit_specificationParser.h"


using namespace antlrcpp;
using namespace antlr4;

circuit_specificationParser::circuit_specificationParser(TokenStream *input) : Parser(input) {
  _interpreter = new atn::ParserATNSimulator(this, _atn, _decisionToDFA, _sharedContextCache);
}

circuit_specificationParser::~circuit_specificationParser() {
  delete _interpreter;
}

std::string circuit_specificationParser::getGrammarFileName() const {
  return "circuit_specification.g4";
}

const std::vector<std::string>& circuit_specificationParser::getRuleNames() const {
  return _ruleNames;
}

dfa::Vocabulary& circuit_specificationParser::getVocabulary() const {
  return _vocabulary;
}


//----------------- ProgContext ------------------------------------------------------------------

circuit_specificationParser::ProgContext::ProgContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

std::vector<circuit_specificationParser::ExprContext *> circuit_specificationParser::ProgContext::expr() {
  return getRuleContexts<circuit_specificationParser::ExprContext>();
}

circuit_specificationParser::ExprContext* circuit_specificationParser::ProgContext::expr(size_t i) {
  return getRuleContext<circuit_specificationParser::ExprContext>(i);
}

std::vector<tree::TerminalNode *> circuit_specificationParser::ProgContext::NEWLINE() {
  return getTokens(circuit_specificationParser::NEWLINE);
}

tree::TerminalNode* circuit_specificationParser::ProgContext::NEWLINE(size_t i) {
  return getToken(circuit_specificationParser::NEWLINE, i);
}


size_t circuit_specificationParser::ProgContext::getRuleIndex() const {
  return circuit_specificationParser::RuleProg;
}

void circuit_specificationParser::ProgContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<circuit_specificationListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterProg(this);
}

void circuit_specificationParser::ProgContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<circuit_specificationListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitProg(this);
}

circuit_specificationParser::ProgContext* circuit_specificationParser::prog() {
  ProgContext *_localctx = _tracker.createInstance<ProgContext>(_ctx, getState());
  enterRule(_localctx, 0, circuit_specificationParser::RuleProg);
  size_t _la = 0;

  auto onExit = finally([=] {
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(9);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == circuit_specificationParser::T__4

    || _la == circuit_specificationParser::INT) {
      setState(4);
      expr(0);
      setState(5);
      match(circuit_specificationParser::NEWLINE);
      setState(11);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ExprContext ------------------------------------------------------------------

circuit_specificationParser::ExprContext::ExprContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* circuit_specificationParser::ExprContext::INT() {
  return getToken(circuit_specificationParser::INT, 0);
}

std::vector<circuit_specificationParser::ExprContext *> circuit_specificationParser::ExprContext::expr() {
  return getRuleContexts<circuit_specificationParser::ExprContext>();
}

circuit_specificationParser::ExprContext* circuit_specificationParser::ExprContext::expr(size_t i) {
  return getRuleContext<circuit_specificationParser::ExprContext>(i);
}


size_t circuit_specificationParser::ExprContext::getRuleIndex() const {
  return circuit_specificationParser::RuleExpr;
}

void circuit_specificationParser::ExprContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<circuit_specificationListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterExpr(this);
}

void circuit_specificationParser::ExprContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<circuit_specificationListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitExpr(this);
}


circuit_specificationParser::ExprContext* circuit_specificationParser::expr() {
   return expr(0);
}

circuit_specificationParser::ExprContext* circuit_specificationParser::expr(int precedence) {
  ParserRuleContext *parentContext = _ctx;
  size_t parentState = getState();
  circuit_specificationParser::ExprContext *_localctx = _tracker.createInstance<ExprContext>(_ctx, parentState);
  circuit_specificationParser::ExprContext *previousContext = _localctx;
  size_t startState = 2;
  enterRecursionRule(_localctx, 2, circuit_specificationParser::RuleExpr, precedence);

    size_t _la = 0;

  auto onExit = finally([=] {
    unrollRecursionContexts(parentContext);
  });
  try {
    size_t alt;
    enterOuterAlt(_localctx, 1);
    setState(18);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case circuit_specificationParser::INT: {
        setState(13);
        match(circuit_specificationParser::INT);
        break;
      }

      case circuit_specificationParser::T__4: {
        setState(14);
        match(circuit_specificationParser::T__4);
        setState(15);
        expr(0);
        setState(16);
        match(circuit_specificationParser::T__5);
        break;
      }

    default:
      throw NoViableAltException(this);
    }
    _ctx->stop = _input->LT(-1);
    setState(28);
    _errHandler->sync(this);
    alt = getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 3, _ctx);
    while (alt != 2 && alt != atn::ATN::INVALID_ALT_NUMBER) {
      if (alt == 1) {
        if (!_parseListeners.empty())
          triggerExitRuleEvent();
        previousContext = _localctx;
        setState(26);
        _errHandler->sync(this);
        switch (getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 2, _ctx)) {
        case 1: {
          _localctx = _tracker.createInstance<ExprContext>(parentContext, parentState);
          pushNewRecursionContext(_localctx, startState, RuleExpr);
          setState(20);

          if (!(precpred(_ctx, 4))) throw FailedPredicateException(this, "precpred(_ctx, 4)");
          setState(21);
          _la = _input->LA(1);
          if (!(_la == circuit_specificationParser::T__0

          || _la == circuit_specificationParser::T__1)) {
          _errHandler->recoverInline(this);
          }
          else {
            _errHandler->reportMatch(this);
            consume();
          }
          setState(22);
          expr(5);
          break;
        }

        case 2: {
          _localctx = _tracker.createInstance<ExprContext>(parentContext, parentState);
          pushNewRecursionContext(_localctx, startState, RuleExpr);
          setState(23);

          if (!(precpred(_ctx, 3))) throw FailedPredicateException(this, "precpred(_ctx, 3)");
          setState(24);
          _la = _input->LA(1);
          if (!(_la == circuit_specificationParser::T__2

          || _la == circuit_specificationParser::T__3)) {
          _errHandler->recoverInline(this);
          }
          else {
            _errHandler->reportMatch(this);
            consume();
          }
          setState(25);
          expr(4);
          break;
        }

        } 
      }
      setState(30);
      _errHandler->sync(this);
      alt = getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 3, _ctx);
    }
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }
  return _localctx;
}

bool circuit_specificationParser::sempred(RuleContext *context, size_t ruleIndex, size_t predicateIndex) {
  switch (ruleIndex) {
    case 1: return exprSempred(dynamic_cast<ExprContext *>(context), predicateIndex);

  default:
    break;
  }
  return true;
}

bool circuit_specificationParser::exprSempred(ExprContext *_localctx, size_t predicateIndex) {
  switch (predicateIndex) {
    case 0: return precpred(_ctx, 4);
    case 1: return precpred(_ctx, 3);

  default:
    break;
  }
  return true;
}

// Static vars and initialization.
std::vector<dfa::DFA> circuit_specificationParser::_decisionToDFA;
atn::PredictionContextCache circuit_specificationParser::_sharedContextCache;

// We own the ATN which in turn owns the ATN states.
atn::ATN circuit_specificationParser::_atn;
std::vector<uint16_t> circuit_specificationParser::_serializedATN;

std::vector<std::string> circuit_specificationParser::_ruleNames = {
  "prog", "expr"
};

std::vector<std::string> circuit_specificationParser::_literalNames = {
  "", "'*'", "'/'", "'+'", "'-'", "'('", "')'"
};

std::vector<std::string> circuit_specificationParser::_symbolicNames = {
  "", "", "", "", "", "", "", "NEWLINE", "INT"
};

dfa::Vocabulary circuit_specificationParser::_vocabulary(_literalNames, _symbolicNames);

std::vector<std::string> circuit_specificationParser::_tokenNames;

circuit_specificationParser::Initializer::Initializer() {
	for (size_t i = 0; i < _symbolicNames.size(); ++i) {
		std::string name = _vocabulary.getLiteralName(i);
		if (name.empty()) {
			name = _vocabulary.getSymbolicName(i);
		}

		if (name.empty()) {
			_tokenNames.push_back("<INVALID>");
		} else {
      _tokenNames.push_back(name);
    }
	}

  _serializedATN = {
    0x3, 0x608b, 0xa72a, 0x8133, 0xb9ed, 0x417c, 0x3be7, 0x7786, 0x5964, 
    0x3, 0xa, 0x22, 0x4, 0x2, 0x9, 0x2, 0x4, 0x3, 0x9, 0x3, 0x3, 0x2, 0x3, 
    0x2, 0x3, 0x2, 0x7, 0x2, 0xa, 0xa, 0x2, 0xc, 0x2, 0xe, 0x2, 0xd, 0xb, 
    0x2, 0x3, 0x3, 0x3, 0x3, 0x3, 0x3, 0x3, 0x3, 0x3, 0x3, 0x3, 0x3, 0x5, 
    0x3, 0x15, 0xa, 0x3, 0x3, 0x3, 0x3, 0x3, 0x3, 0x3, 0x3, 0x3, 0x3, 0x3, 
    0x3, 0x3, 0x7, 0x3, 0x1d, 0xa, 0x3, 0xc, 0x3, 0xe, 0x3, 0x20, 0xb, 0x3, 
    0x3, 0x3, 0x2, 0x3, 0x4, 0x4, 0x2, 0x4, 0x2, 0x4, 0x3, 0x2, 0x3, 0x4, 
    0x3, 0x2, 0x5, 0x6, 0x2, 0x23, 0x2, 0xb, 0x3, 0x2, 0x2, 0x2, 0x4, 0x14, 
    0x3, 0x2, 0x2, 0x2, 0x6, 0x7, 0x5, 0x4, 0x3, 0x2, 0x7, 0x8, 0x7, 0x9, 
    0x2, 0x2, 0x8, 0xa, 0x3, 0x2, 0x2, 0x2, 0x9, 0x6, 0x3, 0x2, 0x2, 0x2, 
    0xa, 0xd, 0x3, 0x2, 0x2, 0x2, 0xb, 0x9, 0x3, 0x2, 0x2, 0x2, 0xb, 0xc, 
    0x3, 0x2, 0x2, 0x2, 0xc, 0x3, 0x3, 0x2, 0x2, 0x2, 0xd, 0xb, 0x3, 0x2, 
    0x2, 0x2, 0xe, 0xf, 0x8, 0x3, 0x1, 0x2, 0xf, 0x15, 0x7, 0xa, 0x2, 0x2, 
    0x10, 0x11, 0x7, 0x7, 0x2, 0x2, 0x11, 0x12, 0x5, 0x4, 0x3, 0x2, 0x12, 
    0x13, 0x7, 0x8, 0x2, 0x2, 0x13, 0x15, 0x3, 0x2, 0x2, 0x2, 0x14, 0xe, 
    0x3, 0x2, 0x2, 0x2, 0x14, 0x10, 0x3, 0x2, 0x2, 0x2, 0x15, 0x1e, 0x3, 
    0x2, 0x2, 0x2, 0x16, 0x17, 0xc, 0x6, 0x2, 0x2, 0x17, 0x18, 0x9, 0x2, 
    0x2, 0x2, 0x18, 0x1d, 0x5, 0x4, 0x3, 0x7, 0x19, 0x1a, 0xc, 0x5, 0x2, 
    0x2, 0x1a, 0x1b, 0x9, 0x3, 0x2, 0x2, 0x1b, 0x1d, 0x5, 0x4, 0x3, 0x6, 
    0x1c, 0x16, 0x3, 0x2, 0x2, 0x2, 0x1c, 0x19, 0x3, 0x2, 0x2, 0x2, 0x1d, 
    0x20, 0x3, 0x2, 0x2, 0x2, 0x1e, 0x1c, 0x3, 0x2, 0x2, 0x2, 0x1e, 0x1f, 
    0x3, 0x2, 0x2, 0x2, 0x1f, 0x5, 0x3, 0x2, 0x2, 0x2, 0x20, 0x1e, 0x3, 
    0x2, 0x2, 0x2, 0x6, 0xb, 0x14, 0x1c, 0x1e, 
  };

  atn::ATNDeserializer deserializer;
  _atn = deserializer.deserialize(_serializedATN);

  size_t count = _atn.getNumberOfDecisions();
  _decisionToDFA.reserve(count);
  for (size_t i = 0; i < count; i++) { 
    _decisionToDFA.emplace_back(_atn.getDecisionState(i), i);
  }
}

circuit_specificationParser::Initializer circuit_specificationParser::_init;
