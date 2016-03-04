%option reentrant
%option bison-bridge
%option noyywrap
%option nounput
%option noinput

%{

#define YYSTYPE mps_formal_polynomial *
#define YY_EXTRA_TYPE _mps_yacc_parser_data*

#include <stdio.h>
#include <mps/mps.h>
#include "yacc-parser.h"

#ifndef _MPS_YACC_PARSER_H
#define _MPS_YACC_PARSER_H
typedef struct {
  void * scanner;
  mps_context * ctx;
  mps_abstract_input_stream * stream;
  mps_formal_polynomial * p;
} _mps_yacc_parser_data;
#endif

void yyerror(yyscan_t scanner, _mps_yacc_parser_data * data, const char * str)
{
  mps_context * ctx = data->ctx;
  mps_error (ctx, str);
}

extern int yylex_init_extra (YY_EXTRA_TYPE, yyscan_t *);
extern int yylex_destroy (yyscan_t);
extern int mps_abstract_input_stream_getchar (mps_abstract_input_stream *);

#define YY_INPUT(buf, result, max_size) { \
  mps_abstract_input_stream * stream = yyextra->stream; \
  int c = mps_abstract_input_stream_getchar(stream); \
  result = (c == EOF) ? YY_NULL : (buf[0] = c, 1); \
}

YY_EXTRA_TYPE yyget_extra (yyscan_t yyscanner);
YYSTYPE * yyget_lval (yyscan_t yyscanner);

mps_polynomial *
mps_monomial_yacc_parser (mps_context * ctx, mps_abstract_input_stream * stream)
{
  _mps_yacc_parser_data data = { NULL, ctx, stream, NULL };
  yylex_init_extra (&data, &data.scanner);
  yyparse (data.scanner, &data);
  yylex_destroy (data.scanner);

  if (! mps_context_has_errors (ctx))
    return MPS_POLYNOMIAL (mps_formal_polynomial_create_monomial_poly (data.p, ctx));
  else
    return NULL;
}

%}

%%

[0-9]+(\/[0-9]+)?		*yylval = (mps_formal_polynomial*) strdup (yytext); return RATIONAL; 
[0-9]+(.[0-9]*)?([eE][0-9]+)?	*yylval = (mps_formal_polynomial*) strdup (yytext); return FLOATING_POINT;
\*?[a-zA-Z](^[0-9]+)?		*yylval = (mps_formal_polynomial*) strdup (yytext); return MONOMIAL;
\+				return PLUS;
\-				return MINUS;
[ \t]+                          ;

%%

