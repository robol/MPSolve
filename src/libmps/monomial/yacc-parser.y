%pure-parser
%parse-param { void * scanner }
%lex-param { void * scanner }
%param { void * data }

%{

#include <mps/mps.h> 
#include <stdio.h>
#include <string.h>

  /*
   * Here we are using a quite hackish way to transfer monomials in the different parts
   * of the parsing. 
   *
   * We use char* as YYSTYPE for the lex tokenizer and here we use a custom struct
   * to better hold coefficient and exponent of every monomial. 
   *
   * So we have an implicit casting that makes the char* passed  from the tokenizer like 
   * the first field in the struct, i.e., the coefficient. 
   *
   */

  #define YYSTYPE mps_formal_polynomial*
  #define YY_EXTRA_TYPE _mps_yacc_parser_data*
  
  #ifndef _MPS_YACC_PARSER_H
  #define _MPS_YACC_PARSER_H
  typedef struct {
    void * scanner;
    mps_context * ctx;
    mps_abstract_input_stream * stream;
    mps_formal_polynomial * p;
  } _mps_yacc_parser_data;
  #endif

  extern int yylex(void * yylval, void * scanner, void * data);
  extern int yyerror(void*,void*,const char*);

%}

%token RATIONAL FLOATING_POINT MONOMIAL PLUS MINUS

%%

polynomial: monomial
	    {
	      $$ = mps_formal_polynomial_new_with_monomial ((mps_formal_monomial*) $1);
	      mps_formal_monomial_free ((mps_formal_monomial *) $1);
	      ((_mps_yacc_parser_data *) data)->p = $$;

#ifdef MPS_PARSER_DEBUG
	      printf ("Created polynomial: "); mps_formal_polynomial_print ($$); printf ("\n");
#endif
	    }
	    | polynomial PLUS monomial
	    {
#ifdef MPS_PARSER_DEBUG
	      printf ("Extending polynomial: "); mps_formal_polynomial_print ($1);
	      printf (" with "); mps_formal_monomial_print ((mps_formal_monomial*) $3); printf ("\n");
#endif

	      $$ = mps_formal_polynomial_sum_eq ($1, (mps_formal_monomial*) $3);
	      mps_formal_monomial_free ((mps_formal_monomial *) $3);

#ifdef MPS_PARSER_DEBUG
	      printf ("Extended polynomial: "); mps_formal_polynomial_print ($$); printf ("\n");
	      printf ("Polynomial: %p\n", $$);
#endif
	    }
	    | polynomial MINUS monomial
	    {
#ifdef MPS_PARSER_DEBUG
	      printf ("Extending polynomial: "); mps_formal_polynomial_print ($1); printf (" with -"); 
	      mps_formal_monomial_print ((mps_formal_monomial*) $3); printf ("\n");
#endif

	      $$ = mps_formal_polynomial_sub_eq ($1, (mps_formal_monomial*) $3);
	      mps_formal_monomial_free ((mps_formal_monomial *) $3);

#ifdef MPS_PARSER_DEBUG
	      printf ("Extended polynomial: "); mps_formal_polynomial_print ($$); printf ("\n");	      
	      printf ("Polynomial: %p\n", $$);
#endif
	    }

monomial: MINUS monomial
            {
              $$ = (mps_formal_polynomial *) mps_formal_monomial_neg ((mps_formal_monomial *) $2);
#ifdef MPS_PARSER_DEBUG
	      printf ("Changed sign to monomial: "); 
	      mps_formal_monomial_print ((mps_formal_monomial*) $$); 
	      printf ("\n");
#endif
              mps_formal_monomial_free ((mps_formal_monomial*) $2);
            }

monomial: RATIONAL
	  {
#ifdef MPS_PARSER_DEBUG
	    printf ("Rational coefficient: %s\n", (const char *) $1);
#endif
	    $$ = (mps_formal_polynomial*) mps_formal_monomial_new_with_string ((const char*) $1, 0);	    
	    free ($1);
	  }		
	  | FLOATING_POINT 
	  {
#ifdef MPS_PARSER_DEBUG
	    printf ("Floating point coefficient: %s\n", (const char *) $1);
#endif
	    mps_formal_monomial * m = mps_formal_monomial_new_with_string ((const char *) $1, 0);
	    free ($1);
	    $$ = (mps_formal_polynomial *) m; 
	  }
	  | MONOMIAL 
	  {
#ifdef MPS_PARSER_DEBUG
 	    printf ("Simple monomial: %s\n", (const char *) $1);
#endif

	    const char * exp = strchr ((const char *) $1, '^');
	    long degree = (exp == NULL) ? 1 : atoi (exp + 1);
	    $$ = (mps_formal_polynomial*) mps_formal_monomial_new_with_string ("1", degree);
	    free ($1);
	  }
	  | RATIONAL MONOMIAL
	  {
#ifdef MPS_PARSER_DEBUG
	    printf ("Coefficient: %s %s\n", (const char*) $1, (const char*) $2);
#endif

	    const char * exp = strchr ((const char *) $2, '^');
	    long degree = (exp == NULL) ? 1 : atoi (exp + 1);
	    $$ = (mps_formal_polynomial*) mps_formal_monomial_new_with_string ((const char *) $1, degree);

	    free ($1);
	    free ($2);
	  }
	  | FLOATING_POINT MONOMIAL
	  {
#ifdef MPS_PARSER_DEBUG
	    printf ("Coefficient: %s %s\n", (const char *) $1, (const char*) $2);
#endif

	    const char * exp = strchr ((const char *) $2, '^');
	    long degree = (exp == NULL) ? 1 : atoi (exp + 1);
	    $$ = (mps_formal_polynomial*) mps_formal_monomial_new_with_string ((const char *) $1, degree);

	    free ($1);
	    free ($2);
	  }


%%
