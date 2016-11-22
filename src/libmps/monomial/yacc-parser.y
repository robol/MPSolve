%pure-parser
%parse-param { void * scanner }
%parse-param { void * data }
%lex-param { void * scanner }
%lex-param { void * data }
%debug

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

%token RATIONAL FLOATING_POINT PLUS MINUS IMAGINARY_UNIT TIMES LEFT_BRACKET RIGHT_BRACKET MONOMIAL

%left PLUS MINUS
%left TIMES
%right IMAGINARY_UNIT

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
            | LEFT_BRACKET polynomial RIGHT_BRACKET
	    {
	      $$ = $2;
	    }
            | polynomial TIMES polynomial
	    {
#ifdef MPS_PARSER_DEBUG
	      printf ("Multiplying polynomials "); mps_formal_polynomial_print ($1);
	      printf (" and "); mps_formal_polynomial_print ($3); printf ("\n");
#endif
	      $$ = mps_formal_polynomial_mul_eq ($1, $3);
	      mps_formal_polynomial_free ($3);
	      
	      ((_mps_yacc_parser_data *) data)->p = $$;
	      
#ifdef MPS_PARSER_DEBUG
	      printf ("Multiplied polynomials:"); mps_formal_polynomial_print ($$);
	      printf ("\n");
#endif	
	    }
            | polynomial PLUS polynomial
            {
#ifdef MPS_PARSER_DEBUG
	      printf ("Summing polynomials: "); mps_formal_polynomial_print ($1);
	      printf (" and ");
	      mps_formal_polynomial_print($3);
	      printf ("\n");
#endif
              $$ = mps_formal_polynomial_sum_eq_p ($1, $3);
	      mps_formal_polynomial_free ($3);
	      ((_mps_yacc_parser_data *) data)->p = $$;
            }
            | polynomial MINUS polynomial
	    {
              $$ = mps_formal_polynomial_sub_eq_p ($1, $3);
	      mps_formal_polynomial_free ($3);
	      ((_mps_yacc_parser_data *) data)->p = $$;
	    }

monomial: MONOMIAL 
	  {
#ifdef MPS_PARSER_DEBUG
 	    printf ("Simple monomial: %s\n", (const char *) $1);
#endif

	    const char * exp = strchr ((const char *) $1, '^');
	    long degree = (exp == NULL) ? 1 : atoi (exp + 1);
	    $$ = (mps_formal_polynomial*) mps_formal_monomial_new_with_string ("1", degree);
	    free ($1);
	  }
	  | number MONOMIAL
	  {
#ifdef MPS_PARSER_DEBUG
	    printf ("Number: "); mps_formal_monomial_print ((mps_formal_monomial *) $1); 
	    printf (" Monomial: %s\n", (const char *) $2); 
#endif
	    const char * exp = strchr ((const char *) $2, '^');
	    long degree = (exp == NULL) ? 1 : atoi (exp + 1);
	    mps_formal_monomial * m = mps_formal_monomial_new_with_string ("1", degree);

	    $$ = (mps_formal_polynomial*) mps_formal_monomial_mul_eq ((mps_formal_monomial *) $1, m);
	    mps_formal_monomial_free ((mps_formal_monomial *) m);
	  }
          | number 
	  {
	    $$ = $1;
	  }
          | MINUS monomial
          {
	    $$ = (mps_formal_polynomial *) mps_formal_monomial_neg ((mps_formal_monomial *) $2);
#ifdef MPS_PARSER_DEBUG
	    printf ("Changed sign to monomial: "); 
	    mps_formal_monomial_print ((mps_formal_monomial*) $$); 
	    printf ("\n");
#endif
	    mps_formal_monomial_free ((mps_formal_monomial*) $2);
	  }

number: RATIONAL
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
          | number IMAGINARY_UNIT
          {
            #ifdef MPS_PARSER_DEBUG
	    printf ("Imaginary rule\n");
	    printf ("Number: "); mps_formal_monomial_print ((mps_formal_monomial *) $1); printf (" i\n");
            #endif

	    mps_formal_monomial * imunit = mps_formal_monomial_new_with_strings ("0", "1", 0);
            $$ = (mps_formal_polynomial *) mps_formal_monomial_mul_eq (imunit, 
								       (mps_formal_monomial*) $1);
	    mps_formal_monomial_free ((mps_formal_monomial*) $1);
          }


%%
