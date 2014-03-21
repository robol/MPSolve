/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2014, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

/**
 * @file
 * @brief Generic parsers for common polynomial types.
 */

#ifndef MPS_PARSER_H_
#define MPS_PARSER_H_

MPS_BEGIN_DECLS

mps_polynomial * mps_parse_stream (mps_context * s, FILE * input_stream);
mps_polynomial * mps_parse_file (mps_context * s, const char * path);
mps_polynomial * mps_parse_string (mps_context * s, const char * c_string);

mps_polynomial * mps_parse_inline_poly (mps_context * ctx, FILE * stream);
mps_polynomial * mps_parse_inline_poly_from_string (mps_context * ctx, const char * input);


MPS_END_DECLS

#endif /* MPS_PARSER_H_ */

