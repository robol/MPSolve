/*
 * This file is part of MPSolve 3.1.7
 *
 * Copyright (C) 2001-2014, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#include <mps/mps.h>

using namespace mps;

extern "C" 
{
  size_t mps_abstract_input_stream_readline (mps_abstract_input_stream * stream, 
					     char ** buffer, size_t * length)
  {
    return reinterpret_cast<AbstractInputStream*>(stream)->readline (buffer, length);
  }

  mps_boolean mps_abstract_input_stream_eof (mps_abstract_input_stream * stream)
  {
    mps_boolean eof = reinterpret_cast<AbstractInputStream*>(stream)->eof (); 
    return eof;
  }

  int mps_abstract_input_stream_getchar (mps_abstract_input_stream * stream)
  {
    return reinterpret_cast<AbstractInputStream*> (stream)->getchar();
  }
}

AbstractInputStream::~AbstractInputStream()
{
}
