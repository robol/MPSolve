/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2014, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#include <mps/mpsxx.h>

using namespace mps;

extern "C"
{
  mps_file_input_stream * mps_file_input_stream_new (FILE * source)
  {
    return reinterpret_cast<mps_file_input_stream*>( new FileInputStream(source) );
  }

  size_t mps_file_input_stream_readline (mps_file_input_stream * stream, 
					 char ** buffer, size_t * length)
  {
    return reinterpret_cast<FileInputStream*>(stream)->readline(buffer, length);
  }

  mps_boolean mps_file_input_stream_eof (mps_file_input_stream * stream)
  {
    return reinterpret_cast<FileInputStream*>(stream)->eof();
  }

  void 
  mps_file_input_stream_free (mps_file_input_stream * stream)
  {
    delete reinterpret_cast<FileInputStream*> (stream);
  }
}

FileInputStream::FileInputStream (FILE * source)
{
  mSource = source;
}

size_t
FileInputStream::readline (char ** buffer, size_t * length)
{
  return getline (buffer, length, mSource);
}

bool
FileInputStream::eof ()
{
  return feof (mSource);
}

FileInputStream::~FileInputStream()
{
}
