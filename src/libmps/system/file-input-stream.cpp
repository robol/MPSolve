/*
 * This file is part of MPSolve 3.1.8
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
  mps_file_input_stream * mps_file_input_stream_new (FILE * source)
  {
    return reinterpret_cast<mps_file_input_stream*>( new FileInputStream(source) );
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

int
FileInputStream::getchar ()
{
  return fgetc (mSource);
}

FileInputStream::~FileInputStream()
{
}
