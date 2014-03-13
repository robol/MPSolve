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
#include <cstring>
#include <mps/private/system/memory-file-stream.h>

using namespace mps;

extern "C"
{
  mps_memory_file_stream * 
  mps_memory_file_stream_new (const char * source)
  {
    return reinterpret_cast<mps_memory_file_stream*> (new MemoryFileStream(source));
  }

  void
  mps_memory_file_stream_free (mps_memory_file_stream * stream)
  {
    delete reinterpret_cast<MemoryFileStream*> (stream);
  }
}

#define MPS_MAXIMUM_LINE_LENGTH (1024 * 1024)

MemoryFileStream::MemoryFileStream(const char * source) : 
  mInputStream (source)
{
}

MemoryFileStream::~MemoryFileStream()
{
}

size_t
MemoryFileStream::readline(char ** buffer, size_t * length)
{
  if (*buffer == NULL)
    {
      *buffer = mps_newv (char, 1024); 
      *length = 1024;
    }

  mInputStream.getline(*buffer, *length);

  /* Try to increase the size of the buffer until we can read the string. */
  while ( (mInputStream.fail() && ! (mInputStream.eof() || mInputStream.bad())) && 
	  (*length < MPS_MAXIMUM_LINE_LENGTH))
    {
      *length *= 2;
      *buffer = (char*) mps_realloc (*buffer, sizeof (char) * *length);

      mInputStream.getline (*buffer, *length);
    }

  return (mInputStream.fail()) ? -1 : strlen (*buffer) + 1;
}

bool
MemoryFileStream::eof()
{
  return mInputStream.eof();
}


