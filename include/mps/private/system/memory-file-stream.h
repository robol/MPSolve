/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2014, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

 /**
  * @file
  * @brief Implementation of a fmemopen-like stream. 
  */

#ifndef MPS_MEMORY_FILE_STREAM_H_
#define MPS_MEMORY_FILE_STREAM_H_

MPS_BEGIN_DECLS

  struct mps_memory_file_stream;
  typedef struct mps_memory_file_stream mps_memory_file_stream;

  mps_memory_file_stream * mps_memory_file_stream_new (const char * source);
  void mps_memory_file_stream_free (mps_memory_file_stream * stream);

MPS_END_DECLS

#ifdef __cplusplus

#include <iostream>
#include <sstream>

#include <mps/mps.h>

namespace mps {

  class MemoryFileStream : AbstractInputStream {

  public: 

    MemoryFileStream(const char * source);
    ~MemoryFileStream();

    size_t readline(char ** buffer, size_t * length);

    bool eof ();

  private:
    std::istringstream mInputStream;

  };

}

#endif

#endif /* MPS_MEMORY_FILE_STREAM_H_ */

