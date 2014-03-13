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
  * @brief
  */

#ifndef MPS_FILE_INPUT_STREAM_H_
#define MPS_FILE_INPUT_STREAM_H_

MPS_BEGIN_DECLS

struct mps_file_input_stream;
typedef struct mps_file_input_stream mps_file_input_stream;

mps_file_input_stream * mps_file_input_stream_new (FILE * source);
size_t mps_file_input_stream_readline (mps_file_input_stream * stream, 
				       char ** buffer, size_t * length);
mps_boolean mps_file_input_stream_eof (mps_file_input_stream * stream);
void mps_file_input_stream_free (mps_file_input_stream * stream);

MPS_END_DECLS

#ifdef __cplusplus

namespace mps {

  class FileInputStream : AbstractInputStream {

  public:
    
    FileInputStream (FILE * source);

    ~FileInputStream ();
    
    size_t readline (char ** buffer, size_t * length);

    bool eof (); 

  private:
    FILE * mSource;

  };

}

#endif /* __cplusplus */

#endif /* MPS_FILE_INPUT_STREAM_H_ */

