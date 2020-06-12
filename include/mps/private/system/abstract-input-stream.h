/*
 * This file is part of MPSolve 3.2.1
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

/**
 * @file
 * @brief Abstract input stream used to parse files.
 *
 * The implementation is C++ but C API is provided.
 */

#ifndef MPS_ABSTRACT_INPUT_STREAM_H_
#define MPS_ABSTRACT_INPUT_STREAM_H_

/* C compatibility layer */
MPS_BEGIN_DECLS
/**
 * @brief This is a C wrapper around the C++ implementation of
 * {@link AbstractInputStream}.
 */
struct mps_abstract_input_stream;

/**
 * @brief C wrapper around {@link AbstractInputStream}.
 */
typedef struct mps_abstract_input_stream mps_abstract_input_stream;

/**
 * @brief Wrapper around {@link AbstractInputStream::readline()}.
 */
size_t mps_abstract_input_stream_readline (mps_abstract_input_stream * stream,
                                           char ** buffer, size_t * length);

/**
 * @brief Wrapper around {@link AbstractInputStream::eof()}.
 */
mps_boolean mps_abstract_input_stream_eof (mps_abstract_input_stream * stream);

/**
 * @brief Wrapper around {@link AbstractInputStream::getchar()}. 
 */
int mps_abstract_input_stream_getchar (mps_abstract_input_stream * stream);

MPS_END_DECLS

/* The following is C++ only */
#ifdef __cplusplus

namespace mps {
  /**
   * @brief Abstract class that represent a generic input stream that can
   * be used by MPSolve to read polynomial files and/or descriptions.
   *
   * @seealso {@link MemoryFileStream}, {@link FileInputStream}
   */
  class AbstractInputStream {
public:

    virtual ~AbstractInputStream() = 0;

    /**
     * @brief Return a new line of the stream or NULL if we are at
     * the end.
     *
     * @return A pointer to a newly allocated line or NULL.
     */
    virtual size_t readline (char ** buffer, size_t * length) = 0;

    /**
     * @brief Check if we are at the end of the stream.
     *
     * @return true if we are at the end of the stream.
     */
    virtual bool eof () = 0;

    /**
     * @brief Obtain a single character. 
     *
     * @return A new character read from the stream. 
     */
    virtual int getchar () = 0;
  };
}

#endif /* __cplusplus */

#endif /* MPS_ABSTRACT_INPUT_STREAM_H_ */
