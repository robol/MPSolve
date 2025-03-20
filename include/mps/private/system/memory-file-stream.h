/*
 * This file is part of MPSolve 3.2.2
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

/**
 * @file
 * @brief Implementation of a fmemopen-like stream.
 */

#ifndef MPS_MEMORY_FILE_STREAM_H_
#define MPS_MEMORY_FILE_STREAM_H_

MPS_BEGIN_DECLS

/**
 * @brief C wrapper around {@link MemoryFileStream}.
 */
struct mps_memory_file_stream;

/**
 * @brief C wrapper around {@link MemoryFileStream}.
 */
typedef struct mps_memory_file_stream mps_memory_file_stream;

/**
 * @brief Allocate a new {@link MemoryFileStream} that will output
 * the data pointed by source.
 *
 * @param source The data that will be provided by the mps_memory_file_stream.
 */
mps_memory_file_stream * mps_memory_file_stream_new (const char * source);

/**
 * @brief Release the resources holded by a {@link MemoryFileStream}.
 *
 * @param stream The {@link MemoryFileStream} to release.
 */
void mps_memory_file_stream_free (mps_memory_file_stream * stream);

MPS_END_DECLS

#ifdef __cplusplus

#include <iostream>
#include <sstream>

#include <mps/mps.h>

namespace mps {
  /**
   * @brief The MemoryFileStream class provides an implementation of
   * the abstract class {@link AbstractInputStream} that will stream
   * the data contained in the area stored in memory.
   */
  class MemoryFileStream : public AbstractInputStream {
public:

    /**
     * @brief Allocate a new MemoryFileStream that wil provide that data
     * stored by the given pointer.
     *
     * @param source A pointer to the data that should be provided
     * by this instance.
     */
    MemoryFileStream(const char * source);

    ~MemoryFileStream();

    /**
     * @brief Implementation of the readline() method of the
     * {@link AbstractInputStream} parent.
     *
     * @param buffer A pointer to the buffer where the line will be stored.
     * @param length A pointer where the length of the allocated buffer at
     * the end will be saved.
     *
     * @return The number of characters that have been stored in
     * buffer.
     */
    size_t readline (char ** buffer, size_t * length);

    /**
     * @brief Implementation of the eof() method of {@link AbstractInputStream}.
     *
     * @return true if the source stream has reached the end.
     */
    bool eof ();

    /**
     * @brief Obtain a single character. 
     *
     * @return A new character read from the stream. 
     */
    int getchar ();

private:
    std::istringstream mInputStream;
  };
}

#endif

#endif /* MPS_MEMORY_FILE_STREAM_H_ */

