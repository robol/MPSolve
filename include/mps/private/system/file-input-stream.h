/*
 * This file is part of MPSolve 3.2.0
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

/**
 * @file
 * @brief
 */

#ifndef MPS_FILE_INPUT_STREAM_H_
#define MPS_FILE_INPUT_STREAM_H_

MPS_BEGIN_DECLS

/**
 * @brief Wrapper around {@link FileInputStream}.
 */
struct mps_file_input_stream;

/**
 * @brief Wrapper around {@link FileInputStream}.
 */
typedef struct mps_file_input_stream mps_file_input_stream;

/**
 * @brief Allocate a new {@link FileInputStream} instane that will
 * stream the given file.
 *
 * @param source A FILE* object returned by a call to fopen() on the
 * desired file.
 */
mps_file_input_stream * mps_file_input_stream_new (FILE * source);

/**
 * @brief Release the resources holded by this {@link FileInputStream}
 * instance.
 *
 * @param stream The {@link FileInputStream} that should be freed.
 */
void mps_file_input_stream_free (mps_file_input_stream * stream);

MPS_END_DECLS

#ifdef __cplusplus

namespace mps {
  class FileInputStream : public AbstractInputStream {
public:

    /**
     * @brief Create a new instance of FileInputStream that will output the
     * content of the FILE opened.
     *
     * @param source A FILE* object returned by a fopen() call on the
     * file whose content should be provided by this stream.
     */
    FileInputStream (FILE * source);

    ~FileInputStream ();

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
    FILE * mSource;
  };
}

#endif /* __cplusplus */

#endif /* MPS_FILE_INPUT_STREAM_H_ */

