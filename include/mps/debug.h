/* 
 * File:   debug.h
 * Author: leonardo
 *
 * Created on 23 aprile 2011, 11.00
 */

/**
 * @file
 * @brief Debugging functions, that honor the status os <code>s->DOLOG</code> and 
 * autodetect if the output stream is or not a tty to select the best output method.
 */

#ifndef MPS_DEBUG_H
#define  MPS_DEBUG_H

#ifndef __WINDOWS
#include <unistd.h>
#else
#include <io.h>
#endif

#include <mps/gmptools.h>
#include <gmp.h>
#include <time.h>

#ifdef  __cplusplus  
extern "C"  
{  
#endif  

  /* Timer functions */
  clock_t * mps_start_timer (void);
  unsigned long int mps_stop_timer (clock_t * my_timer);

  /* Keep away assert() when compiling without debug */
#ifdef DISABLE_DEBUG
#define NDEBUG
#endif

#ifdef NICE_DEBUG

  /**
   * @brief Shorthand for compiling with the MPS_DEBUG_INFO level.
   */
#define MPS_DEBUG_WITH_INFO(s, templ...) {      \
    if (s->debug_level & MPS_DEBUG_INFO) {      \
      MPS_DEBUG(s, templ);                      \
    }                                           \
}

#define MPS_DEBUG(s, templ...) {                \
    __MPS_DEBUG(s, templ) ;                     \
    if (s->DOLOG) {                             \
      fprintf(s->logstr, "\n");                 \
    }                                           \
  }

  /**
   * @brief Debug the value of a complex multiprecision
   * variable.
   */
#define MPS_DEBUG_MPC(s, digits, c, name...) {  \
    __MPS_DEBUG_EQ(s, name);                    \
    if (s->DOLOG) {                             \
      mpc_outln_str(s->logstr, 10, digits, c);  \
    }                                           \
  }

  /**
   * @brief Debug the value of a real multiprecision
   * variable.
   */
#define MPS_DEBUG_MPF(s, digits, c, name...) {  \
    __MPS_DEBUG_EQ(s, name);                    \
    if (s->DOLOG) {                             \
      mpf_out_str(s->logstr, 10, digits, c);    \
      fprintf (s->logstr, "\n");                \
    }                                           \
  }

#define MPS_DEBUG_MPC2(s, radius, c, name...) {         \
  __MPS_DEBUG_EQ(s, name);                              \
  if (s->DOLOG) {                                       \
  int t = -rdpe_log10 (radius) - 1;                     \
  cdpe_t ctmp;                                          \
  rdpe_t mm;                                            \
  mpc_get_cdpe (ctmp, c);                               \
  cdpe_mod (mm, ctmp);                                  \
  t += rdpe_log (mm);                                   \
  mpc_outln_str (s->logstr, 10, t, c);                  \
  }                                                     \
}
  

  /**
   * @brief Debug the value of a rdpe variable.
   */
#define MPS_DEBUG_RDPE(s, r, name...) {         \
    __MPS_DEBUG_EQ(s, name);                    \
    if (s->DOLOG) {                             \
      rdpe_outln_str(s->logstr, r);             \
    }                                           \
  }

  /**
   * @brief Debug the value of a cdpe variable.
   */
#define MPS_DEBUG_CDPE(s, c, name...) {         \
    __MPS_DEBUG_EQ(s, name);                    \
    if (s->DOLOG) {                             \
      cdpe_outln_str(s->logstr, c);             \
    }                                           \
  }

  /**
   * @brief Debug the values of a cplx_t variable
   */
#define MPS_DEBUG_CPLX(s, c, name...) {         \
    __MPS_DEBUG_EQ(s, name);                    \
    if (s->DOLOG) {                             \
      cplx_outln_str(s->logstr, c);             \
    }                                           \
  }

  /**
   * @brief Make some space in the debug stream to make clean that
   * another section is starting.
   */
#define MPS_DEBUG_BREAK(s) if (s->DOLOG) {      \
    fprintf(s->logstr, "\n");                   \
  }

  /**
   * @brief Low-level debug print that appends an " = " sign to the
   * output (useful for debugging values of variables).
   */
#define __MPS_DEBUG_EQ(s, templ...)             \
  __MPS_DEBUG(s, templ);                        \
  if (s->DOLOG) {                               \
    fprintf(s->logstr, " = ");                  \
  }

  /**
   * @brief Debug that a function is going to be called.
   */
#ifndef __WINDOWS
#define MPS_DEBUG_CALL(s, function) if (s->DOLOG && (s->debug_level & MPS_DEBUG_FUNCTION_CALLS)) { \
    if (isatty(fileno(s->logstr))) {                                   \
      __MPS_DEBUG(s, "Called \033[31;1m");                              \
    }                                                                   \
    else {                                                              \
      __MPS_DEBUG(s, "Called ");                                        \
    }                                                                   \
    if (isatty(fileno(s->logstr))) {                                   \
      fprintf(s->logstr, function); fprintf(s->logstr, "()\033[0m\n");  \
    }                                                                   \
    else                                                                \
      {                                                                 \
        fprintf(s->logstr, function); fprintf(s->logstr, "()\n");       \
      }                                                                 \
}
#else
#define MPS_DEBUG_CALL(s, function) if (s->DOLOG && (s->debug_level & MPS_DEBUG_FUNCTION_CALLS)) { \
    if (_isatty(_fileno(s->logstr))) {                                  \
      __MPS_DEBUG(s, "Called \033[31;1m");                              \
    }                                                                   \
    else {                                                              \
      __MPS_DEBUG(s, "Called ");                                        \
    }                                                                   \
    if (_isatty(_fileno(s->logstr))) {                                  \
      fprintf(s->logstr, function); fprintf(s->logstr, "()\033[0m\n");  \
    }                                                                   \
    else                                                                \
      {                                                                 \
        fprintf(s->logstr, function); fprintf(s->logstr, "()\n");       \
      }                                                                 \
  }
#endif

#define MPS_DEBUG_THIS_CALL MPS_DEBUG_CALL(s, __FUNCTION__)


  /**
   * @brief Low-level DEBUG() used by other MPS_DEBUG_* statements.
   */
#if __STDC_VERSION__ >= 199901L
#ifndef __WINDOWS
#define __MPS_DEBUG(s, templ...) if (s->DOLOG) {                \
    if (isatty(fileno(s->logstr))) {                           \
      fprintf(s->logstr, "%s:%d \033[32;1m%s()\033[;0m ",       \
              __FILE__, __LINE__, __FUNCTION__);                \
    }                                                           \
    else {                                                      \
      fprintf(s->logstr, "%s:%d %s() ",                         \
              __FILE__, __LINE__, __FUNCTION__);                \
    }                                                           \
    fprintf(s->logstr, templ);                                  \
  }
#else
#define __MPS_DEBUG(s, templ...) if (s->DOLOG) {                \
    if (_isatty(_fileno(s->logstr))) {                          \
      fprintf(s->logstr, "%s:%d \033[32;1m%s()\033[;0m ",       \
              __FILE__, __LINE__, __FUNCTION__);                \
    }                                                           \
    else {                                                      \
      fprintf(s->logstr, "%s:%d %s() ",                         \
              __FILE__, __LINE__, __FUNCTION__);                \
    }                                                           \
    fprintf(s->logstr, templ);                                  \
  }
#endif

#endif
#else
#define MPS_DEBUG(args...)
#define MPS_DEBUG_MPC(args...)
#define MPS_DEBUG_CPLX(args...)
#define MPS_DEBUG_RDPE(args...)
#define MPS_DEBUG_CDPE(args...)
#define MPS_DEBUG_CALL(args...)
#define MPS_DEBUG_MCLUSTER_ROOTS(args...)
#define MPS_DEBUG_THIS_CALL
#define MPS_DEBUG_WITH_INFO(args...)
#define __MPS_DEBUG(args...)
#endif

  /*
   * Ansi C version implemented without using variadic macros, only for compatibility
   * since could became imprecise when using complex debug instructions.
   */
#ifndef DISABLE_DEBUG
#if __STDC_VERSION__ < 199901L
#include <mps/mps.h>
#ifdef MPS_DEBUG
#undef MPS_DEBUG
#endif
#ifdef __MPS_DEBUG
#undef __MPS_DEBUG
#endif
#define MPS_DEBUG __c_impl__MPS_DEBUG
#define __MPS_DEBUG __c_impl____MPS_DEBUG
  void __c_impl__MPS_DEBUG (mps_context * s, const char *templ, ...);
  void __c_impl____MPS_DEBUG (mps_context * s, const char *templ, ...);
#endif
#endif

  /**
   * @brief This is the flag that enables debugging of general
   * information about the flow of the program.
   */
#define MPS_DEBUG_INFO     (0x0001 << 0)

  /**
   * @brief This is the flag used to debug cluster-related
   * informations.
   */
#define MPS_DEBUG_CLUSTER  (0x0001 << 1)

  /**
   * @brief This is the flag used to debug informations on
   * approximations during the execution of MPSolve.
   */
#define MPS_DEBUG_APPROXIMATIONS (0x0001 << 2)

  /**
   * @brief Flag used to show debug about final approximation
   * of the roots using newton iterations.
   */
#define MPS_DEBUG_IMPROVEMENT (0x0001 << 3)

  /**
   * @brief Flag used to obtain information about timings in 
   * the algorithm.
   */
#define MPS_DEBUG_TIMINGS (0x0001 << 4)

  /**
   * @brief Debug function calls
   */
#define MPS_DEBUG_FUNCTION_CALLS (0x0001 << 5)

  /**
   * @brief Debug I/O informations
   */
#define MPS_DEBUG_IO (0x0001 << 6)

  /**
   * @brief Debug memory management
   */
#define MPS_DEBUG_MEMORY (0x0001 << 7)

  /**
   * @brief Debug checks for the convergence
   * of the various iterations packets.
   */
#define MPS_DEBUG_PACKETS (0x0001 << 8)

  /**
   * @brief Debug the regenration of the coefficients-
   */
#define MPS_DEBUG_REGENERATION (0x0001 << 9)

  /**
   * @brief This is the flag used to debug informations
   * about virtually every step in the program, it enables
   * every debug level.
   */
#define MPS_DEBUG_TRACE    (0xFFFF)

#define MPS_DEBUG_IF(s, debug_level, debug_instruction) if (debug_level & s->debug_level) { \
    debug_instruction;                                                  \
  }

#ifdef  __cplusplus  
} 
#endif  

#endif                          /* DEBUG_H */
