/* 
 * File:   debug.h
 * Author: leonardo
 *
 * Created on 23 aprile 2011, 11.00
 */

#ifndef DEBUG_H
#define	 DEBUG_H

#ifdef	__cplusplus
extern "C" {
#endif

#ifndef __WINDOWS
#include <unistd.h>
#else
#include <io.h>
#endif

#include <mps/gmptools.h>

/* Keep away assert() when compiling without debug */
#ifdef DISABLE_DEBUG
#define NDEBUG
#endif

#ifdef NICE_DEBUG
/**
 * @brief Debug cluster approximations of the roots in the
 * case of multiprecision computation.
 */
#define MPS_DEBUG_MCLUSTER_ROOTS(s, i) if (s->DOLOG) { \
  int __debug_i; \
  __MPS_DEBUG(s, "Dumping cluster approximations:\n"); \
  for(__debug_i = s->punt[i]; __debug_i < s->punt[i+1]; __debug_i++) { \
    __MPS_DEBUG(s, "%d: Approximation: ", __debug_i - s->punt[i]); \
    mpc_out_str(s->logstr, 10, 10, s->mroot[s->clust[__debug_i]]); \
    fprintf(s->logstr, " - Radius: "); \
    rdpe_outln(s->drad[s->clust[__debug_i]]); \
} \
}

/**
 * @brief Print a debug information.
 */
#if __STDC_VERSION__ >= 199901L
#define MPS_DEBUG(s, templ...) __MPS_DEBUG(s,templ) ; \
if (s->DOLOG) { \
    fprintf(s->logstr, "\n"); \
}
#endif

/**
 * @brief Debug the value of a complex multiprecision
 * variable.
 */
#define MPS_DEBUG_MPC(s, digits, c, name...) __MPS_DEBUG_EQ(s, name); \
if (s->DOLOG) { \
    mpc_outln_str(s->logstr, 10, digits, c); \
}

/**
 * @brief Debug the value of a rdpe variable.
 */
#define MPS_DEBUG_RDPE(s, r, name...) __MPS_DEBUG_EQ(s, name); \
if (s->DOLOG) { \
    rdpe_outln_str(s->logstr, r); \
}

/**
 * @brief Debug the value of a cdpe variable.
 */
#define MPS_DEBUG_CDPE(s, c, name...) __MPS_DEBUG_EQ(s, name); \
if (s->DOLOG) { \
    cdpe_outln_str(s->logstr, c); \
}

/**
 * @brief Debug the values of a cplx_t variable
 */
#define MPS_DEBUG_CPLX(s, c, name...)  __MPS_DEBUG_EQ(s, name); \
if (s->DOLOG) { \
    cplx_outln_str(s->logstr, c); \
}

/**
 * @brief Make some space in the debug stream to make clean that
 * another section is starting.
 */
#define MPS_DEBUG_BREAK(s) if (s->DOLOG) { \
    fprintf(s->logstr, "\n"); \
}

/**
 * @brief Low-level debug print that appends an " = " sign to the
 * output (useful for debugging values of variables).
 */
#define __MPS_DEBUG_EQ(s, templ...) \
    __MPS_DEBUG(s, templ); \
    if (s->DOLOG) { \
      fprintf(s->logstr, " = "); \
    }

/**
 * @brief Debug that a function is going to be called.
 */
#ifndef __WINDOWS
#define MPS_DEBUG_CALL(s, function) if (s->DOLOG) { \
  if (isatty(s->logstr->_fileno)) { \
    __MPS_DEBUG(s, "Called \033[31;1m"); \
  } \
else { \
  __MPS_DEBUG(s, "Called "); \
}\
  if (isatty(s->logstr->_fileno)) { \
    fprintf(s->logstr, function); fprintf(s->logstr, "()\033[0m\n"); \
  } \
else \
{ \
	fprintf(s->logstr, function); fprintf(s->logstr, "()\n"); \
}\
}
#else
#define MPS_DEBUG_CALL(s, function) if (s->DOLOG) { \
  if (_isatty(_fileno(s->logstr))) { \
    __MPS_DEBUG(s, "Called \033[31;1m"); \
  } \
else { \
  __MPS_DEBUG(s, "Called "); \
}\
  if (_isatty(_fileno(s->logstr))) { \
    fprintf(s->logstr, function); fprintf(s->logstr, "()\033[0m\n"); \
  } \
else \
{ \
        fprintf(s->logstr, function); fprintf(s->logstr, "()\n"); \
}\
}
#endif

#define MPS_DEBUG_THIS_CALL MPS_DEBUG_CALL(s, __FUNCTION__)


/**
 * @brief Low-level DEBUG() used by other MPS_DEBUG_* statements.
 */
#if __STDC_VERSION__ >= 199901L
#ifndef __WINDOWS
#define __MPS_DEBUG(s, templ...) if (s->DOLOG) {\
if (isatty(s->logstr->_fileno)) { \
    fprintf(s->logstr, "%s:%d \033[32;1m%s()\033[;0m ", \
    __FILE__, __LINE__, __FUNCTION__); \
} \
else { \
    fprintf(s->logstr, "%s:%d %s() ", \
    __FILE__, __LINE__, __FUNCTION__); \
} \
gmp_fprintf(s->logstr, templ); \
}
#else
#define __MPS_DEBUG(s, templ...) if (s->DOLOG) {\
if (_isatty(_fileno(s->logstr))) { \
    fprintf(s->logstr, "%s:%d \033[32;1m%s()\033[;0m ", \
    __FILE__, __LINE__, __FUNCTION__); \
} \
else { \
    fprintf(s->logstr, "%s:%d %s() ", \
    __FILE__, __LINE__, __FUNCTION__); \
} \
gmp_fprintf(s->logstr, templ); \
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
#endif

/*
 * Ansi C version implemented without using variadic macros, only for compatibility
 * since could became imprecise when using complex debug instructions.
 */
#ifndef DISABLE_DEBUG
#if __STDC_VERSION__ < 199901L
#include <mps/interface.h>
#ifdef MPS_DEBUG
#undef MPS_DEBUG
#endif
#ifdef __MPS_DEBUG
#undef __MPS_DEBUG
#endif
#define MPS_DEBUG if (s->DOLOG && mps_is_a_tty(s->logstr))\
    fprintf(s->logstr, "%s:%d \033[32;1m%s()\033[;0m ", __FILE__, __LINE__, __FUNCTION__); \
    if (s->DOLOG && !mps_is_a_tty(s->logstr))\
        fprintf(s->logstr, "%s:%d %s() ", __FILE__, __LINE__, __FUNCTION__); \
    __c_impl__MPS_DEBUG
#define __MPS_DEBUG if (s->DOLOG && mps_is_a_tty(s->logstr))\
    fprintf(s->logstr, "%s:%d \033[32;1m%s()\033[;0m ", __FILE__, __LINE__, __FUNCTION__); \
    if (s->DOLOG && !mps_is_a_tty(s->logstr))\
        fprintf(s->logstr, "%s:%d %s() ", __FILE__, __LINE__, __FUNCTION__); \
    __c_impl____MPS_DEBUG
void __c_impl__MPS_DEBUG(mps_status* s, const char* templ, ...);
void __c_impl____MPS_DEBUG(mps_status* s, const char* templ, ...);
#endif
#endif


#ifdef	__cplusplus
}
#endif

#endif	/* DEBUG_H */

