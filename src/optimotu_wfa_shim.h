#ifndef OPTIMOTU_WFA_SHIM_H
#define OPTIMOTU_WFA_SHIM_H

/*
 * Force-included from PKG_CFLAGS into every C file of the R build.
 * Do not include this from C++: it rewrites stderr and exit, and those
 * macros must not touch libstdc++ or the R error path.
 */
#ifdef __cplusplus
#error "optimotu_wfa_shim.h is C only"
#endif

#include <stdio.h>
#include <stdlib.h>

/* noreturn: a rewritten exit() must not fall through. */
__attribute__((noreturn)) void optimotu_wfa_fail(int code);
extern FILE *optimotu_wfa_null;

#undef stderr
#define stderr optimotu_wfa_null

#undef exit
#define exit(code) optimotu_wfa_fail(code)

#endif
