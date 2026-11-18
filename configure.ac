### Initialize autoconf

# Require the latest version of autoconf
AC_PREREQ([2.69])

# Provide information about package and version numbering
AC_INIT([OptimOTU],m4_bregexp(m4_quote(m4_include([DESCRIPTION])),
                              [Version: \(.*\)], [\1]))

# Supply copyright information
AC_COPYRIGHT(Copyright (C) 2025 Brendan Furneaux)

# Verify the correct src directory exists from where the file is being run.
AC_CONFIG_SRCDIR([src])

### Make sure we have our submodules
# Pins live in tools/submodule-pins (not m4_esyscmd) so they stay in sync
# with the gitlinks. GitHub tarballs ship empty submodule directories, and
# remotes clones the default branch of each submodule before this runs.
AS_IF([sh tools/fetch-submodules.sh],
      [],
      [AC_MSG_ERROR([fetching git submodules failed])])

### Determine subdirectory source files

# Dynamically generate list of sources from subdirectories via shell
CXXSOURCES="$(cd src/ && ls *.cpp | tr '\n' ' ')"
# ksw2: only the affine/dual-affine global+extension implementations
# (exclude cli.c, gg*, splice, and unused helpers)
CSOURCES="$(cd src/ && ls *.c WFA2-lib/**/*.c ksw2/ksw2_extz.c ksw2/ksw2_extd.c ksw2/ksw2_extz2_sse.c ksw2/ksw2_extd2_sse.c | grep -v 'commons\.c$' | tr '\n' ' ')"

# Supply it as a variable
AC_SUBST(CXXSOURCES)
AC_SUBST(CSOURCES)

### Optional CUDA (nvcc). Match the toolkit to the driver, as in
### tools/bench_evo_kernels.sh. /usr/local/cuda may be newer than the
### driver (this machine: driver 12.7, default toolkit 13).
NVCC_BIN=""
CUSOURCES=""
CUDA_CPPFLAGS=""
CUDALIBS=""
CUDA_ARCH_FLAGS=""
DRIVER_MAJOR=`nvidia-smi 2>/dev/null | sed -n 's/.*CUDA Version: \([[0-9]][[0-9]]*\).*/\1/p' | head -n 1`
nvcc_release_major() {
  "[$]1" --version 2>/dev/null | sed -n 's/.*release \([[0-9]][[0-9]]*\).*/\1/p' | head -n 1
}
nvcc_compatible() {
  maj=`nvcc_release_major "[$]1"`
  if test -z "$maj"; then
    return 1
  fi
  if test -n "$DRIVER_MAJOR" && test "$maj" -gt "$DRIVER_MAJOR"; then
    return 1
  fi
  return 0
}
CANDIDATES=""
if test "$DRIVER_MAJOR" = "12"; then
  CANDIDATES="/usr/local/cuda-12.6/bin/nvcc /usr/local/cuda-12/bin/nvcc /usr/local/cuda/bin/nvcc"
elif test "$DRIVER_MAJOR" = "13"; then
  CANDIDATES="/usr/local/cuda-13/bin/nvcc /usr/local/cuda-13.3/bin/nvcc /usr/local/cuda/bin/nvcc"
else
  CANDIDATES="/usr/local/cuda/bin/nvcc"
fi
if test -n "$NVCC" && test -x "$NVCC" && nvcc_compatible "$NVCC"; then
  NVCC_BIN="$NVCC"
else
  if test -n "$NVCC" && test -x "$NVCC"; then
    AC_MSG_NOTICE([CUDA: ignoring NVCC=$NVCC (newer than driver CUDA $DRIVER_MAJOR)])
  fi
  for cand in $CANDIDATES; do
    if test -x "$cand" && nvcc_compatible "$cand"; then
      NVCC_BIN="$cand"
      break
    fi
  done
  if test -z "$NVCC_BIN"; then
    AC_PATH_PROG([NVCC_PATH], [nvcc], [])
    if test -n "$NVCC_PATH" && nvcc_compatible "$NVCC_PATH"; then
      NVCC_BIN="$NVCC_PATH"
    fi
  fi
fi
if test -n "$NVCC_BIN"; then
  cuda_bindir=`dirname "$NVCC_BIN"`
  cuda_home=`dirname "$cuda_bindir"`
  if test -f "$cuda_home/lib64/libcudart.so" || test -f "$cuda_home/lib64/libcudart.dylib"; then
    cuda_libdir="$cuda_home/lib64"
  else
    cuda_libdir="$cuda_home/lib"
  fi
  CUSOURCES="evo_gpu.cu gasal2_gpu.cu"
  CUDA_CPPFLAGS="-DOPTIMOTU_CUDA -I${cuda_home}/include"
  CUDALIBS="-L${cuda_libdir} -lcudart -Wl,-rpath,${cuda_libdir}"
  if test -n "$OPTIMOTU_CUDA_ARCH"; then
    CUDA_ARCH_FLAGS="$OPTIMOTU_CUDA_ARCH"
  else
    CUDA_ARCH_FLAGS="-gencode arch=compute_80,code=sm_80 -gencode arch=compute_86,code=sm_86 -gencode arch=compute_89,code=sm_89 -gencode arch=compute_90,code=sm_90 -gencode arch=compute_90,code=compute_90"
  fi
  AC_MSG_NOTICE([CUDA: using $NVCC_BIN (driver CUDA $DRIVER_MAJOR)])
else
  NVCC_BIN=""
  AC_MSG_NOTICE([CUDA: no driver-compatible nvcc; building CPU-only])
fi
NVCC="$NVCC_BIN"
AC_SUBST([NVCC])
AC_SUBST([CUSOURCES])
AC_SUBST([CUDA_CPPFLAGS])
AC_SUBST([CUDALIBS])
AC_SUBST([CUDA_ARCH_FLAGS])

### Output result

# Write to the file
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
