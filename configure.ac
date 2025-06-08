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
m4_define([WFA2_LIB_VERSION],
          m4_esyscmd([git submodule | awk '$2=="src/WFA2-lib"{printf $1}']))
m4_define([EDLIB_VERSION],
          m4_esyscmd([git submodule | awk '$2=="src/edlib"{printf $1}']))
AC_CHECK_FILE([src/WFA2-lib/Makefile], [],
              [[git clone https://github.com/brendanf/WFA2-lib src/WFA2-lib &&
               cd src/WFA2-lib &&
               git checkout] WFA2_LIB_VERSION [&&
               cd ../..]])
AC_CHECK_FILE([src/edlib/edlib.h], [],
              [[git clone https://github.com/brendanf/edlib src/edlib &&
               cd src/edlib &&
               git checkout] EDLIB_VERSION [&&
               cd ../..]])

### Determine subdirectory source files

# Dynamically generate list of sources from subdirectories via shell
CXXSOURCES="$(cd src/ && ls *.cpp | tr '\n' ' ')"
CSOURCES="$(cd src/ && ls *.c WFA2-lib/**/*.c | tr '\n' ' ')"

# Supply it as a variable
AC_SUBST(CXXSOURCES)
AC_SUBST(CSOURCES)

### Output result

# Write to the file
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
