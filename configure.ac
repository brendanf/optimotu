# Copyright (C) 2018 James J Balamuta
#
# This file is part of SubdirSrc.
#
# SubdirSrc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SubdirSrc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SubdirSrc. If not, see <http://www.gnu.org/licenses/>.

### Initialize autoconf

# Require the latest version of autoconf
AC_PREREQ([2.71])

m4_define([PACKAGE_VERSION],
          m4_esyscmd([awk -e '/^Version:/ {print $2}' DESCRIPTION]))

# Provide information about package and version numbering
AC_INIT([SubdirSrc],m4_defn([PACKAGE_VERSION]))

# Supply copyright information
AC_COPYRIGHT(Copyright (C) 2018 James J Balamuta)

# Verify the correct src directory exists from where the file is being run.
AC_CONFIG_SRCDIR([src])

### Make sure we have our submodules
m4_define([SNEAKY_SNAKE_VERSION],
          m4_esyscmd([git submodule | awk '$2=="src/SneakySnake"{printf $1}']))
m4_define([WFA2_LIB_VERSION],
          m4_esyscmd([git submodule | awk '$2=="src/WFA2-lib"{printf $1}']))
m4_define([EDLIB_VERSION],
          m4_esyscmd([git submodule | awk '$2=="src/edlib"{printf $1}']))
AC_CHECK_FILE([src/SneakySnake/SneakySnake.c], [],
              [[git clone https://github.com/brendanf/SneakySnake src/SneakySnake &&
               cd src/SneakySnake &&
               git checkout] SNEAKY_SNAKE_VERSION [&&
               cd ../..]])
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
CSOURCES="$(cd src/ && ls *.c WFA2-lib/**/*.c SneakySnake/*.c | tr '\n' ' ')"

# Supply it as a variable
AC_SUBST(CXXSOURCES)
AC_SUBST(CSOURCES)

### Output result

# Write to the file
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
