# Original: https://gist.github.com/multiplemonomials/9cc2c4343a4bdd93689520ce64098b84

# CheckLinkerFlag
# ------------------
#
# Checks whether a compiler supports a given linker flag.
#

include(CheckCSourceCompiles)
include(CheckCXXSourceCompiles)
include(CheckFortranSourceCompiles)

include(CMakeCheckCompilerFlagCommonPatterns)

# FLAG: the linker flag to check LANGUAGE: the language to test it on.  C, CXX, or Fortran RESULT: return variable

macro(check_linker_flag FLAG LANGUAGE RESULT)
  if(NOT DEFINED ${RESULT})
    set(OLD_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${FLAG}")
    # Normalize locale during test compilation.
    set(_CheckLinkerFlag_LOCALE_VARS LC_ALL LC_MESSAGES LANG)
    foreach(v ${_CheckLinkerFlag_LOCALE_VARS})
      set(_CheckLinkerFlag_SAVED_${v} "$ENV{${v}}")
      set(ENV{${v}} C)
    endforeach()

    check_compiler_flag_common_patterns(_CheckLinkerFlag_COMMON_PATTERNS)

    if("${LANGUAGE}" STREQUAL C)

      check_c_source_compiles(
        "int main(void) { return 0; }"
        ${RESULT}
        # Some compilers do not fail with a bad flag
        FAIL_REGEX
        "command line option .* is valid for .* but not for C" # GNU
        ${_CheckLinkerFlag_COMMON_PATTERNS})
    elseif("${LANGUAGE}" STREQUAL CXX)
      check_cxx_source_compiles(
        "int main(void) { return 0; }"
        ${RESULT}
        # Some compilers do not fail with a bad flag
        FAIL_REGEX
        "command line option .* is valid for .* but not for C++" # GNU
        ${_CheckLinkerFlag_COMMON_PATTERNS})
    elseif("${LANGUAGE}" STREQUAL Fortran)
      check_fortran_source_compiles(
        "program test
			print *, \'Hello, world\'
			end program test"
        ${RESULT}
        # Some compilers do not fail with a bad flag
        FAIL_REGEX
        "command line option .* is valid for .* but not for Fortran" # GNU
        ${_CheckLinkerFlag_COMMON_PATTERNS})
    else()
      message(FATAL_ERROR "Invalid LANGUAGE argument ${LANGUAGE}")
    endif()

    foreach(v ${_CheckCCompilerFlag_LOCALE_VARS})
      set(ENV{${v}} ${_CheckCCompilerFlag_SAVED_${v}})
      unset(_CheckCCompilerFlag_SAVED_${v})
    endforeach()
    unset(_CheckCCompilerFlag_LOCALE_VARS)
    unset(_CheckCCompilerFlag_COMMON_PATTERNS)

    set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_DEFINITIONS}")
  endif()
endmacro(check_linker_flag)
