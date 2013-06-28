@echo off
setlocal

set args=-MD -MT -debug -lfs -legacy -64 -scope -dll -install -f2c -ifort
set args=%args% -absoft -hdf5 -zlib -szip -mpi
set args=%args% -cgnstools -tcl -tk -nocut -nomesh -winhtml
set copts=
set debug=
set cfgflags=
set dllopts=
set build=
set dolegacy=0
set do64bit=0
set doscope=0
set target=cgns
set make=nmake /nologo
set f2c=
set f77=f77
set hdf5inc=
set hdf5lib=
set zliblib=
set sziplib=
set mpiinc=
set mpilibs=
set mpiexec=
set hdf5dll=#HDF5DLL
set cgnstools=
set tcldir=
set tkdir=
set tclinc=
set tcllib=
set tklib=
set plotopts=
set instdir=c:\cgnslib
set builddir=lib
set winhtml=

set drive=%~d0
set cgnsdir=%~dps0

:next
if "%1" == "" goto doit
if %1 == -help goto usage
for %%a in ( %args% ) do (
  if %1 == %%a goto copts
)
goto badarg

rem ----- compiler options

:copts
if %1 == -MD (
  echo using MD : multi-threaded library with RTL
  set copts=/MD
  shift
  goto next
)
if %1 == -MT (
  echo using MT : multi-threaded library
  set copts=/MT
  shift
  goto next
)
if %1 == -debug (
  echo debug is enabled
  set debug=/Zi
  shift
  goto next
)

rem ----- enable large file support

if %1 == -lfs (
  echo large file support is enabled
  set cfgflags=/DHAVE_LSEEK64
  shift
  goto next
)

rem ----- build legacy code

if %1 == -legacy (
  echo building legacy code version
  set dolegacy=1
  shift
  goto next
)

rem ----- build for 64 bit

if %1 == -64 (
  echo building 64-bit version
  set do64bit=1
  shift
  goto next
)

rem ----- enable enum scoping

if %1 == -scope (
  echo enabling enumeration scoping
  set doscope=1
  shift
  goto next
)

rem ----- build DLL

if %1 == -dll (
  echo building DLL instead of static library
  set target=dll
  set dllopts=/DBUILD_DLL
  shift
  goto next
)

rem ----- ifort Fortran compiler

if %1 == -ifort (
  echo using ifort Fortran compiler with UPPERCASE
  set f77=ifort
  set f2c=/DUPPERCASE
  shift
  goto next
)

rem ----- absoft Fortran compiler

if %1 == -absoft (
  echo using absoft Fortran compiler with LOWERCASE
  set f2c=/DLOWERCASE
  shift
  goto next
)

rem ----- cgnstools options

if %1 == -cgnstools (
  echo building cgsntools
  set cgnstools=cgnstools
  shift
  goto next
)

if %1 == -nocut (
  echo building cgnsplot without cutting plane
  set plotopts=/DNO_CUTTING_PLANE %plotopts%
  shift
  goto next
)

if %1 == -nomesh (
  echo building cgnsplot without structured mesh boundaries
  set plotopts=/DNO_MESH_BOUNDARIES %plotopts%
  shift
  goto next
)

if %1 == -winhtml (
  set winhtml=c:\PROGRA~1\HTMLHE~1
  shift
  goto next
)

rem ----- Fortran to C interface

if not %1 == -f2c goto hdf5
set f2c=/DLOWERCASE_
shift
if "%1" == "" (
  echo using LOWERCASE_ as Fortran to C interface
  goto doit
)
for %%a in ( %args% ) do (
  if %1 == %%a (
    echo using LOWERCASE_ as Fortran to C interface
    goto next
  )
)
if %1 == none (
  echo Fortran interface is disabled
  set f2c=none
  shift
  goto next
)
if %1 == LOWERCASE goto setf2c
if %1 == LOWERCASE_ goto setf2c
if %1 == LOWERCASE__ goto setf2c
if %1 == UPPERCASE goto setf2c
if %1 == UPPERCASE_ goto setf2c
if %1 == UPPERCASE__ goto setf2c
echo ERROR:-f2c argument %1 is invalid
goto usage

:setf2c
echo using %1 as Fortran to C interface
set f2c=/D%1
shift
goto next

rem ----- HDF5 setup

:hdf5
if not %1 == -hdf5 goto zlib
shift

if '%1' == '' goto findhdf5
for %%a in ( %args% ) do if %1 == %%a goto findhdf5

set hdf5dir=%~s1
if not exist %hdf5dir%\nul (
  echo ERROR:HDF5 directory %1 does not exist or is not a directory
  goto done
)
shift
goto gethdf5

:findhdf5
echo checking for HDF5 ...
for /D %%d in ( %drive%\*.* ) do (
  if exist %%d\include\hdf5.h (
    echo %%d
    set hdf5dir=%%d
    goto gethdf5
  )
  if exist %%d\src\hdf5.h (
    echo %%d
    set hdf5dir=%%d
    goto gethdf5
  )
  for /D %%e in ( %%d\*.* ) do (
    if exist %%e\include\hdf5.h (
      echo %%e
      set hdf5dir=%%e
      goto gethdf5
    )
    if exist %%e\src\hdf5.h (
      echo %%e
      set hdf5dir=%%e
      goto gethdf5
    )
    for /D %%f in ( %%e\*.* ) do (
      if exist %%f\include\hdf5.h (
        echo %%f
        set hdf5dir=%%f
        goto gethdf5
      )
      if exist %%f\src\hdf5.h (
        echo %%f
        set hdf5dir=%%f
        goto gethdf5
      )
      for /D %%g in ( %%f\*.* ) do (
        if exist %%g\include\hdf5.h (
          echo %%g
          set hdf5dir=%%g
          goto gethdf5
        )
        if exist %%g\src\hdf5.h (
          echo %%g
          set hdf5dir=%%g
          goto gethdf5
        )
      )
    )
  )
)
echo ERROR:couldn't find hdf5 directory
goto done

:gethdf5
echo checking for hdf5 headers in %hdf5dir% ...
if exist %hdf5dir%\include\hdf5.h (
  echo %hdf5dir%\include
  set hdf5inc=%hdf5dir%\include
  echo checking for hdf5 library in %hdf5dir%\lib ...
  for %%l in ( hdf5dll hdf5 hdf5dlld hdf5d ) do (
    for %%d in ( dll lib ) do (
      if exist %hdf5dir%\%%d\%%l.lib (
        echo %hdf5dir%\%%d\%%l.lib
        set hdf5lib=%hdf5dir%\%%d\%%l.lib
        if %%l == hdf5dll set hdf5dll=HDF5DLL
        if %%l == hdf5dlld set hdf5dll=HDF5DLL
        goto next
      )
    )
  )
  echo ERROR:hdf5 library not found in %hdf5dir%
  goto done
)
if exist %hdf5dir%\src\hdf5.h (
  echo %hdf5dir%\src
  set hdf5inc=%hdf5dir%\src
  echo checking for hdf5 library in %hdf5dir%\proj\ ...
  for %%j in ( Release Debug ) do (
    for %%i in ( hdf5dll hdf5 ) do (
      if exist %hdf5dir%\proj\%%i\%%j\%%i.lib (
        echo %hdf5dir%\proj\%%i\%%j\%%i.lib
        set hdf5lib=%hdf5dir%\proj\%%i\%%j\%%i.lib
        if %%i == hdf5dll set hdf5dll=HDF5DLL
        goto next
      )
      if exist %hdf5dir%\proj\%%i\%%j\%%id.lib (
        echo %hdf5dir%\proj\%%i\%%j\%%id.lib
        set hdf5lib=%hdf5dir%\proj\%%i\%%j\%%id.lib
        if %%i == hdf5dll set hdf5dll=HDF5DLL
        goto next
      )
    )
  )
  echo ERROR:hdf5 library not found in %hdf5dir%\proj
  goto done
)
echo ERROR:hdf5.h not found in %hdf5dir%\include or %hdf5dir%\src
goto done

rem ----- zlib setup

:zlib
if not %1 == -zlib goto szip
echo checking for zlib ...
shift
if '%1' == '' goto findzlib
for %%a in ( %args% ) do if %1 == %%a goto findzlib

if not exist %1 (
  echo ERROR:zlib library %1 doesn't exist
  goto done
)
set zliblib=%~s1
echo %zliblib%
shift
goto next

:findzlib
for %%i in ( zlib zdll ) do (
  if exist %hdf5dir%\lib\%%i.lib (
    echo %hdf5dir%\lib\%%i.lib
    set zliblib=%hdf5dir%\lib\%%i.lib
    goto next
  )
)
for /D %%d in ( %drive%\*.* ) do (
  for %%i in ( zlib zdll ) do (
    if exist %%d\%%i.lib (
      echo %%d\%%i.lib
      set zliblib=%%d\%%i.lib
      goto next
    )
  )
  for /D %%e in ( %%d\*.* ) do (
    for %%i in ( zlib zdll ) do (
      if exist %%e\%%i.lib (
        echo %%e\%%i.lib
        set zliblib=%%e\%%i.lib
        goto next
      )
    )
    for /D %%f in ( %%e\*.* ) do (
      for %%i in ( zlib zdll ) do (
        if exist %%f\%%i.lib (
          echo %%f\%%i.lib
          set zliblib=%%f\%%i.lib
          goto next
        )
      )
      for /D %%g in ( %%f\*.* ) do (
        for %%i in ( zlib zdll ) do (
          if exist %%g\%%i.lib (
            echo %%g\%%i.lib
            set zliblib=%%g\%%i.lib
            goto next
          )
        )
      )
    )
  )
)
echo ERROR:couldn't find zlib or zdll library
goto done

rem ----- szip setup

:szip
if not %1 == -szip goto mpi
echo checking for szip ...
shift
if '%1' == '' goto findszip
for %%a in ( %args% ) do if %1 == %%a goto findszip

if not exist %1 (
  echo ERROR:szip library %1 doesn't exist
  goto done
)
set sziplib=%~s1
echo %sziplib%
shift
goto next

:findszip
for %%i in ( szip szlib szlibdll ) do (
  if exist %hdf5dir%\lib\%%i.lib (
    echo %hdf5dir%\lib\%%i.lib
    set sziplib=%hdf5dir%\lib\%%i.lib
    goto next
)
for /D %%d in ( %drive%\*.* ) do (
  for %%i in ( szip szlib szlibdll ) do (
    if exist %%d\%%i.lib (
      echo %%d\%%i.lib
      set sziplib=%%d\%%i.lib
      goto next
    )
  )
  for /D %%e in ( %%d\*.* ) do (
    for %%i in ( szip szlib szlibdll ) do (
      if exist %%e\%%i.lib (
        echo %%e\%%i.lib
        set sziplib=%%e\%%i.lib
        goto next
      )
    )
    for /D %%f in ( %%e\*.* ) do (
      for %%i in ( szip szlib szlibdll ) do (
        if exist %%f\%%i.lib (
          echo %%f\%%i.lib
          set sziplib=%%f\%%i.lib
          goto next
        )
      )
      for /D %%g in ( %%f\*.* ) do (
        for %%i in ( szip szlib szlibdll ) do (
          if exist %%g\%%i.lib (
            echo %%g\%%i.lib
            set sziplib=%%g\%%i.lib
            goto next
          )
        )
      )
    )
  )
)
echo ERROR:couldn't find szip, szlib or szlibdll library
goto done

rem ----- MPI setup

:mpi
if not %1 == -mpi goto tcl
shift

if '%1' == '' goto findmpi
for %%a in ( %args% ) do if %1 == %%a goto findmpi

set mpidir=%~s1
if not exist %mpidir%\nul (
  echo ERROR:MPI directory %1 does not exist or is not a directory
  goto done
)
shift
goto getmpi

:findmpi
echo checking for MPI ...
for /D %%d in ( c:\PROGRA~1\*.* ) do (
  if exist %%d\include\mpi.h (
    echo %%d
    set mpidir=%%d
    goto getmpi
  )
)
for /D %%d in ( %drive%\*.* ) do (
  if exist %%d\include\mpi.h (
    echo %%d
    set mpidir=%%d
    goto getmpi
  )
  if exist %%d\src\include\mpi.h (
    echo %%d
    set mpidir=%%d
    goto getmpi
  )
  for /D %%e in ( %%d\*.* ) do (
    if exist %%e\include\mpi.h (
      echo %%e
      set mpidir=%%e
      goto getmpi
    )
    if exist %%e\src\include\mpi.h (
      echo %%e
      set mpidir=%%e
      goto getmpi
    )
    for /D %%f in ( %%e\*.* ) do (
      if exist %%f\include\mpi.h (
        echo %%f
        set mpidir=%%f
        goto getmpi
      )
      if exist %%f\src\include\mpi.h (
        echo %%f
        set mpidir=%%f
        goto getmpi
      )
      for /D %%g in ( %%f\*.* ) do (
        if exist %%g\include\mpi.h (
          echo %%g
          set mpidir=%%g
          goto getmpi
        )
        if exist %%g\src\include\mpi.h (
          echo %%g
          set mpidir=%%g
          goto getmpi
        )
      )
    )
  )
)
echo ERROR:couldn't find MPI directory
goto done

:getmpi
echo checking for mpi headers ...
if exist %mpidir%\include\mpi.h (
  echo %mpidir%\include
  set mpiinc=%mpidir%\include
)
if exist %mpidir%\src\include\mpi.h (
  echo %mpidir%\src\include
  set mpiinc=%mpidir%\src\include
)
if "%mpiinc%" == "" (
  echo ERROR:mpi.h not found in %mpidir%\include or %mpidir%\src\include
  goto done
)
rem ----- if using OpenMPI
if exist %mpidir%\share\openmpi\nul set mpiinc=%mpiinc% /DOMPI_IMPORTS
echo checking for mpiexec ...
if exist %mpidir%\bin\mpiexec.exe (
  echo %mpidir%\bin\mpiexec.exe
  set mpiexec=%mpidir%\bin\mpiexec.exe
)
echo checking for MPI library ...
for %%l in ( mpi mpich libmpi mpid mpichd libmpid ) do (
  if exist %mpidir%\lib\%%l.lib (
    echo %mpidir%\lib\%%l.lib
    set mpilibs=%mpidir%\lib\%%l.lib
    goto next
  )
)
echo ERROR:MPI library not found in %mpidir%\lib
goto done

rem ----- tcl directory

:tcl
if not %1 == -tcl goto tk
shift
if '%1' == '' (
  echo ERROR:tcl directory arg to -tcl not given
  goto usage
)
for %%a in ( %args% ) do (
  if %1 == %%a (
    echo ERROR:tcl directory arg to -tcl not given
    goto usage
  )
)
set tcldir=%~s1
if exist %tcldir%\generic\tcl.h goto got_tcldir
if exist %tcldir%\include\tcl.h goto got_tcldir
echo ERROR:can't find tcl.h in %1\include or %1\generic
goto done
:got_tcldir
shift
goto next

rem ----- tk directory

:tk
if not %1 == -tk goto install
shift
if '%1' == '' (
  echo ERROR:tk directory arg to -tk not given
  goto usage
)
for %%a in ( %args% ) do (
  if %1 == %%a (
    echo ERROR:tk directory arg to -tk not given
    goto usage
  )
)
set tkdir=%~s1
if exist %tkdir%\generic\tk.h goto got_tkdir
if exist %tkdir%\include\tk.h goto got_tkdir
echo ERROR:can't find tk.h in %1\include or %1\generic
goto done
:got_tkdir
shift
goto next

rem ----- installation directory

:install
if not %1 == -install goto badarg
shift
if '%1' == '' (
  echo ERROR:installation directory not given
  goto usage
)
for %%a in ( %args% ) do (
  if %1 == %%a (
    echo ERROR:installation directory not given
    goto usage
  )
)
set instdir=%~s1
shift
goto next

rem ----- print usage

:badarg
echo ERROR:unknown argument %1
:usage
echo usage: configure [options]
echo options:
echo   -MT : multi-threaded using libcmt.lib (default)
echo   -MD : multi-threaded using mscvrt.lib and mscvrt.dll
echo   -debug : add debugging to library
echo   -lfs : enable large file support (more than 2Gb)
echo   -legacy : build as legacy code
echo   -64 : build 64-bit version
echo   -scope : enable enumeration scoping by prefixing CG_
echo   -cgnstools : build CGNStools
echo   -dll : build DLL istead of static library
echo   -ifort : use ifort Fortran compiler (implies -f2c UPPERCASE)
echo   -absoft : use the absoft Fortran compiler (implies -f2c LOWERCASE)
echo   -f2c [type] : set Fortran to C interface. "type" is one of LOWERCASE,
echo        LOWERCASE_,LOWERCASE__,UPPERCASE,UPPERCASE_, or UPPERCASE__ If not
echo        given, LOWERCASE_ is used. This option will also use Unix-like
echo        argument passing, instead of Visual Fortran.
echo        If you specify "type" as none, the Fortran interface is disbled.
echo   -hdf5 [hdf5dir] : build HDF5 interface. "hdf5dir" is the HDF5 toplevel
echo        directory. If "hdf5dir" is not given, the current drive is searched.
echo   -zlib [zliblib] : use zlib. "zliblib" is the pathname to the library.
echo        If "zliblib" is not given, the current drive is searched.
echo   -szip [sziplib] : use szip. "sziplib" is the pathname to the library.
echo        If "sziplib" is not given, the current drive is searched.
echo   -mpi [mpidir] : build MPI interface. "mpidir" is the MPI toplevel
echo        directory. If "mpidir" is not given, the current drive is searched.
echo   -tcl tcldir : specify the Tcl source directory
echo   -tk tkdir : specify the Tk source directory
echo   -nocut : build cgnsplot without cutting plane
echo   -nomesh : build cgnsplot without structured mesh boundaries
echo   -install instdir : install to directory "instdir" (default %instdir%)
echo        headers are installed to instdir\include
echo        library is installed to instdir\lib
echo        executables are installed to instdir\bin
goto done

:doit

if "%copts%" == "" set copts=/MT
if %copts% == /MT (
  set mslib=libcmt.lib
) else (
  set mslib=msvcrt.lib
)
if not "%debug%" == "" (
  set copts=%copts%d %debug%
  set lopts=/debug
) else (
  set lopts=/release
)
if %do64bit% == 1 (
  set copts=%copts% /Wp64
  set windir=WIN64
) else (
  set windir=WIN32
)

set libs="ADF"

if "%hdf5inc%" == "" (
  set hdf5def=
  set hdf5lib=
  set zliblib=
  set sziplib=
) else (
  set hdf5def=/I%hdf5inc%
  set libs=%libs HDF5
  set build=%build% /DBUILD_HDF5
)

set mpidef=
if not "%mpiinc%" == "" set mpidef=/I%mpiinc%

if "%f2c%" == "none" (
  set f2cobjs=
  set f2cflags=
) else (
  set f2cobjs=$^(F2COBJS^)
  set f2cflags=%f2c%
)

rem ----- Tcl setup

if not "%tcldir%" == "" goto tclincludes
if not "%tkdir%" == "" (
  if exist %tkdir%\include\tcl.h (
    set tcldir=%tkdir%
    goto tclincludes
  )
)
echo checking for Tcl
if exist c:\Tcl\include\tcl.h (
  set tcldir=c:\Tcl
  goto tclincludes
)
if exist c:\progra~1\Tcl\include\tcl.h (
  set tcldir=c:c:\progra~1\Tcl
  goto tclincludes
)
if exist %drive%\Tcl\include\tcl.h (
  set tcldir=%drive%\Tcl
  goto tclincludes
)
for /D %%d in ( %drive%\*.* ) do (
  if exist %%d\generic\tcl.h (
    set tcldir=%%d
    goto tclincludes
  )
  if exist %%d\include\tcl.h (
    set tcldir=%%d
    goto tclincludes
  )
  for /D %%e in ( %%d\*.* ) do (
    if exist %%e\generic\tcl.h (
      set tcldir=%%e
      goto tclincludes
    )
    if exist %%e\include\tcl.h (
      set tcldir=%%e
      goto tclincludes
    )
    for /D %%f in ( %%e\*.* ) do (
      if exist %%f\generic\tcl.h (
        set tcldir=%%f
        goto tclincludes
      )
      if exist %%f\include\tcl.h (
        set tcldir=%%f
        goto tclincludes
      )
      for /D %%g in ( %%f\*.* ) do (
        if exist %%g\generic\tcl.h (
          set tcldir=%%g
          goto tclincludes
        )
        if exist %%g\include\tcl.h (
          set tcldir=%%g
          goto tclincludes
        )
      )
    )
  )
)
echo couldn't find Tcl directory
echo WARNING:You will not be able to build cgnstools
goto configure

:tclincludes
echo using Tcl directory %tcldir%
if exist %tcldir%\generic\tcl.h (
  set tclinc=/I%tcldir%\generic
) else (
  set tclinc=/I%tcldir%\include
)
for %%i in ( 81 82 83 84 85 86 ) do (
  if exist %tcldir%\lib\tcl%%i.lib (
    set tcllib=%tcldir%\lib\tcl%%i.lib
    goto gettk
  )
  if exist %tcldir%\win\Release\tcl%%i.lib (
    set tcllib=%tcldir%\win\Release\tcl%%i.lib
    goto gettk
  )
  if exist %tcldir%\win\Debug\tcl%%i.lib (
    set tcllib=%tcldir%\win\Debug\tcl%%i.lib
    goto gettk
  )
)
set tcllib=%tcldir%\lib\tcl.lib
echo couldn't find Tcl library - using %tcllib%

:gettk
if "%tkdir%" == "" (
  if exist %tcldir%\include\tk.h set tkdir=%tcldir%
)
if not "%tkdir%" == "" goto tkincludes
echo checking for Tk
for /D %%d in ( %drive%\*.* ) do (
  if exist %%d\generic\tk.h (
    set tkdir=%%d
    goto tkincludes
  )
  if exist %%d\include\tk.h (
    set tkdir=%%d
    goto tkincludes
  )
  for /D %%e in ( %%d\*.* ) do (
    if exist %%e\generic\tk.h (
      set tkdir=%%e
      goto tkincludes
    )
    if exist %%e\include\tk.h (
      set tkdir=%%e
      goto tkincludes
    )
    for /D %%f in ( %%e\*.* ) do (
      if exist %%f\generic\tk.h (
        set tkdir=%%f
        goto tkincludes
      )
      if exist %%f\include\tk.h (
        set tkdir=%%f
        goto tkincludes
      )
      for /D %%g in ( %%f\*.* ) do (
        if exist %%g\generic\tk.h (
          set tkdir=%%g
          goto tkincludes
        )
        if exist %%g\include\tk.h (
          set tkdir=%%g
          goto tkincludes
        )
      )
    )
  )
)
echo couldn't find Tk directory
echo WARNING:You will not be able to build cgnstools
goto configure

:tkincludes
echo using Tk directory %tkdir%
if not %tkdir% == %tcldir% (
  if exist %tkdir%\generic\tk.h (
    set tclinc=%tclinc% /I%tkdir%\generic
  ) else (
    set tclinc=%tclinc% /I%tkdir%\include
  )
)
if exist %tkdir%\include\tkWinInt.h (
  goto gotwinint
)
if exist %tkdir%\win\tkWinInt.h (
  set tclinc=%tclinc% /I%tkdir%\win
  goto gotwinint
)
if exist %tkdir%\include\win\tkWinInt.h (
  set tclinc=%tclinc% /I%tkdir%\include\win
  goto gotwinint
)
echo couldn't find tkWinInt.h in %tkdir%\include or %tkdir%\win or %tkdir%\include\win
echo WARNING:You will not be able to build cgnstools
goto configure

:gotwinint
if exist %tkdir%\xlib\nul set tclinc=%tclinc% /I%tkdir%\xlib

for %%i in ( 81 82 83 84 85 86 ) do (
  if exist %tkdir%\lib\tk%%i.lib (
    set tklib=%tkdir%\lib\tk%%i.lib
    goto configure
  )
  if exist %tkdir%\win\Release\tk%%i.lib (
    set tklib=%tkdir%\win\Release\tk%%i.lib
    goto setopts
  )
  if exist %tkdir%\win\Debug\tk%%i.lib (
    set tklib=%tkdir%\win\Debug\tk%%i.lib
    goto setopts
  )
)
set tklib=%tkdir%\lib\tk.lib
echo couldn't find Tk library - using %tklib%

:configure
echo install library in %instdir%\lib
echo install headers in %instdir%\include
echo install executables in %instdir%\bin

if %dolegacy% == 1 (
  set do64bit=0
  set fbytes=
  if not "%hdf5inc%" == "" (
    set adfinc=install-adf install-adfh
  ) else (
    set adfinc=install-adf
  )
) else (
  set fbytes=*4
  if %do64bit% == 1 set fbytes=*8
  set adfinc=
)

rem ----- create cgnstypes.h

echo creating cgnstypes.h
echo #ifndef CGNSTYPES_H> cgnstypes.h
echo #define CGNSTYPES_H>> cgnstypes.h
echo.>> cgnstypes.h
echo #define CG_BUILD_LEGACY %dolegacy% >> cgnstypes.h
echo #define CG_BUILD_64BIT  %do64bit% >> cgnstypes.h
echo #define CG_BUILD_SCOPE  %doscope% >> cgnstypes.h
echo.>> cgnstypes.h
echo #define CG_MAX_INT32 0x7FFFFFFF>> cgnstypes.h
echo #define CG_LONG_T    __int64>> cgnstypes.h
echo.>> cgnstypes.h
echo #if CG_BUILD_LEGACY>> cgnstypes.h
echo # define CG_SIZEOF_SIZE    32 >> cgnstypes.h
echo # define CG_SIZE_DATATYPE "I4">> cgnstypes.h
echo # define cgerr_t  int>> cgnstypes.h
echo # define cgint_t  int>> cgnstypes.h
echo # define cgsize_t int>> cgnstypes.h
echo # define cgid_t   double>> cgnstypes.h
echo #else>> cgnstypes.h
echo # if CG_BUILD_64BIT>> cgnstypes.h
echo #  define CG_SIZEOF_SIZE    64 >> cgnstypes.h
echo #  define CG_SIZE_DATATYPE "I8">> cgnstypes.h
echo    typedef CG_LONG_T cgsize_t;>> cgnstypes.h
echo # else>> cgnstypes.h
echo #  define CG_SIZEOF_SIZE    32 >> cgnstypes.h
echo #  define CG_SIZE_DATATYPE "I4">> cgnstypes.h
echo    typedef int cgsize_t;>> cgnstypes.h
echo # endif>> cgnstypes.h
echo   typedef int cgerr_t;>> cgnstypes.h
echo   typedef int cgint_t;>> cgnstypes.h
echo   typedef double cgid_t;>> cgnstypes.h
echo #endif>> cgnstypes.h
echo.>> cgnstypes.h
echo typedef CG_LONG_T cglong_t;>> cgnstypes.h
echo typedef unsigned CG_LONG_T cgulong_t;>> cgnstypes.h
echo.>> cgnstypes.h
echo #endif>> cgnstypes.h

rem ----- create cgnstypes_f.h

echo creating cgnstypes_f.h
echo #ifndef CGNSTYPES_F_H> cgnstypes_f.h
echo #define CGNSTYPES_F_H>> cgnstypes_f.h
echo.>> cgnstypes_f.h
echo #define CG_BUILD_64BIT %do64bit% >> cgnstypes_f.h
echo.>> cgnstypes_f.h
echo #if CG_BUILD_64BIT>> cgnstypes_f.h
echo # define cgsize_t integer*8 >> cgnstypes_f.h
echo # define CGSIZE_T integer*8 >> cgnstypes_f.h
echo #else>> cgnstypes_f.h
echo # define cgsize_t integer*4 >> cgnstypes_f.h
echo # define CGSIZE_T integer*4 >> cgnstypes_f.h
echo #endif>> cgnstypes_f.h
echo.>> cgnstypes_f.h
echo #define cglong_t integer*8 >> cgnstypes_f.h
echo #define CGLONG_T integer*8 >> cgnstypes_f.h
echo #define cgid_t   real*8 >> cgnstypes_f.h
echo #define CGID_T   real*8 >> cgnstypes_f.h
echo.>> cgnstypes_f.h
echo #endif>> cgnstypes_f.h

rem ----- create cgnslib_f.h

echo creating cgnslib_f.h
echo ! Fortran version of cgnslib.h> cgnslib_f.h
echo         integer CG_BUILD_64BIT>> cgnslib_f.h
echo         parameter (CG_BUILD_64BIT = %do64bit%)>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      modes for cgns file                                            *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo         integer%fbytes% CG_MODE_READ, CG_MODE_WRITE, CG_MODE_MODIFY>> cgnslib_f.h
echo         parameter (CG_MODE_READ   = 0)>> cgnslib_f.h
echo         parameter (CG_MODE_WRITE  = 1)>> cgnslib_f.h
echo         parameter (CG_MODE_MODIFY = 2)>> cgnslib_f.h
echo !* legacy code support>> cgnslib_f.h
echo         integer%fbytes% MODE_READ, MODE_WRITE, MODE_MODIFY>> cgnslib_f.h
echo         parameter (MODE_READ   = 0)>> cgnslib_f.h
echo         parameter (MODE_WRITE  = 1)>> cgnslib_f.h
echo         parameter (MODE_MODIFY = 2)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      file types                                                     *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo         integer%fbytes% CG_FILE_NONE, CG_FILE_ADF, CG_FILE_HDF5>> cgnslib_f.h
echo         integer%fbytes% CG_FILE_ADF2>> cgnslib_f.h
echo         parameter (CG_FILE_NONE = 0)>> cgnslib_f.h
echo         parameter (CG_FILE_ADF  = 1)>> cgnslib_f.h
echo         parameter (CG_FILE_HDF5 = 2)>> cgnslib_f.h
echo         parameter (CG_FILE_ADF2 = 3)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      some error code                                                *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo         integer%fbytes% CG_OK, CG_ERROR, CG_NODE_NOT_FOUND>> cgnslib_f.h
echo         integer%fbytes% CG_INCORRECT_PATH, CG_NO_INDEX_DIM>> cgnslib_f.h
echo         parameter (CG_OK             = 0)>> cgnslib_f.h
echo         parameter (CG_ERROR          = 1)>> cgnslib_f.h
echo         parameter (CG_NODE_NOT_FOUND = 2)>> cgnslib_f.h
echo         parameter (CG_INCORRECT_PATH = 3)>> cgnslib_f.h
echo         parameter (CG_NO_INDEX_DIM   = 4)>> cgnslib_f.h
echo !* legacy code support>> cgnslib_f.h
echo         integer%fbytes% ALL_OK, ERROR, NODE_NOT_FOUND, INCORRECT_PATH>> cgnslib_f.h
echo         parameter (ALL_OK         = 0)>> cgnslib_f.h
echo         parameter (ERROR          = 1)>> cgnslib_f.h
echo         parameter (NODE_NOT_FOUND = 2)>> cgnslib_f.h
echo         parameter (INCORRECT_PATH = 3)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Dimensional Units                                              *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo         integer%fbytes% CG_UserDefined, CG_Null>> cgnslib_f.h
echo         parameter (CG_Null = 0)>> cgnslib_f.h
echo         parameter (CG_UserDefined = 1)>> cgnslib_f.h
echo !* legacy code support>> cgnslib_f.h
echo         integer%fbytes% Null, UserDefined>> cgnslib_f.h
echo         parameter (Null = 0)>> cgnslib_f.h
echo         parameter (UserDefined = 1)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% Kilogram, Gram, Slug, PoundMass>> cgnslib_f.h
echo         character*32 MassUnitsName(0:5)>> cgnslib_f.h
echo         parameter (Kilogram  = 2)>> cgnslib_f.h
echo         parameter (Gram      = 3)>> cgnslib_f.h
echo         parameter (Slug      = 4)>> cgnslib_f.h
echo         parameter (PoundMass = 5)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% Meter, Centimeter, Millimeter>> cgnslib_f.h
echo         integer%fbytes% Foot, Inch>> cgnslib_f.h
echo         character*32 LengthUnitsName(0:6)>> cgnslib_f.h
echo         parameter (Meter      = 2)>> cgnslib_f.h
echo         parameter (Centimeter = 3)>> cgnslib_f.h
echo         parameter (Millimeter = 4)>> cgnslib_f.h
echo         parameter (Foot       = 5)>> cgnslib_f.h
echo         parameter (Inch       = 6)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% Second>> cgnslib_f.h
echo         character*32 TimeUnitsName(0:2)>> cgnslib_f.h
echo         parameter (Second = 2)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% Kelvin, Celsius, Rankine, Fahrenheit>> cgnslib_f.h
echo         character*32 TemperatureUnitsName(0:5)>> cgnslib_f.h
echo         parameter (Kelvin     = 2)>> cgnslib_f.h
echo         parameter (Celsius    = 3)>> cgnslib_f.h
echo         parameter (Rankine    = 4)>> cgnslib_f.h
echo         parameter (Fahrenheit = 5)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% Degree, Radian>> cgnslib_f.h
echo         character*32 AngleUnitsName(0:3)>> cgnslib_f.h
echo         parameter (Degree = 2)>> cgnslib_f.h
echo         parameter (Radian = 3)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% Ampere, Abampere, Statampere, Edison, auCurrent>> cgnslib_f.h
echo         character*32 ElectricCurrentUnitsName(0:6)>> cgnslib_f.h
echo         parameter (Ampere     = 2)>> cgnslib_f.h
echo         parameter (Abampere   = 3)>> cgnslib_f.h
echo         parameter (Statampere = 4)>> cgnslib_f.h
echo         parameter (Edison     = 5)>> cgnslib_f.h
echo         parameter (auCurrent  = 6)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% Mole, Entities, StandardCubicFoot, StandardCubicMeter>> cgnslib_f.h
echo         character*32 SubstanceAmountUnitsName(0:5)>> cgnslib_f.h
echo         parameter (Mole               = 2)>> cgnslib_f.h
echo         parameter (Entities           = 3)>> cgnslib_f.h
echo         parameter (StandardCubicFoot  = 4)>> cgnslib_f.h
echo         parameter (StandardCubicMeter = 5)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% Candela, Candle, Carcel, Hefner, Violle>> cgnslib_f.h
echo         character*32 LuminousIntensityUnitsName(0:6)>> cgnslib_f.h
echo         parameter (Candela = 2)>> cgnslib_f.h
echo         parameter (Candle  = 3)>> cgnslib_f.h
echo         parameter (Carcel  = 4)>> cgnslib_f.h
echo         parameter (Hefner  = 5)>> cgnslib_f.h
echo         parameter (Violle  = 6)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Data Class                                                     *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo         integer%fbytes% Dimensional, NormalizedByDimensional>> cgnslib_f.h
echo         integer%fbytes% NormalizedByUnknownDimensional>> cgnslib_f.h
echo         integer%fbytes% NondimensionalParameter, DimensionlessConstant>> cgnslib_f.h
echo         character*32 DataClassName(0:6)>> cgnslib_f.h
echo         parameter (Dimensional                    = 2)>> cgnslib_f.h
echo         parameter (NormalizedByDimensional        = 3)>> cgnslib_f.h
echo         parameter (NormalizedByUnknownDimensional = 4)>> cgnslib_f.h
echo         parameter (NondimensionalParameter        = 5)>> cgnslib_f.h
echo         parameter (DimensionlessConstant          = 6)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Grid Location                                                  *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% Vertex, CellCenter, FaceCenter>> cgnslib_f.h
echo         integer%fbytes% IFaceCenter, JFaceCenter, KFaceCenter, EdgeCenter>> cgnslib_f.h
echo         character*32 GridLocationName(0:8)>> cgnslib_f.h
echo         parameter (Vertex      = 2)>> cgnslib_f.h
echo         parameter (CellCenter  = 3)>> cgnslib_f.h
echo         parameter (FaceCenter  = 4)>> cgnslib_f.h
echo         parameter (IFaceCenter = 5)>> cgnslib_f.h
echo         parameter (JFaceCenter = 6)>> cgnslib_f.h
echo         parameter (KFaceCenter = 7)>> cgnslib_f.h
echo         parameter (EdgeCenter  = 8)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Grid Connectivity Types                                        *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% Overset, Abutting, Abutting1to1>> cgnslib_f.h
echo         character*32 GridConnectivityTypeName(0:4)>> cgnslib_f.h
echo         parameter (Overset      = 2)>> cgnslib_f.h
echo         parameter (Abutting     = 3)>> cgnslib_f.h
echo         parameter (Abutting1to1 = 4)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Point Set Types                                                *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% PointList, PointListDonor>> cgnslib_f.h
echo         integer%fbytes% PointRange, PointRangeDonor>> cgnslib_f.h
echo         integer%fbytes% ElementRange, ElementList, CellListDonor>> cgnslib_f.h
echo         character*32 PointSetTypeName(0:8)>> cgnslib_f.h
echo         parameter (PointList       = 2)>> cgnslib_f.h
echo         parameter (PointListDonor  = 3)>> cgnslib_f.h
echo         parameter (PointRange      = 4)>> cgnslib_f.h
echo         parameter (PointRangeDonor = 5)>> cgnslib_f.h
echo         parameter (ElementRange    = 6)>> cgnslib_f.h
echo         parameter (ElementList     = 7)>> cgnslib_f.h
echo         parameter (CellListDonor   = 8)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Governing Equations and Physical Models Types                  *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% FullPotential, Euler>> cgnslib_f.h
echo         integer%fbytes% NSLaminar, NSTurbulent>> cgnslib_f.h
echo         integer%fbytes% NSLaminarIncompressible>> cgnslib_f.h
echo         integer%fbytes% NSTurbulentIncompressible>> cgnslib_f.h
echo         character*32 GoverningEquationsTypeName(0:7)>> cgnslib_f.h
echo         parameter (FullPotential             = 2)>> cgnslib_f.h
echo         parameter (Euler                     = 3)>> cgnslib_f.h
echo         parameter (NSLaminar                 = 4)>> cgnslib_f.h
echo         parameter (NSTurbulent               = 5)>> cgnslib_f.h
echo         parameter (NSLaminarIncompressible   = 6)>> cgnslib_f.h
echo         parameter (NSTurbulentIncompressible = 7)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !** Any model type will accept both ModelTypeNull and ModelTypeUserDefined.>> cgnslib_f.h
echo !** The following models will accept these values as vaild...>> cgnslib_f.h
echo !**>> cgnslib_f.h
echo !** GasModel_t: Ideal, VanderWaals, CaloricallyPerfect, ThermallyPerfect,>> cgnslib_f.h
echo !**    ConstantDensity, RedlichKwong>> cgnslib_f.h
echo !**>> cgnslib_f.h
echo !** ViscosityModel_t: Constant, PowerLaw, SutherlandLaw>> cgnslib_f.h
echo !**>> cgnslib_f.h
echo !** ThermalConductivityModel_t: PowerLaw, SutherlandLaw, ConstantPrandtl>> cgnslib_f.h
echo !**>> cgnslib_f.h
echo !** TurbulenceModel_t: Algebraic_BaldwinLomax, Algebraic_CebeciSmith,>> cgnslib_f.h
echo !**    HalfEquation_JohnsonKing, OneEquation_BaldwinBarth,>> cgnslib_f.h
echo !**    OneEquation_SpalartAllmaras, TwoEquation_JonesLaunder,>> cgnslib_f.h
echo !**    TwoEquation_MenterSST,TwoEquation_Wilcox>> cgnslib_f.h
echo !**>> cgnslib_f.h
echo !** TurbulenceClosure_t: EddyViscosity, ReynoldsStress,>> cgnslib_f.h
echo !**    ReynoldsStressAlgebraic>> cgnslib_f.h
echo !**>> cgnslib_f.h
echo !** ThermalRelaxationModel_t: Frozen, ThermalEquilib, ThermalNonequilib>> cgnslib_f.h
echo !**>> cgnslib_f.h
echo !** ChemicalKineticsModel_t: Frozen, ChemicalEquilibCurveFit,>> cgnslib_f.h
echo !**    ChemicalEquilibMinimization, ChemicalNonequilib>> cgnslib_f.h
echo !**>> cgnslib_f.h
echo !** EMElectricFieldModel_t: Voltage, Interpolated, Constant, Frozen>> cgnslib_f.h
echo !**>> cgnslib_f.h
echo !** EMMagneticFieldModel_t: Interpolated, Constant, Frozen>> cgnslib_f.h
echo !**>> cgnslib_f.h
echo !** EMConductivityModel_t: Constant, Frozen, Equilibrium_LinRessler,>> cgnslib_f.h
echo !**                        Chemistry_LinRessler>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% Ideal, VanderWaals>> cgnslib_f.h
echo         integer%fbytes% Constant>> cgnslib_f.h
echo         integer%fbytes% PowerLaw, SutherlandLaw>> cgnslib_f.h
echo         integer%fbytes% ConstantPrandtl>> cgnslib_f.h
echo         integer%fbytes% EddyViscosity, ReynoldsStress>> cgnslib_f.h
echo         integer%fbytes% ReynoldsStressAlgebraic>> cgnslib_f.h
echo         integer%fbytes% Algebraic_BaldwinLomax, Algebraic_CebeciSmith>> cgnslib_f.h
echo         integer%fbytes% HalfEquation_JohnsonKing, OneEquation_BaldwinBarth>> cgnslib_f.h
echo         integer%fbytes% OneEquation_SpalartAllmaras>> cgnslib_f.h
echo         integer%fbytes% TwoEquation_JonesLaunder>> cgnslib_f.h
echo         integer%fbytes% TwoEquation_MenterSST, TwoEquation_Wilcox>> cgnslib_f.h
echo         integer%fbytes% CaloricallyPerfect, ThermallyPerfect>> cgnslib_f.h
echo         integer%fbytes% ConstantDensity, RedlichKwong>> cgnslib_f.h
echo         integer%fbytes% Frozen, ThermalEquilib, ThermalNonequilib>> cgnslib_f.h
echo         integer%fbytes% ChemicalEquilibCurveFit>> cgnslib_f.h
echo         integer%fbytes% ChemicalEquilibMinimization>> cgnslib_f.h
echo         integer%fbytes% ChemicalNonequilib>> cgnslib_f.h
echo         integer%fbytes% EMElectricField, EMMagneticField, Voltage>> cgnslib_f.h
echo         integer%fbytes% Interpolated>> cgnslib_f.h
echo         integer%fbytes% EMConductivity, Equilibrium_LinRessler>> cgnslib_f.h
echo         integer%fbytes% Chemistry_LinRessler>> cgnslib_f.h
echo         character*32 ModelTypeName(0:35)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         parameter (Ideal                       = 2)>> cgnslib_f.h
echo         parameter (VanderWaals                 = 3)>> cgnslib_f.h
echo         parameter (Constant                    = 4)>> cgnslib_f.h
echo         parameter (PowerLaw                    = 5)>> cgnslib_f.h
echo         parameter (SutherlandLaw               = 6)>> cgnslib_f.h
echo         parameter (ConstantPrandtl             = 7)>> cgnslib_f.h
echo         parameter (EddyViscosity               = 8)>> cgnslib_f.h
echo         parameter (ReynoldsStress              = 9)>> cgnslib_f.h
echo         parameter (ReynoldsStressAlgebraic     = 10)>> cgnslib_f.h
echo         parameter (Algebraic_BaldwinLomax      = 11)>> cgnslib_f.h
echo         parameter (Algebraic_CebeciSmith       = 12)>> cgnslib_f.h
echo         parameter (HalfEquation_JohnsonKing    = 13)>> cgnslib_f.h
echo         parameter (OneEquation_BaldwinBarth    = 14)>> cgnslib_f.h
echo         parameter (OneEquation_SpalartAllmaras = 15)>> cgnslib_f.h
echo         parameter (TwoEquation_JonesLaunder    = 16)>> cgnslib_f.h
echo         parameter (TwoEquation_MenterSST       = 17)>> cgnslib_f.h
echo         parameter (TwoEquation_Wilcox          = 18)>> cgnslib_f.h
echo         parameter (CaloricallyPerfect          = 19)>> cgnslib_f.h
echo         parameter (ThermallyPerfect            = 20)>> cgnslib_f.h
echo         parameter (ConstantDensity             = 21)>> cgnslib_f.h
echo         parameter (RedlichKwong                = 22)>> cgnslib_f.h
echo         parameter (Frozen                      = 23)>> cgnslib_f.h
echo         parameter (ThermalEquilib              = 24)>> cgnslib_f.h
echo         parameter (ThermalNonequilib           = 25)>> cgnslib_f.h
echo         parameter (ChemicalEquilibCurveFit     = 26)>> cgnslib_f.h
echo         parameter (ChemicalEquilibMinimization = 27)>> cgnslib_f.h
echo         parameter (ChemicalNonequilib          = 28)>> cgnslib_f.h
echo         parameter (EMElectricField             = 29)>> cgnslib_f.h
echo         parameter (EMMagneticField             = 30)>> cgnslib_f.h
echo         parameter (EMConductivity              = 31)>> cgnslib_f.h
echo         parameter (Voltage                     = 32)>> cgnslib_f.h
echo         parameter (Interpolated                = 33)>> cgnslib_f.h
echo         parameter (Equilibrium_LinRessler      = 34)>> cgnslib_f.h
echo         parameter (Chemistry_LinRessler        = 35)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Boundary Condition Types                                       *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% BCAxisymmetricWedge, BCDegenerateLine>> cgnslib_f.h
echo         integer%fbytes% BCDegeneratePoint>> cgnslib_f.h
echo         integer%fbytes% BCDirichlet, BCExtrapolate, BCFarfield, BCGeneral>> cgnslib_f.h
echo         integer%fbytes% BCInflow, BCInflowSubsonic,  BCInflowSupersonic>> cgnslib_f.h
echo         integer%fbytes% BCNeumann>> cgnslib_f.h
echo         integer%fbytes% BCOutflow, BCOutflowSubsonic, BCOutflowSupersonic>> cgnslib_f.h
echo         integer%fbytes% BCSymmetryPlane, BCSymmetryPolar>> cgnslib_f.h
echo         integer%fbytes% BCTunnelInflow, BCTunnelOutflow>> cgnslib_f.h
echo         integer%fbytes% BCWall, BCWallInviscid, BCWallViscous>> cgnslib_f.h
echo         integer%fbytes% BCWallViscousHeatFlux, BCWallViscousIsothermal>> cgnslib_f.h
echo         integer%fbytes% FamilySpecified>> cgnslib_f.h
echo         character*32 BCTypeName(0:25)>> cgnslib_f.h
echo         parameter (BCAxisymmetricWedge     = 2)>> cgnslib_f.h
echo         parameter (BCDegenerateLine        = 3)>> cgnslib_f.h
echo         parameter (BCDegeneratePoint       = 4)>> cgnslib_f.h
echo         parameter (BCDirichlet             = 5)>> cgnslib_f.h
echo         parameter (BCExtrapolate           = 6)>> cgnslib_f.h
echo         parameter (BCFarfield              = 7)>> cgnslib_f.h
echo         parameter (BCGeneral               = 8)>> cgnslib_f.h
echo         parameter (BCInflow                = 9)>> cgnslib_f.h
echo         parameter (BCInflowSubsonic        = 10)>> cgnslib_f.h
echo         parameter (BCInflowSupersonic      = 11)>> cgnslib_f.h
echo         parameter (BCNeumann               = 12)>> cgnslib_f.h
echo         parameter (BCOutflow               = 13)>> cgnslib_f.h
echo         parameter (BCOutflowSubsonic       = 14)>> cgnslib_f.h
echo         parameter (BCOutflowSupersonic     = 15)>> cgnslib_f.h
echo         parameter (BCSymmetryPlane         = 16)>> cgnslib_f.h
echo         parameter (BCSymmetryPolar         = 17)>> cgnslib_f.h
echo         parameter (BCTunnelInflow          = 18)>> cgnslib_f.h
echo         parameter (BCTunnelOutflow         = 19)>> cgnslib_f.h
echo         parameter (BCWall                  = 20)>> cgnslib_f.h
echo         parameter (BCWallInviscid          = 21)>> cgnslib_f.h
echo         parameter (BCWallViscous           = 22)>> cgnslib_f.h
echo         parameter (BCWallViscousHeatFlux   = 23)>> cgnslib_f.h
echo         parameter (BCWallViscousIsothermal = 24)>> cgnslib_f.h
echo         parameter (FamilySpecified         = 25)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Data types                                                     *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% Integer, RealSingle, RealDouble, Character>> cgnslib_f.h
echo         integer%fbytes% LongInteger>> cgnslib_f.h
echo         character*32 DataTypeName(0:6)>> cgnslib_f.h
echo         parameter (Integer     = 2)>> cgnslib_f.h
echo         parameter (RealSingle  = 3)>> cgnslib_f.h
echo         parameter (RealDouble  = 4)>> cgnslib_f.h
echo         parameter (Character   = 5)>> cgnslib_f.h
echo         parameter (LongInteger = 6)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      BCData_t types                                                 *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% Dirichlet, Neumann>> cgnslib_f.h
echo         character*32 BCDataTypeName(0:3)>> cgnslib_f.h
echo         parameter (Dirichlet = 2)>> cgnslib_f.h
echo         parameter (Neumann   = 3)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Element types                                                  *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% NODE, BAR_2, BAR_3, TRI_3, TRI_6>> cgnslib_f.h
echo         integer%fbytes% QUAD_4, QUAD_8, QUAD_9>> cgnslib_f.h
echo         integer%fbytes% TETRA_4, TETRA_10, PYRA_5, PYRA_14>> cgnslib_f.h
echo         integer%fbytes% PENTA_6, PENTA_15, PENTA_18>> cgnslib_f.h
echo         integer%fbytes% HEXA_8, HEXA_20, HEXA_27>> cgnslib_f.h
echo         integer%fbytes% MIXED, PYRA_13, NGON_n, NFACE_n>> cgnslib_f.h
echo         character*32 ElementTypeName(0:23)>> cgnslib_f.h
echo         parameter (NODE     =  2)>> cgnslib_f.h
echo         parameter (BAR_2    =  3)>> cgnslib_f.h
echo         parameter (BAR_3    =  4)>> cgnslib_f.h
echo         parameter (TRI_3    =  5)>> cgnslib_f.h
echo         parameter (TRI_6    =  6)>> cgnslib_f.h
echo         parameter (QUAD_4   =  7)>> cgnslib_f.h
echo         parameter (QUAD_8   =  8)>> cgnslib_f.h
echo         parameter (QUAD_9   =  9)>> cgnslib_f.h
echo         parameter (TETRA_4  = 10)>> cgnslib_f.h
echo         parameter (TETRA_10 = 11)>> cgnslib_f.h
echo         parameter (PYRA_5   = 12)>> cgnslib_f.h
echo         parameter (PYRA_14  = 13)>> cgnslib_f.h
echo         parameter (PENTA_6  = 14)>> cgnslib_f.h
echo         parameter (PENTA_15 = 15)>> cgnslib_f.h
echo         parameter (PENTA_18 = 16)>> cgnslib_f.h
echo         parameter (HEXA_8   = 17)>> cgnslib_f.h
echo         parameter (HEXA_20  = 18)>> cgnslib_f.h
echo         parameter (HEXA_27  = 19)>> cgnslib_f.h
echo         parameter (MIXED    = 20)>> cgnslib_f.h
echo         parameter (PYRA_13  = 21)>> cgnslib_f.h
echo         parameter (NGON_n   = 22)>> cgnslib_f.h
echo         parameter (NFACE_n  = 23)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Zone types                                                     *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% Structured, Unstructured>> cgnslib_f.h
echo         character*32 ZoneTypeName(0:3)>> cgnslib_f.h
echo         parameter (Structured   =  2)>> cgnslib_f.h
echo         parameter (Unstructured =  3)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Rigid Grid Motion types                                        *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% ConstantRate, VariableRate>> cgnslib_f.h
echo         character*32 RigidGridMotionTypeName(0:3)>> cgnslib_f.h
echo         parameter (ConstantRate = 2)>> cgnslib_f.h
echo         parameter (VariableRate = 3)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Arbitrary Grid Motion types                                    *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% NonDeformingGrid, DeformingGrid>> cgnslib_f.h
echo         character*32 ArbitraryGridMotionTypeName(0:3)>> cgnslib_f.h
echo         parameter (NonDeformingGrid = 2)>> cgnslib_f.h
echo         parameter (DeformingGrid = 3)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Simulation type                                                *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% TimeAccurate, NonTimeAccurate>> cgnslib_f.h
echo         character*32 SimulationTypeName(0:3)>> cgnslib_f.h
echo         parameter (TimeAccurate = 2)>> cgnslib_f.h
echo         parameter (NonTimeAccurate = 3)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      BC Property types                                              *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% Generic>> cgnslib_f.h
echo         character*32 WallFunctionTypeName(0:2)>> cgnslib_f.h
echo         parameter (Generic = 2)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% BleedArea, CaptureArea>> cgnslib_f.h
echo         character*32 AreaTypeName(0:3)>> cgnslib_f.h
echo         parameter (BleedArea = 2)>> cgnslib_f.h
echo         parameter (CaptureArea = 3)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Grid Connectivity Property types                               *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         integer%fbytes% AverageAll, AverageCircumferential, AverageRadial>> cgnslib_f.h
echo         integer%fbytes% AverageI, AverageJ, AverageK>> cgnslib_f.h
echo         character*32 AverageInterfaceTypeName(0:7)>> cgnslib_f.h
echo         parameter (AverageAll = 2)>> cgnslib_f.h
echo         parameter (AverageCircumferential = 3)>> cgnslib_f.h
echo         parameter (AverageRadial = 4)>> cgnslib_f.h
echo         parameter (AverageI = 5)>> cgnslib_f.h
echo         parameter (AverageJ = 6)>> cgnslib_f.h
echo         parameter (AverageK = 7)>> cgnslib_f.h
echo.>> cgnslib_f.h
echo ! For portability to Linux Absoft, all data statements were moved after the>> cgnslib_f.h
echo ! variables and parametres declarations>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Dimensional Units                                              *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo         data MassUnitsName /'Null','UserDefined','Kilogram','Gram',     ^&>> cgnslib_f.h
echo      ^&                      'Slug','PoundMass'/>> cgnslib_f.h
echo         data LengthUnitsName / 'Null', 'UserDefined',                   ^&>> cgnslib_f.h
echo      ^&         'Meter','Centimeter','Millimeter','Foot','Inch'/>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data TimeUnitsName /'Null','UserDefined','Second'/>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data TemperatureUnitsName /'Null','UserDefined',                ^&>> cgnslib_f.h
echo      ^&         'Kelvin','Celsius','Rankine','Fahrenheit'/>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data AngleUnitsName /'Null','UserDefined','Degree','Radian'/>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data ElectricCurrentUnitsName /'Null', 'UserDefined', 'Ampere', ^&>> cgnslib_f.h
echo      ^&         'Abampere', 'Statampere', 'Edison', 'a.u.'/>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data SubstanceAmountUnitsName /'Null', 'UserDefined', 'Mole',   ^&>> cgnslib_f.h
echo      ^&         'Entities', 'StandardCubicFoot', 'StandardCubicMeter'/>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data LuminousIntensityUnitsName /'Null', 'UserDefined',         ^&>> cgnslib_f.h
echo      ^&         'Candela', 'Candle', 'Carcel', 'Hefner', 'Violle'/>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Data Class                                                     *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo         data DataClassName / 'Null','UserDefined',                      ^&>> cgnslib_f.h
echo      ^&          'Dimensional','NormalizedByDimensional',                ^&>> cgnslib_f.h
echo      ^&          'NormalizedByUnknownDimensional',                       ^&>> cgnslib_f.h
echo      ^&          'NondimensionalParameter','DimensionlessConstant'/>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Grid Location                                                  *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data GridLocationName / 'Null','UserDefined',                   ^&>> cgnslib_f.h
echo      ^&          'Vertex','CellCenter','FaceCenter','IFaceCenter',       ^&>> cgnslib_f.h
echo      ^&          'JFaceCenter','KFaceCenter','EdgeCenter' />> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Grid Connectivity Types                                        *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data GridConnectivityTypeName / 'Null','UserDefined',           ^&>> cgnslib_f.h
echo      ^&          'Overset','Abutting','Abutting1to1'/>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Point Set Types                                                *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data PointSetTypeName / 'Null','UserDefined',                   ^&>> cgnslib_f.h
echo      ^&          'PointList','PointListDonor',                           ^&>> cgnslib_f.h
echo      ^&          'PointRange','PointRangeDonor',                         ^&>> cgnslib_f.h
echo      ^&          'ElementRange','ElementList','CellListDonor'/>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Governing Equations and Physical Models Types                  *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data GoverningEquationsTypeName / 'Null','UserDefined',         ^&>> cgnslib_f.h
echo      ^&          'FullPotential','Euler', 'NSLaminar', 'NSTurbulent',    ^&>> cgnslib_f.h
echo      ^&          'NSLaminarIncompressible', 'NSTurbulentIncompressible'/>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data ModelTypeName / 'Null','UserDefined',                      ^&>> cgnslib_f.h
echo      ^&        'Ideal','VanderWaals', 'Constant','PowerLaw',             ^&>> cgnslib_f.h
echo      ^&        'SutherlandLaw','ConstantPrandtl','EddyViscosity',        ^&>> cgnslib_f.h
echo      ^&        'ReynoldsStress','ReynoldsStressAlgebraic',               ^&>> cgnslib_f.h
echo      ^&        'Algebraic_BaldwinLomax','Algebraic_CebeciSmith',         ^&>> cgnslib_f.h
echo      ^&        'HalfEquation_JohnsonKing','OneEquation_BaldwinBarth',    ^&>> cgnslib_f.h
echo      ^&        'OneEquation_SpalartAllmaras','TwoEquation_JonesLaunder', ^&>> cgnslib_f.h
echo      ^&        'TwoEquation_MenterSST','TwoEquation_Wilcox',             ^&>> cgnslib_f.h
echo      ^&        'CaloricallyPerfect', 'ThermallyPerfect',                 ^&>> cgnslib_f.h
echo      ^&        'ConstantDensity', 'RedlichKwong', 'Frozen',              ^&>> cgnslib_f.h
echo      ^&        'ThermalEquilib', 'ThermalNonequilib',                    ^&>> cgnslib_f.h
echo      ^&        'ChemicalEquilibCurveFit', 'ChemicalEquilibMinimization', ^&>> cgnslib_f.h
echo      ^&        'ChemicalNonequilib', 'EMElectricField',                  ^&>> cgnslib_f.h
echo      ^&        'EMMagneticField', 'EMConductivity', 'Voltage',           ^&>> cgnslib_f.h
echo      ^&        'Interpolated', 'Equilibrium_LinRessler',                 ^&>> cgnslib_f.h
echo      ^&        'Chemistry_LinRessler'/>> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Boundary Condition Types                                       *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data BCTypeName / 'Null','UserDefined',                         ^&>> cgnslib_f.h
echo      ^&          'BCAxisymmetricWedge','BCDegenerateLine',               ^&>> cgnslib_f.h
echo      ^&          'BCDegeneratePoint','BCDirichlet','BCExtrapolate',      ^&>> cgnslib_f.h
echo      ^&          'BCFarfield','BCGeneral','BCInflow','BCInflowSubsonic', ^&>> cgnslib_f.h
echo      ^&          'BCInflowSupersonic','BCNeumann','BCOutflow',           ^&>> cgnslib_f.h
echo      ^&          'BCOutflowSubsonic','BCOutflowSupersonic',              ^&>> cgnslib_f.h
echo      ^&          'BCSymmetryPlane','BCSymmetryPolar','BCTunnelInflow',   ^&>> cgnslib_f.h
echo      ^&          'BCTunnelOutflow','BCWall','BCWallInviscid',            ^&>> cgnslib_f.h
echo      ^&          'BCWallViscous','BCWallViscousHeatFlux',                ^&>> cgnslib_f.h
echo      ^&          'BCWallViscousIsothermal','FamilySpecified' />> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Data types                                                     *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data DataTypeName / 'Null','UserDefined',                       ^&>> cgnslib_f.h
echo      ^&          'Integer','RealSingle','RealDouble','Character',        ^&>> cgnslib_f.h
echo      ^&          'LongInteger' />> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      BCData_t types                                                 *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data BCDataTypeName / 'Null','UserDefined',                     ^&>> cgnslib_f.h
echo      ^&          'Dirichlet', 'Neumann' />> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Element types                                                  *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data ElementTypeName / 'Null','UserDefined',                    ^&>> cgnslib_f.h
echo      ^&      'NODE', 'BAR_2', 'BAR_3', 'TRI_3', 'TRI_6',                 ^&>> cgnslib_f.h
echo      ^&      'QUAD_4', 'QUAD_8', 'QUAD_9', 'TETRA_4', 'TETRA_10',        ^&>> cgnslib_f.h
echo      ^&      'PYRA_5', 'PYRA_14', 'PENTA_6', 'PENTA_15',                 ^&>> cgnslib_f.h
echo      ^&      'PENTA_18', 'HEXA_8', 'HEXA_20', 'HEXA_27', 'MIXED',        ^&>> cgnslib_f.h
echo      ^&      'PYRA_13', 'NGON_n', 'NFACE_n' />> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Zone types                                                     *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data ZoneTypeName / 'Null','UserDefined',                       ^&>> cgnslib_f.h
echo      ^&      'Structured', 'Unstructured' />> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Rigid Grid Motion types                                        *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data RigidGridMotionTypeName / 'Null','UserDefined',            ^&>> cgnslib_f.h
echo      ^&       'ConstantRate', 'VariableRate' />> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Arbitrary Grid Motion types                                    *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data ArbitraryGridMotionTypeName / 'Null','UserDefined',        ^&>> cgnslib_f.h
echo      ^&       'NonDeformingGrid', 'DeformingGrid' />> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Simulation type                                                *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data SimulationTypeName / 'Null','UserDefined',                 ^&>> cgnslib_f.h
echo      ^&       'TimeAccurate', 'NonTimeAccurate' />> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      BC Property types                                              *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data WallFunctionTypeName / 'Null','UserDefined',               ^&>> cgnslib_f.h
echo      ^&       'Generic' />> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data AreaTypeName / 'Null','UserDefined',                       ^&>> cgnslib_f.h
echo      ^&       'BleedArea', 'CaptureArea' />> cgnslib_f.h
echo.>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo !*      Grid Connectivity Property types                               *>> cgnslib_f.h
echo !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *>> cgnslib_f.h
echo.>> cgnslib_f.h
echo         data AverageInterfaceTypeName / 'Null','UserDefined',           ^&>> cgnslib_f.h
echo      ^&       'AverageAll', 'AverageCircumferential', 'AverageRadial',   ^&>> cgnslib_f.h
echo      ^&       'AverageI', 'AverageJ', 'AverageK' />> cgnslib_f.h

rem ----- create make.defs

echo creating make.defs
echo # makefile include for %windir%> make.defs
echo.>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo # CC      - C compiler>> make.defs
echo # CFLAGS  - compiler flags>> make.defs
echo # COOUT   - flag to name object output file (typically /Fo)>> make.defs
echo # CEOUT   - flag to name the output executable (typically /Fe)>> make.defs
echo # LDFLAGS - any additional linker options>> make.defs
echo # CLIBS   - any additional libraries needed to link a CGNS application>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo.>> make.defs
echo CC      = cl>> make.defs
echo CFLAGS  = /nologo %copts% /D_CRT_SECURE_NO_WARNINGS>> make.defs
echo COOUT   = /Fo>> make.defs
echo CEOUT   = /Fe>> make.defs
echo LDFLAGS = /nologo %lopts%>> make.defs
echo CLIBS   =>> make.defs
echo.>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo # SPACE  - used to force a space in the compiler executable output flag>> make.defs
echo # O       - object file extension>> make.defs
echo # A       - library file extension>> make.defs
echo # EXE     - executable extension>> make.defs
echo # LIBCGNS - CGNS library name>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo.>> make.defs
echo SPACE   =>> make.defs
echo O       = obj>> make.defs
echo A       = lib>> make.defs
echo EXE     = .exe>> make.defs
echo.>> make.defs
echo LIBCGNS = %builddir%\cgns.lib>> make.defs
echo.>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo # F2CFLAGS defines the type of Fortran to C interface.>> make.defs
echo #>> make.defs
echo # CFGFLAGS defines any additional compiler options needed to build>> make.defs
echo # the CGNS library. This is typically set by the configure script.>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo.>> make.defs
echo F2CFLAGS = %f2cflags%>> make.defs
echo CFGFLAGS = %cfgflags%>> make.defs
echo.>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo # These are not used to build the CGNS library>> make.defs
echo # Fortran compiler (F77) and options (FFLAGS).>> make.defs
echo # FEOUT is the flag to name the output executable (typically /exe:).>> make.defs
echo # FLIBS lists any additional libraries needed to link a CGNS application>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo.>> make.defs
echo F77    = %f77%>> make.defs
if %f77% == ifort (
echo FFLAGS = /nologo /extfpp:F /fpp %copts%>> make.defs
) else (
echo FFLAGS = /nologo /extfpp:F /fpp:"/DWINNT" %copts%>> make.defs
)
echo FEOUT  = /exe:>> make.defs
echo FLIBS  =>> make.defs
echo.>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo # strip command for executables - set to true if not used>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo.>> make.defs
echo STRIP  = :>> make.defs
echo.>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo # library archiver and ranlib>> make.defs
echo # AROUT may be used to set a library output flag as:>> make.defs
echo #    $(AR) $(AROUT)library_name objects>> make.defs
echo # Set RANLIB to true if not used>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo.>> make.defs
echo AR     = link>> make.defs
echo AROUT  = /out:>> make.defs
echo RANLIB = :>> make.defs
echo.>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo # these should only be set if building HDF5 interface>> make.defs
echo # HDF5INC - path to HDF5 header files>> make.defs
echo # HDF5LIB - HDF5 library>> make.defs
echo # SZIPLIB - szip library (if needed)>> make.defs
echo # ZLIBLIB - zlib library (if needed)>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo.>> make.defs
echo HDF5INC = %hdf5def%>> make.defs
echo HDF5LIB = %hdf5lib%>> make.defs
echo SZIPLIB = %sziplib%>> make.defs
echo ZLIBLIB = %zliblib%>> make.defs
echo.>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo # these should only be set if building with HDF5 and MPI>> make.defs
echo # MPIINC  - path to MPI header files>> make.defs
echo # MPILIBS - MPI libraries>> make.defs
echo # MPIEXEC - MPI executor>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo.>> make.defs
echo MPIINC  = %mpidef%>> make.defs
echo MPILIBS = %mpilibs%>> make.defs
echo MPIEXEC = %mpiexec%>> make.defs
echo.>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo # BUILDLIBS contains the list of additional libraries>> make.defs
echo #           with which a CGNS application needs to link>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo.>> make.defs
echo BUILDLIBS = $(HDF5LIB) $(SZIPLIB) $(ZLIBLIB) $(MPILIBS)>> make.defs
echo.>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo # commands for removing files and creating/deleting directory>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo.>> make.defs
echo RM     = del /q>> make.defs
echo RMDIR  = rd /s/q>> make.defs
echo MKDIR  = md>> make.defs
echo.>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo # installation library name and directories>> make.defs
echo #>> make.defs
echo # INSTALL      - install command>> make.defs
echo # INSTALL_PROG - install executable>> make.defs
echo # INSTALL_DATA - install data>> make.defs
echo # LIBDIR       - installation directory for CGNS library>> make.defs
echo # INCLUDEDIR   - installation directory for CGNS headers>> make.defs
echo # BINDIR       - installation directory for CGNS executables>> make.defs
echo #------------------------------------------------------------------------>> make.defs
echo.>> make.defs
echo INSTALL      = copy /b>> make.defs
echo INSTALL_PROG = $(INSTALL)>> make.defs
echo INSTALL_DATA = $(INSTALL)>> make.defs
echo INSTALLDIR   = %instdir%>> make.defs
echo LIBDIR       = $(INSTALLDIR)\lib>> make.defs
echo INCLUDEDIR   = $(INSTALLDIR)\include>> make.defs
echo BINDIR       = $(INSTALLDIR)\bin>> make.defs

rem ----- create cgnsBuild.defs

set dodebug=0
if not "%debug%" == "" set dodebug=1
set dofortran=0
if not "%f2c%" == "none" set dofortran=1

echo creating cgnsBuild.defs
echo # makefile include for building CGNS code under %windir%> cgnsBuild.defs
echo # this file contains the options and libraries used for>> cgnsBuild.defs
echo # building and linking CGNS code, and is intended to be>> cgnsBuild.defs
echo # included in a user's Makefile from the installation>> cgnsBuild.defs
echo.>> cgnsBuild.defs
echo #----------------------------------------------------------------------->> cgnsBuild.defs
echo # CGNS library build options. A 1 indicates that the library>> cgnsBuild.defs
echo # was built with that option, a 0 indicates without>> cgnsBuild.defs
echo # CGNS_DEBUG   = debug>> cgnsBuild.defs
echo # CGNS_LEGACY  = legacy code (prior to 3.0)>> cgnsBuild.defs
echo # CGNS_SCOPING = scoping of enums>> cgnsBuild.defs
echo # CGNS_64BIT   = 64 bit support>> cgnsBuild.defs
echo # CGNS_FORTRAN = Fortran interface>> cgnsBuild.defs
echo #----------------------------------------------------------------------->> cgnsBuild.defs
echo.>> cgnsBuild.defs
echo CGNS_DEBUG   = %dodebug% >> cgnsBuild.defs
echo CGNS_LEGACY  = %dolegacy% >> cgnsBuild.defs
echo CGNS_SCOPING = %doscope% >> cgnsBuild.defs
echo CGNS_64BIT   = %do64bit% >> cgnsBuild.defs
echo CGNS_FORTRAN = %dofortran% >> cgnsBuild.defs
echo.>> cgnsBuild.defs
echo #------------------------------------------------------------------------>> cgnsBuild.defs
echo # CGNS_LIBDIR     - installation directory for CGNS library>> cgnsBuild.defs
echo # CGNS_INCLUDEDIR - installation directory for CGNS headers>> cgnsBuild.defs
echo #------------------------------------------------------------------------>> cgnsBuild.defs
echo.>> cgnsBuild.defs
echo CGNS_LIBDIR     = %instdir%\lib>> cgnsBuild.defs
echo CGNS_INCLUDEDIR = %instdir%\include>> cgnsBuild.defs
echo.>> cgnsBuild.defs
echo #------------------------------------------------------------------------>> cgnsBuild.defs
echo # CGNS_CC      - C compiler used to build library>> cgnsBuild.defs
echo # CGNS_CFLAGS  - compiler flags used to build library>> cgnsBuild.defs
echo # CGNS_LDFLAGS - any additional linker options>> cgnsBuild.defs
echo #------------------------------------------------------------------------>> cgnsBuild.defs
echo.>> cgnsBuild.defs
echo CGNS_CC      = cl>> cgnsBuild.defs
echo CGNS_CFLAGS  = /nologo %copts% /D_CRT_SECURE_NO_WARNINGS>> cgnsBuild.defs
echo CGNS_LDFLAGS = /nologo %lopts%>> cgnsBuild.defs
echo.>> cgnsBuild.defs
echo #------------------------------------------------------------------------>> cgnsBuild.defs
echo # CGNS_LIB - CGNS library name>> cgnsBuild.defs
echo #------------------------------------------------------------------------>> cgnsBuild.defs
echo.>> cgnsBuild.defs
echo CGNS_LIB = $(CGNS_LIBDIR)\cgns.lib>> cgnsBuild.defs
echo.>> cgnsBuild.defs
echo #------------------------------------------------------------------------>> cgnsBuild.defs
echo # CGNS_HDF5INC - path to HDF5 header files>> cgnsBuild.defs
echo # CGNS_HDF5LIB - HDF5 library>> cgnsBuild.defs
echo # CGNS_SZIPLIB - szip library (if needed)>> cgnsBuild.defs
echo # CGNS_ZLIBLIB - zlib library (if needed)>> cgnsBuild.defs
echo #------------------------------------------------------------------------>> cgnsBuild.defs
echo.>> cgnsBuild.defs
echo CGNS_HDF5INC = %hdf5inc%>> cgnsBuild.defs
echo CGNS_HDF5LIB = %hdf5lib%>> cgnsBuild.defs
echo CGNS_SZIPLIB = %sziplib%>> cgnsBuild.defs
echo CGNS_ZLIBLIB = %zliblib%>> cgnsBuild.defs
echo.>> cgnsBuild.defs
echo #------------------------------------------------------------------------>> cgnsBuild.defs
echo # CGNS_MPIINC  - path to MPI header files>> cgnsBuild.defs
echo # CGNS_MPILIBS - MPI libraries>> cgnsBuild.defs
echo #------------------------------------------------------------------------>> cgnsBuild.defs
echo.>> cgnsBuild.defs
echo CGNS_MPIINC  = %mpiinc%>> cgnsBuild.defs
echo CGNS_MPILIBS = %mpilibs%>> cgnsBuild.defs
echo.>> cgnsBuild.defs
echo #------------------------------------------------------------------------>> cgnsBuild.defs
echo # CGNS_LINKLIBS contains the list of libraries>> cgnsBuild.defs
echo #               with which a CGNS application needs to link>> cgnsBuild.defs
echo #------------------------------------------------------------------------>> cgnsBuild.defs
echo.>> cgnsBuild.defs
echo CGNS_LINKLIBS = $(CGNS_LIB) $(CGNS_HDF5LIB) $(CGNS_SZIPLIB) \>> cgnsBuild.defs
echo 	$(CGNS_ZLIBLIB) $(CGNS_MPILIBS)>> cgnsBuild.defs

rem ----- create Makefile

echo creating Makefile
echo # nmake makefile for the CGNS Library under Windows> Makefile
echo.>> Makefile
echo !include make.defs>> Makefile
echo.>> Makefile
echo .SUFFIXES :>> Makefile
echo .SUFFIXES : .c .$(O) $(EXE)>> Makefile
echo.>> Makefile
echo OBJDIR  = %builddir%>> Makefile
echo CGNSLIB = $(LIBCGNS)>> Makefile
echo INSTLIB = cgns.$(A)>> Makefile
if %target% == dll (
  echo CGNSDLL = $^(OBJDIR^)\cgns.dll>> Makefile
  echo INSTDLL = cgns.dll>> Makefile
)
echo.>> Makefile
echo COPTS   = $(CFLAGS) $(CFGFLAGS) /I. %build% %dllopts%>> Makefile
if not "%hdf5inc%" == "" (
  echo # uncomment the following when using HDF5 DLL>> Makefile
  echo %hdf5dll% = /DWIN32 /D_HDF5USEDLL_>> Makefile
)
echo.>> Makefile
echo #---------->> Makefile
echo.>> Makefile
echo CGNSOBJS=\>> Makefile
echo 	$(OBJDIR)\cgns_error.$(O) \>> Makefile
echo 	$(OBJDIR)\cgns_internals.$(O) \>> Makefile
echo 	$(OBJDIR)\cgns_io.$(O) \>> Makefile
echo 	$(OBJDIR)\cgnslib.$(O)>> Makefile
echo.>> Makefile
echo # ADF/ADFH routines>> Makefile
echo.>> Makefile
echo ADFOBJS=\>>Makefile
if not "%hdf5inc%" == "" echo 	$(OBJDIR)\ADFH.$(O) \>>Makefile
echo 	$(OBJDIR)\ADF_interface.$(O) \>>Makefile
echo 	$(OBJDIR)\ADF_internals.$(O)>> Makefile
echo.>> Makefile
echo F2COBJS= $(OBJDIR)\cg_ftoc.$(O) $(OBJDIR)\cgio_ftoc.$(O)>> Makefile
echo.>> Makefile
echo #---------->> Makefile
echo.>> Makefile
echo all     : %target% tools %cgnstools%>> Makefile
if %target% == dll echo dll     : $^(CGNSDLL^)>> Makefile
echo cgns    : $(CGNSLIB)>> Makefile
echo test    : runtests>> Makefile
if "%cgnstools%" == "cgnstools" (
  echo install : install-cgns install-tools install-cgnstools>> Makefile
) else (
  echo install : install-cgns>> Makefile
)
echo.>> Makefile
echo #---------->> Makefile
echo.>> Makefile
echo $(CGNSLIB) : $(OBJDIR) $(CGNSOBJS) $(ADFOBJS) %f2cobjs%>> Makefile
echo 	-@$(RM) $@>> Makefile
echo 	$(AR) /lib /nologo $(AROUT)$@ $(CGNSOBJS) $(ADFOBJS) %f2cobjs%>> Makefile
if not %target% == dll goto make2
echo.>> Makefile
echo #---------->> Makefile
echo.>> Makefile
echo $(CGNSDLL) : $(OBJDIR) $(CGNSOBJS) $(ADFOBJS) %f2cobjs%>> Makefile
echo 	-@$(RM) $@ $(CGNSLIB)>> Makefile
echo 	$(AR) /dll /nologo $(AROUT)$@ $(CGNSOBJS) $(ADFOBJS) %f2cobjs% $(BUILDLIBS)>> Makefile
:make2
echo.>> Makefile
echo #---------->> Makefile
echo.>> Makefile
echo $(OBJDIR) :>> Makefile
echo 	-$(MKDIR) $(OBJDIR)>> Makefile
echo.>> Makefile
echo #---------->> Makefile
echo.>> Makefile
echo tools : %target%>> Makefile
echo 	-cd tools ^&^& %make%>> Makefile
echo.>> Makefile
echo tests : %target%>> Makefile
echo 	-cd tests ^&^& %make%>> Makefile
echo.>> Makefile
echo cgnstools : %target%>> Makefile
echo 	-cd cgnstools ^&^& %make%>> Makefile
echo.>> Makefile
echo runtests : %target%>> Makefile
echo 	-cd tests ^&^& %make% test>> Makefile
echo.>> Makefile
echo #---------->> Makefile
echo.>> Makefile
echo clean :>> Makefile
echo 	-cd $(OBJDIR) ^&^& $(RM) *.$(O)>> Makefile
echo 	-cd tools ^&^& %make% clean>> Makefile
echo 	-cd tests ^&^& %make% clean>> Makefile
echo 	-cd cgnstools ^&^& %make% clean>> Makefile
echo.>> Makefile
echo allclean : distclean>> Makefile
echo.>> Makefile
echo distclean : clean>> Makefile
echo 	-cd tools ^&^& %make% allclean>> Makefile
echo 	-cd tests ^&^& %make% allclean>> Makefile
echo 	-cd tools ^&^& $(RM) Makefile>> Makefile
echo 	-cd tests ^&^& $(RM) Makefile>> Makefile
echo 	-cd cgnstools ^&^& %make% distclean>> Makefile
echo 	-$(RM) $(CGNSLIB)>> Makefile
if %target% == dll echo 	-$(RM) $(CGNSDLL)>> Makefile
echo 	-$(RMDIR) $(OBJDIR)>> Makefile
echo 	-$(RM) cgnstypes.h cgnstypes_f.h cgnslib_f.h>> Makefile
echo 	-$(RM) *.pdb>> Makefile
echo 	-$(RM) make.defs cgnsBuild.defs Makefile>> Makefile
echo.>> Makefile
echo install-cgns : %target% $(INCLUDEDIR) $(LIBDIR) %adfinc%>> Makefile
echo 	$(INSTALL_DATA) cgnstypes.h $(INCLUDEDIR)\cgnstypes.h>> Makefile
echo 	$(INSTALL_DATA) cgnstypes_f.h $(INCLUDEDIR)\cgnstypes_f.h>> Makefile
echo 	$(INSTALL_DATA) cgnslib.h $(INCLUDEDIR)\cgnslib.h>> Makefile
echo 	$(INSTALL_DATA) cgnslib_f.h $(INCLUDEDIR)\cgnslib_f.h>> Makefile
echo 	$(INSTALL_DATA) cgnswin_f.h $(INCLUDEDIR)\cgnswin_f.h>> Makefile
echo 	$(INSTALL_DATA) cgns_io.h $(INCLUDEDIR)\cgns_io.h>> Makefile
echo 	$(INSTALL_DATA) cgnsBuild.defs $(INCLUDEDIR)\cgnsBuild.defs>> Makefile
echo 	$(INSTALL_DATA) $(CGNSLIB) $(LIBDIR)\$(INSTLIB)>> Makefile
if %target% == dll echo 	$(INSTALL_DATA) $(CGNSDLL) $(LIBDIR)\$(INSTDLL)>> Makefile
echo.>> Makefile
echo install-all : install>> Makefile
echo.>> Makefile
echo install-tools :>> Makefile
echo 	-cd tools ^&^& %make% install>> Makefile
echo.>> Makefile
echo install-cgnstools :>> Makefile
echo 	-cd cgnstools ^&^& %make% install>> Makefile
echo.>> Makefile
echo $(INCLUDEDIR) : $(INSTALLDIR)>> Makefile
echo 	-$(MKDIR) $(INCLUDEDIR)>> Makefile
echo.>> Makefile
echo $(LIBDIR) : $(INSTALLDIR)>> Makefile
echo 	-$(MKDIR) $(LIBDIR)>> Makefile
echo.>> Makefile
echo $(INSTALLDIR) :>> Makefile
echo 	-$(MKDIR) $(INSTALLDIR)>> Makefile
echo.>> Makefile
echo install-adf : $(INCLUDEDIR)\adf>> Makefile
echo 	$(INSTALL_DATA) adf\ADF.h $(INCLUDEDIR)\adf\ADF.h>> Makefile
echo.>> Makefile
echo $(INCLUDEDIR)\adf : $(INCLUDEDIR)>> Makefile
echo 	-$(MKDIR) $(INCLUDEDIR)\adf>> Makefile
echo.>> Makefile
echo install-adfh : $(INCLUDEDIR)\adfh>> Makefile
echo 	$(INSTALL_DATA) adfh\ADFH.h $(INCLUDEDIR)\adfh\ADFH.h>> Makefile
echo.>> Makefile
echo $(INCLUDEDIR)\adfh : $(INCLUDEDIR)>> Makefile
echo 	-$(MKDIR) $(INCLUDEDIR)\adfh>> Makefile
echo.>> Makefile
echo #---------- mid-level library>> Makefile
echo.>> Makefile
echo $(OBJDIR)\cgns_error.$(O) : cgns_error.c cgnslib.h cgns_header.h cgns_io.h>> Makefile
echo 	$(CC) $(COPTS) $(COOUT)$@ /c cgns_error.c>> Makefile
echo.>> Makefile
echo $(OBJDIR)\cgns_internals.$(O) : cgns_internals.c cgnslib.h cgns_header.h cgns_io.h>> Makefile
echo 	$(CC) $(COPTS) $(COOUT)$@ /c cgns_internals.c>> Makefile
echo.>> Makefile
echo $(OBJDIR)\cgns_io.$(O) : cgns_io.c cgnslib.h cgns_io.h \>> Makefile
set includes=adf\ADF.h
if not "%hdf5inc%" == "" set includes=%includes% adfh\ADFH.h
echo 	%includes%>> Makefile
echo 	$(CC) $(COPTS) $(COOUT)$@ /c cgns_io.c>> Makefile
echo.>> Makefile
echo $(OBJDIR)\cgnslib.$(O) : cgnslib.c cgnslib.h cgns_header.h cgns_io.h>> Makefile
echo 	$(CC) $(COPTS) $(HDF5INC) $(MPIINC) $(COOUT)$@ /c cgnslib.c>> Makefile
echo.>> Makefile
echo $(OBJDIR)\cg_ftoc.$(O) : cg_ftoc.c fortran_macros.h cgnslib.h cgns_header.h cgns_io.h>> Makefile
echo 	$(CC) $(COPTS) $(F2CFLAGS) $(COOUT)$@ /c cg_ftoc.c>> Makefile
echo $(OBJDIR)\cgio_ftoc.$(O) : cgio_ftoc.c fortran_macros.h cgns_io.h>> Makefile
echo 	$(CC) $(COPTS) $(F2CFLAGS) $(COOUT)$@ /c cgio_ftoc.c>> Makefile
echo.>> Makefile
echo cgnslib.h : cgnstypes.h>> Makefile
echo cgns_header.h : cgnstypes.h>> Makefile
echo cgns_io.h : cgnstypes.h>> Makefile
echo.>> Makefile
echo #---------- ADF>> Makefile
echo.>> Makefile
echo $(OBJDIR)\ADF_interface.$(O) : adf\ADF_interface.c \>> Makefile
echo 	adf\ADF.h adf\ADF_internals.h>> Makefile
echo 	$(CC) $(COPTS) /Iadf $(COOUT)$@ /c adf\ADF_interface.c>> Makefile
echo.>> Makefile
echo $(OBJDIR)\ADF_internals.$(O) : adf\ADF_internals.c \>> Makefile
echo 	adf\ADF.h adf\ADF_internals.h>> Makefile
echo 	$(CC) $(COPTS) /Iadf $(COOUT)$@ /c adf\ADF_internals.c>> Makefile
echo.>> Makefile
echo adf\ADF.h : cgnstypes.h>> Makefile
echo adf\ADF_internals.h : cgnstypes.h>> Makefile
echo.>> Makefile
echo #---------- HDF5>> Makefile
echo.>> Makefile
echo $(OBJDIR)\ADFH.$(O) : adfh\ADFH.c adfh\ADFH.h>> Makefile
echo 	$(CC) $(COPTS) /Iadfh $(HDF5INC) $(HDF5DLL) $(MPIINC) $(COOUT)$@ /c adfh\ADFH.c>> Makefile
echo.>> Makefile
echo adfh\ADFH.h : cgnstypes.h>> Makefile

rem ----- create tools/Makefile

if not exist tools\nul goto tests

echo creating tools\Makefile
echo # nmake makefile for Windows> tools\Makefile
echo.>> tools\Makefile
echo CGNSDIR = ..>> tools\Makefile
echo !include $(CGNSDIR)\make.defs>> tools\Makefile
echo.>> tools\Makefile
echo CGNSLIB = $(CGNSDIR)\$(LIBCGNS)>> tools\Makefile
echo.>> tools\Makefile
echo COPTS  = $(CFLAGS) /I$(CGNSDIR)>> tools\Makefile
echo LDLIBS = $(CGNSLIB) $(BUILDLIBS)>> tools\Makefile
echo.>> tools\Makefile
echo #---------->> tools\Makefile
echo.>> tools\Makefile
echo ALL =	cgnslist$(EXE) \>> tools\Makefile
echo 	cgnscheck$(EXE) \>> tools\Makefile
echo 	cgnsversion$(EXE) \>> tools\Makefile
echo 	cgnsconvert$(EXE) \>> tools\Makefile
echo 	cgnscompress$(EXE) \>> tools\Makefile
echo 	cgnsdiff$(EXE) \>> tools\Makefile
echo 	cgnsnames$(EXE)>> tools\Makefile
echo.>> tools\Makefile
echo all : $(ALL)>> tools\Makefile
echo.>> tools\Makefile
echo #---------->> tools\Makefile
echo.>> tools\Makefile
echo cgnslist$(EXE) : cgnslist.$(O) getargs.$(O) $(CGNSLIB)>> tools\Makefile
echo 	$(CC) $(COPTS) $(CEOUT)$@ cgnslist.$(O) getargs.$(O) $(LDLIBS) $(CLIBS)>> tools\Makefile
echo cgnslist.$(O) : cgnslist.c getargs.h>> tools\Makefile
echo 	$(CC) $(COPTS) /c cgnslist.c>> tools\Makefile
echo.>> tools\Makefile
echo #---------->> tools\Makefile
echo.>> tools\Makefile
echo cgnscheck$(EXE) : cgnscheck.$(O) getargs.$(O) hash.$(O) cgnames.$(O) $(CGNSLIB)>> tools\Makefile
echo 	$(CC) $(COPTS) $(CEOUT)$@ cgnscheck.$(O) getargs.$(O) hash.$(O) cgnames.$(O) $(LDLIBS) $(CLIBS)>> tools\Makefile
echo cgnscheck.$(O) : cgnscheck.c getargs.h hash.h cgnames.h>> tools\Makefile
echo 	$(CC) $(COPTS) /c cgnscheck.c>> tools\Makefile
echo.>> tools\Makefile
echo #---------->> tools\Makefile
echo.>> tools\Makefile
echo cgnsversion$(EXE) : cgnsversion.$(O) getargs.$(O) $(CGNSLIB)>> tools\Makefile
echo 	$(CC) $(CFLAGS) $(CEOUT)$@ cgnsversion.$(O) getargs.$(O) $(LDLIBS) $(CLIBS)>> tools\Makefile
echo cgnsversion.$(O) : cgnsversion.c getargs.h>> tools\Makefile
echo 	$(CC) $(COPTS) /c cgnsversion.c>> tools\Makefile
echo.>> tools\Makefile
echo #---------->> tools\Makefile
echo.>> tools\Makefile
echo cgnsconvert$(EXE) : cgnsconvert.$(O) getargs.$(O) $(CGNSLIB)>> tools\Makefile
echo 	$(CC) $(CFLAGS) $(CEOUT)$@ cgnsconvert.$(O) getargs.$(O) $(LDLIBS) $(CLIBS)>> tools\Makefile
echo cgnsconvert.$(O) : cgnsconvert.c getargs.h>> tools\Makefile
echo 	$(CC) $(COPTS) /c cgnsconvert.c>> tools\Makefile
echo.>> tools\Makefile
echo #---------->> tools\Makefile
echo.>> tools\Makefile
echo cgnscompress$(EXE) : cgnscompress.$(O) $(CGNSLIB)>> tools\Makefile
echo 	$(CC) $(CFLAGS) $(CEOUT)$@ cgnscompress.$(O) $(LDLIBS) $(CLIBS)>> tools\Makefile
echo cgnscompress.$(O) : cgnscompress.c>> tools\Makefile
echo 	$(CC) $(COPTS) /c cgnscompress.c>> tools\Makefile
echo.>> tools\Makefile
echo #---------->> tools\Makefile
echo.>> tools\Makefile
echo cgnsdiff$(EXE) : cgnsdiff.$(O) getargs.$(O) $(CGNSLIB)>> tools\Makefile
echo 	$(CC) $(CFLAGS) $(CEOUT)$@ cgnsdiff.$(O) getargs.$(O) $(LDLIBS) $(CLIBS)>> tools\Makefile
echo cgnsdiff.$(O) : cgnsdiff.c getargs.h>> tools\Makefile
echo 	$(CC) $(COPTS) /c cgnsdiff.c>> tools\Makefile
echo.>> tools\Makefile
echo #---------->> tools\Makefile
echo.>> tools\Makefile
echo cgnsnames$(EXE) : cgnsnames.$(O) cgnames.$(O) $(CGNSLIB)>> tools\Makefile
echo 	$(CC) $(CFLAGS) $(CEOUT)$@ cgnsnames.$(O) cgnames.$(O) $(LDLIBS) $(CLIBS)>> tools\Makefile
echo cgnsnames.$(O) : cgnsnames.c cgnames.h>> tools\Makefile
echo 	$(CC) $(COPTS) /c cgnsnames.c>> tools\Makefile
echo.>> tools\Makefile
echo #---------->> tools\Makefile
echo.>> tools\Makefile
echo getargs.$(O) : getargs.c getargs.h>> tools\Makefile
echo 	$(CC) $(COPTS) /c getargs.c>> tools\Makefile
echo.>> tools\Makefile
echo hash.$(O) : hash.c hash.h>> tools\Makefile
echo 	$(CC) $(COPTS) /c hash.c>> tools\Makefile
echo.>> tools\Makefile
echo cgnames.$(O) : cgnames.c cgnames.h>> tools\Makefile
echo 	$(CC) $(COPTS) /c cgnames.c>> tools\Makefile
echo.>> tools\Makefile
echo install : all $(BINDIR)>> tools\Makefile
echo 	$(INSTALL_PROG) cgnslist$(EXE) $(BINDIR)>> tools\Makefile
echo 	$(INSTALL_PROG) cgnscheck$(EXE) $(BINDIR)>> tools\Makefile
echo 	$(INSTALL_PROG) cgnsversion$(EXE) $(BINDIR)>> tools\Makefile
echo 	$(INSTALL_PROG) cgnsconvert$(EXE) $(BINDIR)>> tools\Makefile
echo 	$(INSTALL_PROG) cgnscompress$(EXE) $(BINDIR)>> tools\Makefile
echo 	$(INSTALL_PROG) cgnsdiff$(EXE) $(BINDIR)>> tools\Makefile
echo 	$(INSTALL_PROG) cgnsnames$(EXE) $(BINDIR)>> tools\Makefile
echo 	$(INSTALL_PROG) adf2hdf.bat $(BINDIR)>> tools\Makefile
echo 	$(INSTALL_PROG) hdf2adf.bat $(BINDIR)>> tools\Makefile
echo 	$(INSTALL_PROG) cgnsupdate.bat $(BINDIR)>> tools\Makefile
echo.>> tools\Makefile
echo $(BINDIR) : $(INSTALLDIR)>> tools\Makefile
echo 	-$(MKDIR) $(BINDIR)>> tools\Makefile
echo.>> tools\Makefile
echo $(INSTALLDIR) :>> tools\Makefile
echo 	-$(MKDIR) $(INSTALLDIR)>> tools\Makefile
echo.>> tools\Makefile
echo clean :>> tools\Makefile
echo 	-$(RM) *.$(O)>> tools\Makefile
echo.>> tools\Makefile
echo allclean : clean>> tools\Makefile
echo 	-$(RM) *.exe>> tools\Makefile
echo 	-$(RM) *.pdb *.ilk>> tools\Makefile

rem ----- create tests/Makefile

:tests
if not exist tests\nul goto cgnstools
echo creating tests\Makefile
echo # nmake makefile for Windows> tests\Makefile
echo.>> tests\Makefile
echo CGNSDIR = ..>> tests\Makefile
echo !include $(CGNSDIR)\make.defs>> tests\Makefile
echo.>> tests\Makefile
echo CGNSLIB = $(CGNSDIR)\$(LIBCGNS)>> tests\Makefile
echo.>> tests\Makefile
echo COPTS  = $(CFLAGS) /I$(CGNSDIR)>> tests\Makefile
echo FOPTS  = $(FFLAGS) /I$(CGNSDIR)>> tests\Makefile
echo LDLIBS = $(CGNSLIB) $(BUILDLIBS)>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo CALL = \>> tests\Makefile
echo 	elemtest$(EXE) \>> tests\Makefile
echo 	test_exts$(EXE) \>> tests\Makefile
echo 	test_partial$(EXE) \>> tests\Makefile
echo 	test_goto$(EXE) \>> tests\Makefile
echo 	test_ver31$(EXE) \>> tests\Makefile
echo 	write_array$(EXE) \>> tests\Makefile
echo 	write_links$(EXE) \>> tests\Makefile
echo 	write_bcdata$(EXE) \>> tests\Makefile
echo 	write_test$(EXE) \>> tests\Makefile
echo 	write_zones$(EXE) \>> tests\Makefile
echo 	write_rind$(EXE)>> tests\Makefile
echo FALL =	cgwrite$^(EXE^) \>> tests\Makefile
echo 	cgread$^(EXE^) \>> tests\Makefile
echo 	cgzconn$(EXE) \>> tests\Makefile
echo 	cgsubreg$(EXE)>> tests\Makefile
echo CALL64 = test64c$^(EXE^)>> tests\Makefile
echo FALL64 = test64f$^(EXE^)>> tests\Makefile
echo.>> tests\Makefile
if not "%f2c%" == "none" (
  echo TESTS = $^(CALL^) $^(FALL^)>> tests\Makefile
  echo ALL64 = $^(CALL64^) $^(FALL64^)>> tests\Makefile
) else (
  echo TESTS = $^(CALL^)>> tests\Makefile
  echo ALL64 = $^(CALL64^)>> tests\Makefile
)
echo ALL   = dbtest$(EXE) open_cgns$(EXE) $(TESTS)>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo all : $(ALL)>> tests\Makefile
echo fortran : $(FALL)>> tests\Makefile
echo test64 : $(ALL64)>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo test : $(TESTS)>> tests\Makefile
echo 	@echo === running tests ===>> tests\Makefile
echo 	-@runtest.bat elemtest$(EXE)>> tests\Makefile
echo 	-@runtest.bat test_exts$(EXE)>> tests\Makefile
echo 	-@runtest.bat test_partial$(EXE)>> tests\Makefile
echo 	-@runtest.bat test_goto$(EXE)>> tests\Makefile
echo 	-@runtest.bat test_ver31$(EXE)>> tests\Makefile
echo 	-@runtest.bat write_array$(EXE)>> tests\Makefile
echo 	-@runtest.bat write_links$(EXE)>> tests\Makefile
echo 	-@runtest.bat write_bcdata$(EXE)>> tests\Makefile
echo 	-@runtest.bat write_test$(EXE)>> tests\Makefile
echo 	-@runtest.bat write_zones$(EXE)>> tests\Makefile
echo 	-@runtest.bat write_rind$(EXE)>> tests\Makefile
if not "%f2c%" == "none" (
  echo 	-@runtest.bat cgwrite$^(EXE^)>> tests\Makefile
  echo 	-@runtest.bat cgread$^(EXE^)>> tests\Makefile
  echo 	-@runtest.bat cgzconn$^(EXE^)>> tests\Makefile
  echo 	-@runtest.bat cgsubreg$^(EXE^)>> tests\Makefile
)
echo 	@echo === finished ===>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo dbtest$(EXE) : dbtest.$(O) utils.$(O) $(CGNSLIB)>> tests\Makefile
echo 	$(CC) $(COPTS) $(CEOUT)$@ dbtest.$(O) utils.$(O) $(LDLIBS) $(CLIBS)>> tests\Makefile
echo dbtest.$(O) : dbtest.c utils.h>> tests\Makefile
echo 	$(CC) $(COPTS) /c dbtest.c>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo elemtest$(EXE) : elemtest.c $(CGNSLIB)>> tests\Makefile
echo 	$(CC) $(COPTS) $(CEOUT)$@ elemtest.c $(LDLIBS) $(CLIBS)>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo open_cgns$(EXE) : open_cgns.$(O) utils.$(O) $(CGNSLIB)>> tests\Makefile
echo 	$(CC) $(CFLAGS) $(CEOUT)$@ open_cgns.$(O) utils.$(O) $(LDLIBS) $(CLIBS)>> tests\Makefile
echo open_cgns.$(O) : open_cgns.c utils.h>> tests\Makefile
echo 	$(CC) $(COPTS) /c open_cgns.c>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo test_exts$(EXE) : test_exts.c $(CGNSLIB)>> tests\Makefile
echo 	$(CC) $(COPTS) $(CEOUT)$@ test_exts.c $(LDLIBS) $(CLIBS)>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo test_partial$(EXE) : test_partial.c $(CGNSLIB)>> tests\Makefile
echo 	$(CC) $(COPTS) $(CEOUT)$@ test_partial.c $(LDLIBS) $(CLIBS)>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo test_goto$(EXE) : test_goto.c $(CGNSLIB)>> tests\Makefile
echo 	$(CC) $(COPTS) $(CEOUT)$@ test_goto.c $(LDLIBS) $(CLIBS)>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo test_ver31$(EXE) : test_ver31.c $(CGNSLIB)>> tests\Makefile
echo 	$(CC) $(COPTS) $(CEOUT)$@ test_ver31.c $(LDLIBS) $(CLIBS)>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo write_array$(EXE) : write_array.$(O) utils.$(O) $(CGNSLIB)>> tests\Makefile
echo 	$(CC) $(CFLAGS) $(CEOUT)$@ write_array.$(O) utils.$(O) $(LDLIBS) $(CLIBS)>> tests\Makefile
echo write_array.$(O) : write_array.c utils.h>> tests\Makefile
echo 	$(CC) $(COPTS) /c write_array.c>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo write_links$(EXE) : write_links.$(O) utils.$(O) $(CGNSLIB)>> tests\Makefile
echo 	$(CC) $(CFLAGS) $(CEOUT)$@ write_links.$(O) utils.$(O) $(LDLIBS) $(CLIBS)>> tests\Makefile
echo write_links.$(O) : write_links.c utils.h>> tests\Makefile
echo 	$(CC) $(COPTS) /c write_links.c>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo write_bcdata$(EXE) : write_bcdata.$(O) utils.$(O) $(CGNSLIB)>> tests\Makefile
echo 	$(CC) $(CFLAGS) $(CEOUT)$@ write_bcdata.$(O) utils.$(O) $(LDLIBS) $(CLIBS)>> tests\Makefile
echo write_bcdata.$(O) : write_bcdata.c utils.h>> tests\Makefile
echo 	$(CC) $(COPTS) /c write_bcdata.c>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo write_test$(EXE) : write_test.c $(CGNSLIB)>> tests\Makefile
echo 	$(CC) $(COPTS) $(CEOUT)$@ write_test.c $(LDLIBS) $(CLIBS)>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo write_zones$(EXE) : write_zones.$(O) utils.$(O) $(CGNSLIB)>> tests\Makefile
echo 	$(CC) $(CFLAGS) $(CEOUT)$@ write_zones.$(O) utils.$(O) $(LDLIBS) $(CLIBS)>> tests\Makefile
echo write_zones.$(O) : write_zones.c utils.h>> tests\Makefile
echo 	$(CC) $(COPTS) /c write_zones.c>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo write_rind$(EXE) : write_rind.c $(CGNSLIB)>> tests\Makefile
echo 	$(CC) $(COPTS) $(CEOUT)$@ write_rind.c $(LDLIBS) $(CLIBS)>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo cgwrite$(EXE) : cgwrite.F $(CGNSLIB)>> tests\Makefile
echo 	$(F77) $(FOPTS) $(FEOUT)$@ cgwrite.F $(LDLIBS) $(FLIBS)>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo cgread$(EXE) : cgread.F $(CGNSLIB)>> tests\Makefile
echo 	$(F77) $(FOPTS) $(FEOUT)$@ cgread.F $(LDLIBS) $(FLIBS)>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo cgzconn$(EXE) : cgzconn.F $(CGNSLIB)>> tests\Makefile
echo 	$(F77) $(FOPTS) $(FEOUT)$@ cgzconn.F $(LDLIBS) $(FLIBS)>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo cgsubreg$(EXE) : cgsubreg.F $(CGNSLIB)>> tests\Makefile
echo 	$(F77) $(FOPTS) $(FEOUT)$@ cgsubreg.F $(LDLIBS) $(FLIBS)>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo test64c$(EXE) : test64c.$(O) utils.$(O) $(CGNSLIB)>> tests\Makefile
echo 	$(CC) $(CFLAGS) $(CEOUT)$@ test64c.$(O) utils.$(O) $(LDLIBS) $(CLIBS)>> tests\Makefile
echo 	$(STRIP) $@>> tests\Makefile
echo test64c.$(O) : test64c.c utils.h>> tests\Makefile
echo 	$(CC) $(COPTS) /c test64c.c>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo test64f$(EXE) : test64f.F $(CGNSLIB)>> tests\Makefile
echo 	$(F77) $(FOPTS) $(FEOUT)$@ test64f.F $(LDLIBS) $(FLIBS)>> tests\Makefile
echo 	$(STRIP) $@>> tests\Makefile
echo.>> tests\Makefile
echo #---------->> tests\Makefile
echo.>> tests\Makefile
echo utils.$(O) : utils.c utils.h>> tests\Makefile
echo 	$(CC) $(COPTS) /c utils.c>> tests\Makefile
echo.>> tests\Makefile
echo clean :>> tests\Makefile
echo 	-$(RM) *.$(O)>> tests\Makefile
echo.>> tests\Makefile
echo allclean : clean>> tests\Makefile
echo 	-$(RM) *.exe>> tests\Makefile
echo 	-$(RM) *.pdb *.ilk>> tests\Makefile
echo 	-$(RM) *.cgns *.cgio>> tests\Makefile

rem ----- create cgsntools\Makefile

:cgnstools
if not exist cgnstools\nul goto done

echo creating cgsntools\Makefile
echo # nmake makefile for the CGNS tools under Windows> cgnstools\Makefile
echo.>> cgnstools\Makefile
echo DOMAKE = nmake /nologo /f Makefile.win>> cgnstools\Makefile
echo !include ..\make.defs>> cgnstools\Makefile
echo.>> cgnstools\Makefile
echo defaults : cgnsview cgnscalc cgnsplot utilities>> cgnstools\Makefile
echo install : install-config install-cgnsview install-cgnscalc \>> cgnstools\Makefile
echo 	install-cgnsplot install-utilities>> cgnstools\Makefile
echo.>> cgnstools\Makefile
echo all : cgnsview cgnscalc cgnsplot utilities>> cgnstools\Makefile
echo install-all : install-config install-cgnsview install-cgnscalc \>> cgnstools\Makefile
echo 	install-cgnsplot install-utilities>> cgnstools\Makefile
echo.>> cgnstools\Makefile
echo clean :>> cgnstools\Makefile
echo 	cd cgnsview ^&^& $(DOMAKE) clean>> cgnstools\Makefile
echo 	cd cgnscalc ^&^& $(DOMAKE) clean>> cgnstools\Makefile
echo 	cd cgnsplot ^&^& $(DOMAKE) clean>> cgnstools\Makefile
echo 	cd utilities ^&^& $(DOMAKE) clean>> cgnstools\Makefile
echo 	cd calclib ^&^& $(DOMAKE) clean>> cgnstools\Makefile
echo 	cd tkogl ^&^& $(DOMAKE) clean>> cgnstools\Makefile
echo.>> cgnstools\Makefile
echo distclean : clean>> cgnstools\Makefile
echo 	-del make.win Makefile cgconfig.bat>> cgnstools\Makefile
echo.>> cgnstools\Makefile
echo cgnsview  : prog-cgnsview>> cgnstools\Makefile
echo cgnscalc  : prog-cgnscalc>> cgnstools\Makefile
echo cgnsplot  : prog-cgnsplot>> cgnstools\Makefile
echo utilities : prog-utilities>> cgnstools\Makefile
echo.>> cgnstools\Makefile
echo prog-cgnsview :>> cgnstools\Makefile
echo 	cd cgnsview ^&^& $(DOMAKE)>> cgnstools\Makefile
echo.>> cgnstools\Makefile
echo prog-cgnscalc : lib-calclib>> cgnstools\Makefile
echo 	cd cgnscalc ^&^& $(DOMAKE)>> cgnstools\Makefile
echo.>> cgnstools\Makefile
echo prog-cgnsplot : lib-tkogl>> cgnstools\Makefile
echo 	cd cgnsplot ^&^& $(DOMAKE)>> cgnstools\Makefile
echo.>> cgnstools\Makefile
echo prog-utilities : lib-calclib>> cgnstools\Makefile
echo 	cd utilities ^&^& $(DOMAKE)>> cgnstools\Makefile
echo.>> cgnstools\Makefile
echo lib-calclib :>> cgnstools\Makefile
echo 	cd calclib ^&^& $(DOMAKE)>> cgnstools\Makefile
echo.>> cgnstools\Makefile
echo lib-tkogl :>> cgnstools\Makefile
echo 	cd tkogl ^&^& $(DOMAKE)>> cgnstools\Makefile
echo.>> cgnstools\Makefile
echo install-config : $(BINDIR)>> cgnstools\Makefile
echo 	$(INSTALL) cgconfig.bat $(BINDIR)>> cgnstools\Makefile
echo.>> cgnstools\Makefile
echo $(BINDIR) : $(INSTALLDIR)>> cgnstools\Makefile
echo 	-mkdir $(BINDIR)>> cgnstools\Makefile
echo.>> cgnstools\Makefile
echo $(INSTALLDIR) :>> cgnstools\Makefile
echo 	-mkdir $(INSTALLDIR)>> cgnstools\Makefile
echo.>> cgnstools\Makefile
echo install-cgnsview :>> cgnstools\Makefile
echo 	cd cgnsview ^&^& $(DOMAKE) install>> cgnstools\Makefile
echo.>> cgnstools\Makefile
echo install-cgnscalc : lib-calclib>> cgnstools\Makefile
echo 	cd cgnscalc ^&^& $(DOMAKE) install>> cgnstools\Makefile
echo.>> cgnstools\Makefile
echo install-cgnsplot : lib-tkogl>> cgnstools\Makefile
echo 	cd cgnsplot ^&^& $(DOMAKE) install>> cgnstools\Makefile
echo.>> cgnstools\Makefile
echo install-utilities : lib-calclib>> cgnstools\Makefile
echo 	cd utilities ^&^& $(DOMAKE) install>> cgnstools\Makefile

rem ----- create cgnstools\make.win

echo creating cgnstools\make.win
echo # makefile include for %windir%> cgnstools\make.win
echo.>> cgnstools\make.win
echo #------------------------------------------------------->> cgnstools\make.win
echo # CGNS setup>> cgnstools\make.win
echo #------------------------------------------------------->> cgnstools\make.win
echo.>> cgnstools\make.win
echo CGNSDIR  = ..\..>> cgnstools\make.win
echo !include $(CGNSDIR)\make.defs>> cgnstools\make.win
echo CGNSLIB  = $(CGNSDIR)\$(LIBCGNS)>> cgnstools\make.win
echo.>> cgnstools\make.win
echo SHAREDIR = $(INSTALLDIR)\share>> cgnstools\make.win
echo.>> cgnstools\make.win
echo #------------------------------------------------------->> cgnstools\make.win
echo # linker options>> cgnstools\make.win
echo #------------------------------------------------------->> cgnstools\make.win
echo.>> cgnstools\make.win
echo LINK   = link>> cgnstools\make.win
echo LFLAGS = /nologo %lopts%>> cgnstools\make.win
echo.>> cgnstools\make.win
echo #------------------------------------------------------->> cgnstools\make.win
echo # path to the standard Tcl/Tk includes and libraries>> cgnstools\make.win
echo # the include path needs to include tcl.h, tk.h,>> cgnstools\make.win
echo # tkWinInt.h and the X11 include subdirectory>> cgnstools\make.win
echo #------------------------------------------------------->> cgnstools\make.win
echo.>> cgnstools\make.win
echo TCLINC = %tclinc%>> cgnstools\make.win
echo TCLLIB = %tcllib%>> cgnstools\make.win
echo TKLIB  = %tklib%>> cgnstools\make.win
echo.>> cgnstools\make.win
echo #------------------------------------------------------->> cgnstools\make.win
echo # TKOGLINCS give the include directories>> cgnstools\make.win
echo # TKOGLLIB is the library relative to cgnsplot directory>> cgnstools\make.win
echo #------------------------------------------------------->> cgnstools\make.win
echo.>> cgnstools\make.win
echo TKOGLINCS = $(TCLINC)>> cgnstools\make.win
echo TKOGLLIB  = ..\tkogl\tkogl.lib>> cgnstools\make.win
echo.>> cgnstools\make.win
echo #---------------------------------------------------------->> cgnstools\make.win
echo # compile options for cgnsplot>> cgnstools\make.win
echo #    -DNO_MESH_BOUNDARIES - no structured mesh boundaries>> cgnstools\make.win
echo #    -DNO_CUTTING_PLANE   - no cutting plane>> cgnstools\make.win
echo #---------------------------------------------------------->> cgnstools\make.win
echo.>> cgnstools\make.win
echo PLOTOPTS = %plotopts%>> cgnstools\make.win
echo.>> cgnstools\make.win
echo #------------------------------------------------------->> cgnstools\make.win
echo # set to trap math errors in calculator>> cgnstools\make.win
echo #------------------------------------------------------->> cgnstools\make.win
echo.>> cgnstools\make.win
echo MATHERR = /DUSE_MATHERR>> cgnstools\make.win
echo.>> cgnstools\make.win
echo #------------------------------------------------------->> cgnstools\make.win
echo # windows libraries>> cgnstools\make.win
echo #------------------------------------------------------->> cgnstools\make.win
echo.>> cgnstools\make.win
echo dlllibs = gdi32.lib comdlg32.lib>> cgnstools\make.win
echo.>> cgnstools\make.win
echo guilibs	= %mslib% oldnames.lib kernel32.lib advapi32.lib \>> cgnstools\make.win
echo 	user32.lib gdi32.lib comdlg32.lib winspool.lib>> cgnstools\make.win
echo.>> cgnstools\make.win
echo ogllibs = opengl32.lib glu32.lib>> cgnstools\make.win
if not "%winhtml%" == "" (
echo.>> cgnstools\make.win
echo #---------------------------------------------------------->> cgnstools\make.win
echo # build tools with HTMLhelp>> cgnstools\make.win
echo #---------------------------------------------------------->> cgnstools\make.win
echo.>> cgnstools\make.win
echo WINHTML_OPT = /DUSE_HTMLHELP>> cgnstools\make.win
echo WINHTML_OBJ = winhtml.obj>> cgnstools\make.win
echo WINHTML_INC = /I%winhtml%\include>> cgnstools\make.win
echo WINHTML_LIB = %winhtml%\lib\htmlhelp.lib>> cgnstools\make.win
)

rem ----- create cgconfig.bat

echo creating cgnstools\cgconfig.bat
echo rem configuration file for windows> cgnstools\cgconfig.bat
echo.>> cgnstools\cgconfig.bat
echo set CG_BIN_DIR=%instdir%\bin>> cgnstools\cgconfig.bat
echo set CG_LIB_DIR=%instdir%\share>> cgnstools\cgconfig.bat
if %target% == dll echo set PATH=%instdir%\lib;%%PATH%%>> cgnstools\cgconfig.bat
if exist %tcldir%\bin\nul (
echo set PATH=%tcldir%\bin;%%PATH%%>> cgnstools\cgconfig.bat
) else (
echo rem may need to add location of tcl/tk dlls to PATH>> cgnstools\cgconfig.bat
echo rem PATH=%tcldir%\bin;%%PATH%%>> cgnstools\cgconfig.bat
)
echo.>> cgnstools\cgconfig.bat
echo rem may need to set these if system can't find Tcl and TK scripts>> cgnstools\cgconfig.bat
echo rem set TCL_LIBRARY=c:/lib/tcl8.4.13/library>> cgnstools\cgconfig.bat
echo rem set TK_LIBRARY=c:/lib/tk8.4.13/library>> cgnstools\cgconfig.bat

:done
endlocal

