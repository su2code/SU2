@echo off
@setlocal
rem ########################################################################
rem # Pointwise - Proprietary software product of Pointwise, Inc.
rem #           Copyright (c) 1995-2009 Pointwise, Inc.
rem #           All rights reserved.
rem #
rem # Plugin SDK, create plugin project bat-script
rem ########################################################################

if NOT "%PWI_TCLSH%" == "" goto DO_IT

rem --------------------------------------------------------------------
rem Search for tclsh in the user-specified runtime distribution dir...
rem --------------------------------------------------------------------
if NOT "%PWI_TCLSH_SEARCHDIR%" == "" (
    if EXIST %PWI_TCLSH_SEARCHDIR%\win64\bin\tclsh.exe (
        set PWI_TCLSH=%PWI_TCLSH_SEARCHDIR%\win64\bin\tclsh.exe
    ) else if EXIST %PWI_TCLSH_SEARCHDIR%\win64\bin\tclshd.exe (
        set PWI_TCLSH=%PWI_TCLSH_SEARCHDIR%\win64\bin\tclshd.exe
    ) else if EXIST %PWI_TCLSH_SEARCHDIR%\win32\bin\tclsh.exe (
        set PWI_TCLSH=%PWI_TCLSH_SEARCHDIR%\win32\bin\tclsh.exe
    ) else if EXIST %PWI_TCLSH_SEARCHDIR%\win32\bin\tclshd.exe (
        set PWI_TCLSH=%PWI_TCLSH_SEARCHDIR%\win32\bin\tclshd.exe
    ) else (
        echo WARNING: tclsh not found in PWI_TCLSH_SEARCHDIR="%PWI_TCLSH_SEARCHDIR%"
    )
)

if NOT "%PWI_TCLSH%" == "" goto DO_IT

rem --------------------------------------------------------------------
rem Search for tclsh in the user-specified runtime install dir...
rem --------------------------------------------------------------------
if NOT "%PWI_INSTALLDIR%" == "" (
    if EXIST %PWI_INSTALLDIR%\win64\bin\tclsh.exe (
        set PWI_TCLSH=%PWI_INSTALLDIR%\win64\bin\tclsh.exe
    ) else if EXIST %PWI_INSTALLDIR%\win64\bin\tclshd.exe (
        set PWI_TCLSH=%PWI_INSTALLDIR%\win64\bin\tclshd.exe
    ) else if EXIST %PWI_INSTALLDIR%\win32\bin\tclsh.exe (
        set PWI_TCLSH=%PWI_INSTALLDIR%\win32\bin\tclsh.exe
    ) else if EXIST %PWI_INSTALLDIR%\win32\bin\tclshd.exe (
        set PWI_TCLSH=%PWI_INSTALLDIR%\win32\bin\tclshd.exe
    ) else (
        echo WARNING: tclsh not found in PWI_INSTALLDIR="%PWI_INSTALLDIR%"
    )
)

if NOT "%PWI_TCLSH%" == "" goto DO_IT

rem --------------------------------------------------------------------
rem Search for tclsh in the user-specified runtime bin dir...
rem --------------------------------------------------------------------
if NOT "%PWI_BINDIR%" == "" (
    if EXIST %PWI_BINDIR%\win64\bin\tclsh.exe (
        set PWI_TCLSH=%PWI_BINDIR%\tclsh.exe
    ) else if EXIST %PWI_BINDIR%\tclshd.exe (
        set PWI_TCLSH=%PWI_BINDIR%\tclshd.exe
    ) else if EXIST %PWI_BINDIR%\tclsh.exe (
        set PWI_TCLSH=%PWI_BINDIR%\tclsh.exe
    ) else if EXIST %PWI_BINDIR%\tclshd.exe (
        set PWI_TCLSH=%PWI_BINDIR%\tclshd.exe
    ) else (
        echo WARNING: tclsh not found in "%PWI_BINDIR%"
    )
)

rem --------------------------------------------------------------------
rem Could not find anything, anywhere. Use tclsh!
rem --------------------------------------------------------------------
if "%PWI_TCLSH%" == "" set PWI_TCLSH=tclsh


rem --------------------------------------------------------------------
rem attempt the command
rem --------------------------------------------------------------------
:DO_IT
set THECMD=%PWI_TCLSH% mkplugin.tcl %*
echo %THECMD%
%THECMD%

endlocal
