@echo off
setlocal
if "%1" == "" (
  echo usage: runtest exename
  goto done
)
set PATH=..\lib;%PATH%
%1 > nul 2>&1
if %ERRORLEVEL% == 0 (
  echo %1 ... passed
) else (
  echo %1 ... *** FAILED ***
)
:done
endlocal

