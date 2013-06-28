@echo off
setlocal

set config=0
set help=0
set system=WIN32

:next
if "%1" == "" goto doit
if "%1" == "conf" set config=1
if "%1" == "-conf" set config=1
if "%1" == "configure" set config=1
if "%1" == "-configure" set config=1
if "%1" == "64" set system=WIN64
if "%1" == "-64" set system=WIN64
if "%1" == "64bit" set system=WIN64
if "%1" == "-64bit" set system=WIN64
if "%1" == "help" set help=1
if "%1" == "-help" set help=1
shift
goto next
:doit

if %help% == 0 goto getsys
echo usage: cgsystem [[-]conf[igure]] [[-]64[bit]]
goto done

:getsys
if %config% == 1 goto config
echo %system%
goto done

:config
echo modifying make.system to set SYSTEM=%systen%
echo SYSTEM = %system%> make.system

:done
endlocal

