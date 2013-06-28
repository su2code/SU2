@echo off
setlocal

set dir=%~dps0
if not exist %dir%cgnsconvert.exe goto notfound

rem -- for path setting to dlls
if exist %dir%..\cgconfig.bat call %dir%..\cgconfig.bat
if exist %dir%cgconfig.bat call %dir%cgconfig.bat

set links=
if "%1" == "-links" (
  set links=-l
  shift
)

%dir%cgnsconvert -a %links% %1 %2 %3
goto done

:notfound
echo cgnsconvert.exe not found in %dir%
:done
endlocal

