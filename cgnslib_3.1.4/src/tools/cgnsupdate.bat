@echo off
setlocal

set dir=%~dps0
if not exist %dir%cgnscheck.exe goto notfound

rem -- for path setting to dlls
if exist %dir%..\cgconfig.bat call %dir%..\cgconfig.bat
if exist %dir%cgconfig.bat call %dir%cgconfig.bat

%dir%cgnscheck -U %1 %2
goto done

:notfound
echo cgnscheck.exe not found in %dir%
:done
endlocal

