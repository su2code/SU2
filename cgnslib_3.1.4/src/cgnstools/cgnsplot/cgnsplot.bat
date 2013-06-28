@echo off
setlocal

set dir=%~dps0
if exist %dir%cgconfig.bat (
  call %dir%cgconfig.bat
  goto getwish
)
if exist %dir%..\cgconfig.bat call %dir%..\cgconfig.bat

:getwish
if "%CG_BIN_DIR%" == "" set CG_BIN_DIR=%dir%

if exist %dir%plotwish.exe (
  set wish=%dir%plotwish.exe
  goto getscript
)
if exist %dir%cgnstools\plotwish.exe (
  set wish=%dir%cgnstools\plotwish.exe
  goto getscript
)
echo plotwish.exe not found
pause
goto done

:getscript
if exist %dir%cgnsplot.tcl (
  set script=%dir%cgnsplot.tcl
  goto run
)
if not "%CG_LIB_DIR%" == "" (
  if exist %CG_LIB_DIR%\cgnsplot.tcl (
    set script=%CG_LIB_DIR%\cgnsplot.tcl
    goto run
  )
)
if exist %dir%..\share\cgnsplot.tcl (
  set CG_LIB_DIR=%dir%..\share
  set script=%dir%..\share\cgnsplot.tcl
  goto run
)
echo cgnsplot.tcl not found
pause
goto done

:run
start /b %wish% %script% %1

:done
endlocal
