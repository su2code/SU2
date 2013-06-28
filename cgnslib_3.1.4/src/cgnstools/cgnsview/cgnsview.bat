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

if exist %dir%cgiowish.exe (
  set wish=%dir%cgiowish.exe
  goto getscript
)
if exist %dir%cgnstools\cgiowish.exe (
  set wish=%dir%cgnstools\cgiowish.exe
  goto getscript
)
echo cgiowish.exe not found
pause
goto done

:getscript
if exist %dir%cgnsview.tcl (
  set script=%dir%cgnsview.tcl
  goto run
)
if not "%CG_LIB_DIR%" == "" (
  if exist %CG_LIB_DIR%\cgnsview.tcl (
    set script=%CG_LIB_DIR%\cgnsview.tcl
    goto run
  )
)
if exist %dir%..\share\cgnsview.tcl (
  if "%CG_LIB_DIR%" == "" set CG_LIB_DIR=%dir%..\share
  set script=%dir%..\share\cgnsview.tcl
  goto run
)
echo cgnsview.tcl not found
pause
goto done

:run
start /b %wish% %script% %1

:done
endlocal
