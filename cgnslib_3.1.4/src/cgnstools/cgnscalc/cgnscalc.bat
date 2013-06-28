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

if exist %dir%calcwish.exe (
  set wish=%dir%calcwish.exe
  goto getscript
)
if exist %dir%cgnstools\calcwish.exe (
  set wish=%dir%cgnstools\calcwish.exe
  goto getscript
)
echo calcwish.exe not found
pause
goto done

:getscript
if exist %dir%cgnscalc.tcl (
  set script=%dir%cgnscalc.tcl
  goto run
)
if not "%CG_LIB_DIR%" == "" (
  if exist %CG_LIB_DIR%\cgnscalc.tcl (
    set script=%CG_LIB_DIR%\cgnscalc.tcl
    goto run
  )
)
if exist %dir%..\share\cgnscalc.tcl (
  if "%CG_LIB_DIR%" == "" set CG_LIB_DIR=%dir%..\share
  set script=%dir%..\share\cgnscalc.tcl
  goto run
)
echo cgnscalc.tcl not found
pause
goto done

:run
start /b %wish% %script% %1

:done
endlocal
