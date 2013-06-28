@echo off
setlocal

rem the standard wish command will work for this
rem set wish=c:\progra~1\tcl\bin\wish83.exe

set dir=%~dps0
if exist %dir%cgconfig.bat (
  call %dir%cgconfig.bat
  goto getwish
)
if exist %dir%..\cgconfig.bat call %dir%..\cgconfig.bat

:getwish
if "%CG_BIN_DIR%" == "" set CG_BIN_DIR=%dir%

if not "%wish%" == "" goto getscript
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
if exist %dir%cgnsnodes.tcl (
  set script=%dir%cgnsnodes.tcl
  goto run
)
if not "%CG_LIB_DIR%" == "" (
  if exist %CG_LIB_DIR%\cgnsnodes.tcl (
    set script=%CG_LIB_DIR%\cgnsnodes.tcl
    goto run
  )
)
if exist %dir%..\share\cgnsnodes.tcl (
  if "%CG_LIB_DIR%" == "" set CG_LIB_DIR=%dir%..\share
  set script=%dir%..\share\cgnsnodes.tcl
  goto run
)
echo cgnsnodes.tcl not found
pause
goto done

:run
start /b %wish% %script% %1

:done
endlocal
