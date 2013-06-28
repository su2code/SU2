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
if not "%wish%" == "" goto getscript
if exist %dir%calcwish.exe (
  set wish=%dir%calcwish.exe
  goto getscript
)
if exist %dir%cgnstools\calcwish.exe (
  set wish=%dir%cgnstools\calcwish.exe
  goto getscript
)
if exist %dir%cgiowish.exe (
  set wish=%dir%cgiowish.exe
  goto getscript
)
if exist %dir%cgnstools\cgiowish.exe (
  set wish=%dir%cgnstools\cgiowish.exe
  goto getscript
)
echo calcwish.exe or cgiowish.exe not found
pause
goto done

:getscript
if exist %dir%unitconv.tcl (
  set script=%dir%unitconv.tcl
  goto run
)
if not "%CG_LIB_DIR%" == "" (
  if exist %CG_LIB_DIR%\unitconv.tcl (
    set script=%CG_LIB_DIR%\unitconv.tcl
    goto run
  )
)
if exist %dir%..\share\unitconv.tcl (
  set CG_LIB_DIR=%dir%..\share
  set script=%dir%..\share\unitconv.tcl
  goto run
)
echo unitconv.tcl not found
pause
goto done

:run
start /b %wish% %script% %1

:done
endlocal
