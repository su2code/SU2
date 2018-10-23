for /f %%i in ('dir /ad /b') do devenv %%i\%%i.sln /Build Release /project %%i
