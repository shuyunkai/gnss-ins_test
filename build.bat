@echo off
echo =========================================
echo Compiling GNSS-INS Project...
echo =========================================

cd /d "%~dp0"

cl.exe /W4 /O2 /EHsc /std:c++14 .\main\main.cpp /Fegnss_ins.exe

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo [ERROR] Build failed! Please check the errors above.
    exit /b %ERRORLEVEL%
)

echo.
echo [SUCCESS] Build completed successfully.
echo Running gnss_ins.exe ...
echo =========================================
.\gnss_ins.exe

