@echo off
echo =========================================
echo  GNSS-INS Project - Build ^& Run
echo =========================================

cd /d "%~dp0"

:: 自动探测并载入 Visual Studio 编译环境
set "VCVARS="
if exist "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" (
    set "VCVARS=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
)
if exist "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat" (
    set "VCVARS=C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat"
)
if exist "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvars64.bat" (
    set "VCVARS=C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvars64.bat"
)
if exist "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat" (
    set "VCVARS=C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat"
)

if "%VCVARS%"=="" (
    echo [WARNING] VS vcvars64.bat not found, trying cl.exe directly...
) else (
    call "%VCVARS%" >nul 2>&1
)

echo Compiling...
cl.exe /W4 /O2 /EHsc /std:c++14 .\main\main.cpp /Fegnss_ins.exe

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo [ERROR] Build failed! Please check the errors above.
    pause
    exit /b %ERRORLEVEL%
)

echo.
echo [SUCCESS] Build completed.
echo =========================================
echo.
.\gnss_ins.exe
pause
