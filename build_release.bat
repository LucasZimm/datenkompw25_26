@echo off 
REM ================================================ 
REM Build script (auto reconfigure for new cpp files) 
REM ================================================ 
set PROJECT_ROOT=%~dp0 
set BUILD_DIR=%~dp0/build
set BIN_DIR=%~dp0/bin

echo Project Root: %PROJECT_ROOT%

echo 🔧 Deleting old build

if exist "%BUILD_DIR%" ( 
    rmdir /s /q "%BUILD_DIR%" 
)
if exist "%BIN_DIR%" ( 
    rmdir /s /q "%BIN_DIR%" 
)

echo 🔄 Configuring CMake...
echo PROJECT_ROOT: %PROJECT_ROOT%
cmake -S %~dp0 -B build

if not exist "%BUILD_DIR%" ( 
    mkdir "%BUILD_DIR%" 
) 
echo 🔄 Reconfiguring CMake... 
cmake -S %~dp0 -B build -DCMAKE_BUILD_TYPE=Release 


echo 🔨 Building... 
cmake --build "%BUILD_DIR%" --config Release --clean-first 

echo. 
echo ✅ Build completed successfully! 
pause