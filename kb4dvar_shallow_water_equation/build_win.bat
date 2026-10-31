@echo off
setlocal enabledelayedexpansion

set FC=gfortran
set SRC_DIR=seafield
set SOURCE=%SRC_DIR%\CalcSeaFieldMod_Alter.f90
set EXIT_CODE=0

echo ===== Fortran DLL Build Start =====

%FC% -O3 -fopenmp -shared -fPIC -J %SRC_DIR% -o %SRC_DIR%\CalcSeaFieldMod.dll %SOURCE% || set EXIT_CODE=1

del /q "%SRC_DIR%\*.o" 2>nul
del /q "%SRC_DIR%\*.obj" 2>nul
del /q "%SRC_DIR%\*.mod" 2>nul
del /q "*.mod" 2>nul

if !EXIT_CODE! == 0 (
    echo.
    echo All Fortran DLLs compiled successfully.
) else (
    echo.
    echo Compilation failed.
    exit /b 1
)

endlocal
