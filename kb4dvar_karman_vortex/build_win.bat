@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0"

set FC=gfortran
set SRC_DIR=ns_model
set FFLAGS=-O3 -fPIC -cpp -fno-range-check -I %SRC_DIR% -J %SRC_DIR%
set EXIT_CODE=0

echo ===== Fortran DLL Build Start =====

%FC% -O3 -fPIC -cpp -fno-range-check -c %SRC_DIR%\mod_variables.f90 -J %SRC_DIR% -o %SRC_DIR%\mod_variables.o || set EXIT_CODE=1

if !EXIT_CODE! neq 0 goto :end
%FC% %FFLAGS% -c %SRC_DIR%\sub_utils.f90 -o %SRC_DIR%\sub_utils.o || set EXIT_CODE=1
if !EXIT_CODE! neq 0 goto :end
%FC% %FFLAGS% -c %SRC_DIR%\sub_initial.f90 -o %SRC_DIR%\sub_initial.o || set EXIT_CODE=1
if !EXIT_CODE! neq 0 goto :end
%FC% %FFLAGS% -c %SRC_DIR%\sub_bc_outer.f90 -o %SRC_DIR%\sub_bc_outer.o || set EXIT_CODE=1
if !EXIT_CODE! neq 0 goto :end
%FC% %FFLAGS% -c %SRC_DIR%\sub_bc_wall.f90 -o %SRC_DIR%\sub_bc_wall.o || set EXIT_CODE=1
if !EXIT_CODE! neq 0 goto :end
%FC% %FFLAGS% -c %SRC_DIR%\sub_rhs3rd_fwd.f90 -o %SRC_DIR%\sub_rhs3rd_fwd.o || set EXIT_CODE=1
if !EXIT_CODE! neq 0 goto :end
%FC% %FFLAGS% -c %SRC_DIR%\sub_hsmac_fwd.f90 -o %SRC_DIR%\sub_hsmac_fwd.o || set EXIT_CODE=1
if !EXIT_CODE! neq 0 goto :end
%FC% %FFLAGS% -c %SRC_DIR%\calc_forward_ns.f90 -o %SRC_DIR%\calc_forward_ns.o || set EXIT_CODE=1
if !EXIT_CODE! neq 0 goto :end

%FC% -fopenmp -shared -o %SRC_DIR%\calc_forward_ns.dll ^
  %SRC_DIR%\mod_variables.o ^
  %SRC_DIR%\sub_utils.o ^
  %SRC_DIR%\sub_initial.o ^
  %SRC_DIR%\sub_bc_outer.o ^
  %SRC_DIR%\sub_bc_wall.o ^
  %SRC_DIR%\sub_rhs3rd_fwd.o ^
  %SRC_DIR%\sub_hsmac_fwd.o ^
  %SRC_DIR%\calc_forward_ns.o || set EXIT_CODE=1

:end
del /q "%SRC_DIR%\*.o" 2>nul
del /q "%SRC_DIR%\*.mod" 2>nul

if !EXIT_CODE! == 0 (
    echo.
    echo All Fortran DLLs compiled successfully.
) else (
    echo.
    echo Compilation failed.
    exit /b 1
)

endlocal
