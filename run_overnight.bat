@echo off
REM ========================================
REM BTEH Overnight Pipeline - Windows HPC
REM Runs H2O and SSDM training for all species
REM ========================================

echo.
echo ╔════════════════════════════════════════════════════════════════╗
echo ║  BTEH Pipeline - Overnight Run                                 ║
echo ╚════════════════════════════════════════════════════════════════╝
echo.
echo Start time: %date% %time%
echo.

REM Navigate to scripts directory
cd /d "%~dp0BTEH\scripts_RemoteSensing for E and C_"
if %errorlevel% neq 0 (
    echo ERROR: Could not find scripts directory
    echo Current directory: %cd%
    pause
    exit /b 1
)

echo Current directory: %cd%
echo.

REM ========================================
REM STAGE 1: H2O Training
REM ========================================
echo ════════════════════════════════════════════════════════════════
echo   STAGE 1: H2O Model Training
echo ════════════════════════════════════════════════════════════════
echo.

echo [1/4] H2O Training - Run A (3 species)...
echo Start: %time%
Rscript 03_h2o_train.R --run A --mode REPRO
if %errorlevel% neq 0 (
    echo.
    echo ✗ ERROR in H2O Run A
    echo Check logs in: ../../logs/03_h2o_train_A.log
    pause
    exit /b %errorlevel%
)
echo ✓ Completed: %time%
echo.

echo [2/4] H2O Training - Run B (6 species)...
echo Start: %time%
Rscript 03_h2o_train.R --run B --mode REPRO
if %errorlevel% neq 0 (
    echo.
    echo ✗ ERROR in H2O Run B
    echo Check logs in: ../../logs/03_h2o_train_B.log
    pause
    exit /b %errorlevel%
)
echo ✓ Completed: %time%
echo.

REM ========================================
REM STAGE 2: SSDM Training
REM ========================================
echo ════════════════════════════════════════════════════════════════
echo   STAGE 2: SSDM Model Training
echo ════════════════════════════════════════════════════════════════
echo.

echo [3/4] SSDM Training - Run A (3 species)...
echo Start: %time%
Rscript 04_ssdm_train.R --run A --mode REPRO
if %errorlevel% neq 0 (
    echo.
    echo ✗ ERROR in SSDM Run A
    echo Check logs in: ../../logs/04_ssdm_train_A.log
    pause
    exit /b %errorlevel%
)
echo ✓ Completed: %time%
echo.

echo [4/4] SSDM Training - Run B (6 species)...
echo Start: %time%
Rscript 04_ssdm_train.R --run B --mode REPRO
if %errorlevel% neq 0 (
    echo.
    echo ✗ ERROR in SSDM Run B
    echo Check logs in: ../../logs/04_ssdm_train_B.log
    pause
    exit /b %errorlevel%
)
echo ✓ Completed: %time%
echo.

REM ========================================
REM SUCCESS!
REM ========================================
echo.
echo ╔════════════════════════════════════════════════════════════════╗
echo ║  Pipeline Complete!                                            ║
echo ╚════════════════════════════════════════════════════════════════╝
echo.
echo End time: %date% %time%
echo.
echo Results saved in:
echo   - H2O models:  ../../results/H2O/
echo   - SSDM models: ../../results/SSDM/
echo.
echo Next steps:
echo   1. Check results folders
echo   2. Run comparison scripts (05_*.R)
echo   3. Generate figures (07_appendix.R)
echo.
pause
