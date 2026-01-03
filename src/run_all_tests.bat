@echo off
echo ===================================
echo Area Light Sampling Test Suite
echo Running all tests for all scenes
echo ===================================

:: Define scenes to test
set SCENES=flux bunny quad occlusion

for %%S in (%SCENES%) do (
    echo.
    echo ===================================
    echo Testing Scene: %%S
    echo ===================================

    :: Create scene directory if it doesn't exist
    if not exist "tests\%%S" mkdir "tests\%%S"

    echo Generating ground truth for %%S...
    Prometheus.exe --gt --scene %%S

    if errorlevel 1 (
        echo ERROR: Failed to generate ground truth for %%S
        pause
        exit /b 1
    )

    echo.
    echo [1/4] Testing Uniform Sampling...
    Prometheus.exe -t -s uniform --scene %%S

    echo.
    echo [2/4] Testing Area Importance Sampling...
    Prometheus.exe -t -s areaimportance --scene %%S

    echo.
    echo [3/4] Testing Hierarchical Flux Sampling...
    Prometheus.exe -t -s hierarchical --scene %%S

    echo.
    echo [4/4] Testing Visibility Aware Sampling...
    Prometheus.exe -t -s visibilityaware --scene %%S

    echo.
    echo Completed tests for scene: %%S
    echo Results saved to: tests\%%S\
)

echo.
echo ===================================
echo All tests complete!
echo.
echo Results organized by scene:
echo   tests\flux\
echo   tests\bunny\
echo   tests\quad\
echo   tests\occlusion\
echo ===================================
pause