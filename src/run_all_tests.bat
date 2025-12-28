@echo off
echo ===================================
echo Area Light Sampling Test Suite
echo ===================================

echo.
echo [1/3] Testing Uniform Sampling...
Prometheus.exe -t -s uniform

echo.
echo [2/3] Testing Area Importance Sampling...
Prometheus.exe -t -s areaimportance

echo.
echo [3/3] Testing Hierarchical Flux Sampling...
Prometheus.exe -t -s hierarchical

echo.
echo All tests complete!
echo Check the generated folders for results.
pause