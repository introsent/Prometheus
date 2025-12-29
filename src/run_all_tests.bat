@echo off
echo ===================================
echo Area Light Sampling Test Suite
echo ===================================

echo.
echo Generate ground truth...
Prometheus.exe --gt

echo.
echo [1/4] Testing Uniform Sampling...
Prometheus.exe -t -s uniform

echo.
echo [2/4] Testing Area Importance Sampling...
Prometheus.exe -t -s areaimportance

echo.
echo [3/4] Testing Hierarchical Flux Sampling...
Prometheus.exe -t -s hierarchical

echo.
echo [4/4] Testing Visibility Aware Sampling...
Prometheus.exe -t -s visibilityaware


echo.
echo All tests complete!
echo Check the generated folders for results.
pause