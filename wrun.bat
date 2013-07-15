@echo off
cls

main.exe %*

if %errorlevel% neq 0 (
	echo "Aborting script: program returned an error."
	exit
)
