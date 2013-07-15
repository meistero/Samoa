@echo off
cls

make -j -f WMakefile

if %errorlevel% neq 0 (
	echo "Aborting script: make returned an error."
	exit
)

wrun.bat %*
