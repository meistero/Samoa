# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#!/bin/bash
echo "Plotting ASAGI results..."

cd $1

echo "#Element throughput by sections and threads" > "darcy.plt"
echo "#Threads ET per (initialization, time steps) per sections" >> "darcy.plt"

for file in darcy*.log; do
	asagimode=$(echo $file | grep -oE "_a[0-9]+" | grep -oE "[0-9]+")
	processes=$(echo $file | grep -oE "_p[0-9]+" | grep -oE "[0-9]+")
	sections=$(echo $file | grep -oE "_s[0-9]+" | grep -oE "[0-9]+")

	asagimode=${asagimode:0}
	processes=${processes:0}
	sections=${sections:-1}

	echo -n $asagimode $sections $processes" " >> "darcy.plt"

	grep -E "r0.*Element throughput" $file | grep -oE "[0-9]+\.[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	echo "" >> "darcy.plt"
done


echo "#Element update throughput and flux solver throughput by sections and threads" > "swe.plt"
echo "#Threads ET per sections" >> "swe.plt"

for file in swe*.log; do
	asagimode=$(echo $file | grep -oE "_a[0-9]+" | grep -oE "[0-9]+")
	processes=$(echo $file | grep -oE "_p[0-9]+" | grep -oE "[0-9]+")
	sections=$(echo $file | grep -oE "_s[0-9]+" | grep -oE "[0-9]+")
	
	asagimode=${asagimode:0}
	processes=${processes:0}
	sections=${sections:-1}

	echo -n $asagimode $sections $processes" "  >> "swe.plt"
	    
	grep -E "r0.*Element throughput" $file | grep -oE "[0-9]+\.[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	echo ""  >> "swe.plt"
done

sort -t" " -n -k 2,2 -k 3,3 -k 1,1 darcy.plt -o darcy.plt
sort -t" " -n -k 2,2 -k 3,3 -k 1,1 -k 2,2 swe.plt -o swe.plt

gnuplot &> /dev/null << EOT

set terminal postscript enhanced color font ',20'
set xlabel "ASAGI mode"
set ylabel "M/s"
set xtics ("default" 0, "pass through" 1, "no mpi" 2, "no mpi, small cache" 3, "large grid" 4)

set style line 1 lt 2 lw 8 lc rgb "black"

set for [n=2:64] style line n lt 1 lw 8
set style line 2 lc rgb "orange"
set style line 4 lc rgb "magenta"
set style line 8 lc rgb "red"
set style line 16 lc rgb "blue"
set style line 32 lc rgb "green"

init_title(n) = sprintf("init %d section(s)", n)
ts_title(n) = sprintf("ts %d section(s)", n)

set title "Darcy - element throughput"
set output '| ps2pdf - darcy_asagi.pdf'

plot for [n=1:64] "darcy.plt" u 1:(\$2 == n ? \$4 : 1/0) ls 8 w linespoints t init_title(n), \
	for [n=1:64] "darcy.plt" u 1:(\$2 == n ? \$5 : 1/0) ls 16 w linespoints t ts_title(n)

set title "SWE - element throughput"
set output '| ps2pdf - swe_asagi.pdf'

plot for [n=1:64] "swe.plt" u 1:(\$2 == n ? \$4 : 1/0) ls 8 w linespoints t init_title(n), \
	for [n=1:64] "swe.plt" u 1:(\$2 == n ? \$5 : 1/0)  ls 16 w linespoints t ts_title(n)

EOT
