# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#!/bin/bash
echo "Plotting results..."

cd $1

echo "#Element throughput by sections and threads" > "darcy.plt"
echo "#Threads ET per (initialization, time steps) per sections" >> "darcy.plt"

for file in darcy*.log; do
	processes=$(echo $file | grep -oE "_p[0-9]+" | grep -oE "[0-9]+")
	threads=$(echo $file | grep -oE "_t[0-9]+" | grep -oE "[0-9]+")
	sections=$(echo $file | grep -oE "_s[0-9]+" | grep -oE "[0-9]+")

	processes=${processes:-1}
	threads=${threads:-1}
	sections=${sections:-1}

	echo -n $processes $threads $sections" " >> "darcy.plt"
	grep -E "r0.*Element throughput" $file | grep -oE "[0-9]+\.[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	echo "" >> "darcy.plt"
done


echo "#Cell update throughput and flux solver throughput by sections and threads" > "swe.plt"
echo "#Threads CUT FST per sections" >> "swe.plt"

for file in swe*.log; do
	processes=$(echo $file | grep -oE "_p[0-9]+" | grep -oE "[0-9]+")
	threads=$(echo $file | grep -oE "_t[0-9]+" | grep -oE "[0-9]+")
	sections=$(echo $file | grep -oE "_s[0-9]+" | grep -oE "[0-9]+")
	
	processes=${processes:-1}
	threads=${threads:-1}
	sections=${sections:-1}

	echo -n $processes $threads $sections" "  >> "swe.plt"
	    
	grep -E "r0.*Cell update throughput" $file | grep -oE "[0-9]+\.[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	grep -E "r0.*Flux solver throughput" $file | grep -oE "[0-9]+\.[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	echo ""  >> "swe.plt"
done

sort -t" " -n -k 3,3 -k 1,1 -k 2,2 darcy.plt -o darcy.plt
sort -t" " -n -k 3,3 -k 1,1 -k 2,2 swe.plt -o swe.plt

gnuplot &> /dev/null << EOT

set terminal postscript enhanced color font ',30'
set xlabel "Threads"
set ylabel "Mio. Elements per sec."
set key left top
set logscale xy 2

title(n) = sprintf("%d section(s)", n)

set style line 1 lt 2 lw 8 lc rgb "black"

set for [n=2:64] style line n lt 1 lw 8
set style line 2 lc rgb "orange"
set style line 4 lc rgb "magenta"
set style line 8 lc rgb "red"
set style line 16 lc rgb "blue"
set style line 32 lc rgb "green"

#*******
# Darcy
#*******

set title "Darcy element throughput - initialization"
unset output
set xrange [*:*]
set yrange [*:*]

plot for [n=1:64] "darcy.plt" u (\$1*\$2):(\$3 == n ? \$4 : 1/0) ls n w linespoints t title(n)

set output '| ps2pdf - darcy_elem_init_log.pdf'

replot log(GPVAL_DATA_Y_MIN) / log(GPVAL_DATA_X_MIN) * x ls 1 w lines title "reference"

#********

set title "Darcy element throughput - time steps"
unset output
set xrange [*:*]
set yrange [*:*]

plot for [n=1:64] "darcy.plt" u (\$1*\$2):(\$3 == n ? \$5 : 1/0) ls n w linespoints t title(n)
		
set output '| ps2pdf - darcy_elem_log.pdf'

replot log(GPVAL_DATA_Y_MIN) / log(GPVAL_DATA_X_MIN) * x ls 1 w lines title "reference"

#*****
# SWE
#*****

set title "SWE flux solver throughput"
unset output
set xrange [*:*]
set yrange [*:*]

plot for [n=1:64] "swe.plt" u (\$1*\$2):(\$3 == n ? \$5 : 1/0) ls n w linespoints t title(n)

set output '| ps2pdf - swe_flux_log.pdf'

replot log(GPVAL_DATA_Y_MIN) / log(GPVAL_DATA_X_MIN) * x ls 1 w lines title "reference"


#********

set title "SWE cell update throughput"
unset output
set xrange [*:*]
set yrange [*:*]

plot for [n=1:64] "swe.plt" u (\$1*\$2):(\$3 == n ? \$4 : 1/0) ls n w linespoints t title(n)

set output '| ps2pdf - swe_cells_log.pdf'

replot log(GPVAL_DATA_Y_MIN) / log(GPVAL_DATA_X_MIN) * x ls 1 w lines title "reference"


EOT
