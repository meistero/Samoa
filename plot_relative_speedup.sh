# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#!/bin/bash
echo "Plotting results..."

cd $1

echo "#Element throughput by sections and threads" > "darcy.plt"
echo "#Threads ET per (initialization, time steps) per sections" >> "darcy.plt"

for file in darcy*.log; do
	flags=$(echo $file | grep -oE "(_no[a-zA-Z0-9]+)+")
	processes=$(echo $file | grep -oE "_p[0-9]+" | grep -oE "[0-9]+")
	threads=$(echo $file | grep -oE "_t[0-9]+" | grep -oE "[0-9]+")
	sections=$(echo $file | grep -oE "_s[0-9]+" | grep -oE "[0-9]+")

	processes=${processes:-1}
	threads=${threads:-1}
	sections=${sections:-1}

	echo -n $processes $threads $sections" " >> "darcy.plt"
	grep -E "r0.*Element throughput" $file | grep -oE "[0-9]+\.[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	echo -n $flags >> "darcy.plt"
	echo "" >> "darcy.plt"
done


echo "#Cell update throughput and flux solver throughput by sections and threads" > "swe.plt"
echo "#Threads CUT FST per sections" >> "swe.plt"

for file in swe*.log; do
	flags=$(echo $file | grep -oE "(_no[a-zA-Z0-9]+)+")
	processes=$(echo $file | grep -oE "_p[0-9]+" | grep -oE "[0-9]+")
	threads=$(echo $file | grep -oE "_t[0-9]+" | grep -oE "[0-9]+")
	sections=$(echo $file | grep -oE "_s[0-9]+" | grep -oE "[0-9]+")
	
	processes=${processes:-1}
	threads=${threads:-1}
	sections=${sections:-1}

	echo -n $processes $threads $sections" "  >> "swe.plt"
	    
	grep -E "r0.*Cell update throughput" $file | grep -oE "[0-9]+\.[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	grep -E "r0.*Flux solver throughput" $file | grep -oE "[0-9]+\.[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	echo -n $flags >> "swe.plt"
	echo ""  >> "swe.plt"
done

sort -t" " -n -k 3,3 -k 1,1 -k 2,2 darcy.plt -o darcy.plt
sort -t" " -n -k 3,3 -k 1,1 -k 2,2 swe.plt -o swe.plt

gnuplot &> /dev/null << EOT

set terminal postscript enhanced color font ',30'
set xlabel "Cores"
set ylabel "Mio. Elements per sec. per core"
unset key
set log x 2

title(n) = sprintf("%d section(s)", n)

set style line 999 lt 2 lw 8 lc rgb "black"

set for [n=1:64] style line n lt 1 lw 8
set style line 1 lc rgb "cyan"
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
set yrange [0:*]

plot for [n=1:64] "darcy.plt" u (\$1*\$2):(\$3 == n ? \$4 / (\$1*\$2) : 1/0) ls n w linespoints t title(n)

set output '| ps2pdf - darcy_elem_init_rel.pdf'

replot GPVAL_DATA_Y_MAX w lines ls 999 title "reference"

#********

set title "Darcy element throughput - time steps"
unset output
set log x 2
set yrange [0:*]

plot for [n=1:64] "darcy.plt" u (\$1*\$2):(\$3 == n ? \$5/ (\$1*\$2) : 1/0) ls n w linespoints t title(n)
		
set output '| ps2pdf - darcy_elem_rel.pdf'

replot GPVAL_DATA_Y_MAX w lines ls 999 title "reference"

#*****
# SWE
#*****

set title "SWE flux solver throughput"
unset output
set log x 2
set yrange [0:*]

plot for [n=1:64] "swe.plt" u (\$1*\$2):(\$3 == n? \$5 / (\$1*\$2) : 1/0) ls n w linespoints t title(n)

set output '| ps2pdf - swe_flux_rel.pdf'

replot GPVAL_DATA_Y_MAX w lines ls 999 title "reference"


#********

set title "SWE cell update throughput"
unset output
set log x 2
set yrange [0:*]

plot for [n=1:64] "swe.plt" u (\$1*\$2):(\$3 == n? \$4 / (\$1*\$2) : 1/0) ls n w linespoints t title(n)

set output '| ps2pdf - swe_cells_rel.pdf'

replot GPVAL_DATA_Y_MAX w lines ls 999 title "reference"

EOT
