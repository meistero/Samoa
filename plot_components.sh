# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#!/bin/bash
echo "Plotting results..."

cd $1
echo "#Component breakdown" > "darcy.plt"

for file in darcy*.log; do
	flags=$(echo $file | grep -oE "(_no[a-zA-Z0-9]+)+")
	processes=$(echo $file | grep -oE "_p[0-9]+" | grep -oE "[0-9]+")
	threads=$(echo $file | grep -oE "_t[0-9]+" | grep -oE "[0-9]+")
	sections=$(echo $file | grep -oE "_s[0-9]+" | grep -oE "[0-9]+")

	processes=${processes:-1}
	threads=${threads:-1}
	sections=${sections:-1}

	echo -n $(($processes * $threads)) \"$processes-$threads\"" " >> "darcy.plt"
	grep -E "r0.*Adaptions" $file | grep -oE "(travs|time): [0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	grep -E "r0.*Adaptions" $file | grep -oE "integrity: [0-9]*\.[0-9]+" | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	grep -E "r0.*Adaptions" $file | grep -oE "load balancing: [0-9]*\.[0-9]+" | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	grep -E "r0.*Adaptions" $file | grep -oE "update neighbors: [0-9]*\.[0-9]+" | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	grep -E "r0.*Transport" $file | grep -oE "(travs|time): [0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	grep -E "r0.*Gradient" $file | grep -oE "(travs|time): [0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	grep -E "r0.*Permeability" $file | grep -oE "(travs|time): [0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	grep -E "r0.*Pressure Solver" $file | grep -oE "(travs|time): [0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	grep -E "r0.*Time step phase time" $file | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	grep -E "r0.*Cells[ ]+:[ ]+[0-9]+" $file | grep -oE "Cells[ ]+:[ ]+[0-9]+" | grep -oE "[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	echo "" >> "darcy.plt"
done

echo "#Component breakdown" > "swe.plt"

for file in swe*.log; do
	flags=$(echo $file | grep -oE "(_no[a-zA-Z0-9]+)+")
	processes=$(echo $file | grep -oE "_p[0-9]+" | grep -oE "[0-9]+")
	threads=$(echo $file | grep -oE "_t[0-9]+" | grep -oE "[0-9]+")
	sections=$(echo $file | grep -oE "_s[0-9]+" | grep -oE "[0-9]+")
	
	processes=${processes:-1}
	threads=${threads:-1}
	sections=${sections:-1}

	echo -n $(($processes * $threads)) \"$processes-$threads\"" " >> "swe.plt"
	grep -E "r0.*Adaptions" $file | grep -oE "(travs|time): [0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	grep -E "r0.*Adaptions" $file | grep -oE "integrity: [0-9]*\.[0-9]+" | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	grep -E "r0.*Adaptions" $file | grep -oE "load balancing: [0-9]*\.[0-9]+" | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	grep -E "r0.*Adaptions" $file | grep -oE "update neighbors: [0-9]*\.[0-9]+" | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	grep -E "r0.*Time steps" $file | grep -oE "(travs|time): [0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	grep -E "r0.*Displace" $file | grep -oE "(travs|time): [0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	grep -E "r0.*Time Step phase time" $file | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	grep -E "r0.*Cells[ ]+:[ ]+[0-9]+" $file | grep -oE "Cells[ ]+:[ ]+[0-9]+" | grep -oE "[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	echo ""  >> "swe.plt"
done

sort -t" " -n -k 1,1 darcy.plt -o darcy.plt
sort -t" " -n -k 1,1 swe.plt -o swe.plt

if [ -n "$2" ]; then
    echo "#Component breakdown filtered using sed -n $2" > darcy.plt.tmp
    echo "#Component breakdown filtered using sed -n $2" > swe.plt.tmp
    
    sed -n $2 darcy.plt >> darcy.plt.tmp 
    sed -n $2 swe.plt >> swe.plt.tmp 

    mv darcy.plt.tmp darcy.plt
    mv swe.plt.tmp swe.plt
fi

sort -t" " -n -k 1,1 darcy.plt -o darcy.plt
sort -t" " -n -k 1,1 swe.plt -o swe.plt

#gnuplot &> /dev/null << EOT
gnuplot << EOT

set terminal postscript enhanced color font ',20'
set xlabel "Cores"
set key below font ",20" spacing 1.0 width -2
set xtics rotate
set yrange [0:*]
set auto x

set style line 999 lt 2 lw 8 lc rgb "black"

set for [n=1:64] style line n lt 1 lw 2
set style line 1 lc rgb "cyan"
set style line 2 lc rgb "orange"
set style line 3 lc rgb "magenta"
set style line 4 lc rgb "red"
set style line 5 lc rgb "blue"
set style line 6 lc rgb "green"
set style line 7 lc rgb "brown"
set style line 8 lc rgb "purple"

set style data histogram
set style histogram rowstacked
set style fill solid border 0
set boxwidth 0.75

#*******
# Darcy
#*******

set title "Darcy component breakdown"
set ylabel "Sec. per core (wall clock time)"
set output '| ps2pdf - darcy_components.pdf'
	
plot "darcy.plt" u (\$8) ls 2 t "Conformity", \
    '' u (\$6 - \$12):xtic(2) ls 1 t "Adaption", \
	'' u (\$12) ls 8 t "Neighbor search", \
	'' u (\$10) ls 3 t "Load Balancing", \
	'' u (\$16) ls 5 t "Gradient", \
	'' u (\$14) ls 4 t "Transport", \
	'' u (\$18) ls 6 t "Permeability", \
	'' u (\$22) ls 7 t "Pressure Solver"

set title "Darcy component breakdown - normalized"
set ylabel "Sec. per element (CPU time)"
set output '| ps2pdf - darcy_components_norm.pdf'
	
plot "darcy.plt" u (10.0 * \$22/\$21 * \$1/\$25) ls 7 t "Pressure Solver", \
	'' u (\$16/\$15 * \$1/\$25) ls 5 t "Gradient", \
	'' u (\$14/\$13 * \$1/\$25) ls 4 t "Transport", \
	'' u (\$18/\$15 * \$1/\$25) ls 6 t "Permeability", \
	'' u (\$8/\$5 * \$1/\$25) ls 2 t "Conformity", \
    '' u ((\$6 - \$12)/\$5 * \$1/\$25):xtic(2) ls 1 t "Adaption", \
	'' u (\$12/\$5 * \$1/\$25) ls 8 t "Neighbor search", \
	'' u (\$10/\$5 * \$1/\$25) ls 3 t "Load Balancing"

#*****
# SWE
#*****

set title "SWE component breakdown"
set ylabel "Sec. per core (wall clock time)"
set output '| ps2pdf - swe_components.pdf'

plot    "swe.plt" u (\$14) ls 4 t "Time step", \
	    '' u (\$16) ls 5 t "Displace", \
	    '' u (\$8) ls 2 t "Conformity", \
        '' u (\$6 - \$12):xtic(2) ls 1 t "Adaption", \
	    '' u (\$12) ls 8 t "Neighbor search", \
	    '' u (\$10) ls 3 t "Load Balancing"

set title "SWE component breakdown - normalized"
set ylabel "Sec. per element (CPU time)"
set output '| ps2pdf - swe_components_norm.pdf'
	
plot "swe.plt" u (\$14/\$13 * \$1/\$19) ls 4 t "Time step", \
    '' u (\$8/\$5 * \$1/\$19) ls 2 t "Conformity", \
    '' u ((\$6 - \$12)/\$5 * \$1/\$19):xtic(2) ls 1 t "Adaption", \
    '' u (\$12/\$5 * \$1/\$19) ls 8 t "Neighbor search", \
	'' u (\$10/\$5 * \$1/\$19) ls 3 t "Load Balancing"
#	'' u (\$16/\$15 * \$1/\$19) ls 5 t "Displace", \


EOT
