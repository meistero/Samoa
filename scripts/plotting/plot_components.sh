# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#!/bin/bash
echo "Plotting results..."

cd $1

rm -f darcy*.plt swe*.plt

for file in darcy*.log ; do
    flags=$(echo $file | grep -oE "(_no[a-zA-Z0-9]+)+")
    processes=$(echo $file | grep -oE "_p[0-9]+" | grep -oE "[0-9]+")
    threads=$(echo $file | grep -oE "_t[0-9]+" | grep -oE "[0-9]+")
    sections=$(echo $file | grep -oE "_s[0-9]+" | grep -oE "[0-9]+")
	layers=$(echo $file | grep -oE "_l[0-9]+" | grep -oE "[0-9]+")

    processes=${processes:-1}
    threads=${threads:-1}
    sections=${sections:-1}
    layers=${layers:-0}

    csplit $file "/Phase statistics:/" {*} &>/dev/null

    i=0
    for phase in xx* ; do
	    echo -n $(($processes * $threads * ($layers + 1))) \"$processes{/Symbol \\264}$threads{/Symbol \\264}$layers\"" " >> "darcy"$i".plt"
	    grep -E "r0.*Adaptions" $phase | grep -oE "(ET|time): [-]*[0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "darcy"$i".plt"
	    grep -E "r0.*Adaptions" $phase | grep -oE "integrity: [-]*[0-9]*\.[0-9]+" | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "darcy"$i".plt"
	    grep -E "r0.*Adaptions" $phase | grep -oE "load balancing: [-]*[0-9]*\.[0-9]+" | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "darcy"$i".plt"
	    grep -E "r0.*Adaptions" $phase | grep -oE "update neighbors: [-]*[0-9]*\.[0-9]+" | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "darcy"$i".plt"
	    grep -E "r0.*Transport" $phase | grep -oE "(ET|time): [-]*[0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "darcy"$i".plt"
	    grep -E "r0.*Gradient" $phase | grep -oE "(ET|time): [-]*[0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "darcy"$i".plt"
	    grep -E "r0.*Permeability" $phase | grep -oE "(ET|time): [-]*[0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "darcy"$i".plt"
	    grep -E "r0.*Pressure Solver" $phase | grep -oE "(ET|time): [-]*[0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "darcy"$i".plt"
	    grep -E "r0.*Error Estimate" $phase | grep -oE "(ET|time): [-]*[0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "darcy"$i".plt"
	    grep -E "r0.*Phase time" $phase | grep -oE "[-]*[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "darcy"$i".plt"
        grep -E "r0.*Element throughput" $phase | grep -oE "[-]*[0-9]+\.[0-9]+" | tr "\n" " " | cat >> "darcy"$i".plt"
        echo -n $layers >> "darcy"$i".plt"
	    echo "" >> "darcy"$i".plt"

        i=$(( $i + 1 ))
    done

    rm -f xx*
done

#Darcy: Cores Tics adap_time adap_ET integ_time lb_time updnb_time transp_time transp_ET grad_time grad_ET perm_time perm_ET pres_time pres_ET err_time err_ET phase_time ET    layers
#       1     2    3         4       5          6       7          8           9         10        11      12        13      14        15      16       17     18         19    20

for file in swe*.log ; do
	flags=$(echo $file | grep -oE "(_no[a-zA-Z0-9]+)+")
	processes=$(echo $file | grep -oE "_p[0-9]+" | grep -oE "[0-9]+")
	threads=$(echo $file | grep -oE "_t[0-9]+" | grep -oE "[0-9]+")
	sections=$(echo $file | grep -oE "_s[0-9]+" | grep -oE "[0-9]+")
	layers=$(echo $file | grep -oE "_l[0-9]+" | grep -oE "[0-9]+")

	processes=${processes:-1}
	threads=${threads:-1}
	sections=${sections:-1}
	layers=${layers:-0}

    csplit $file "/Phase statistics:/" {*} &>/dev/null

    i=0
    for phase in xx* ; do
	    echo -n $(($processes * $threads * ($layers + 1))) \"$processes{/Symbol \\264}$threads\"" " >> "swe"$i".plt"
	    grep -E "r0.*Adaptions" $phase | grep -oE "(ET|time): [-]*[0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "swe"$i".plt"
	    grep -E "r0.*Adaptions" $phase | grep -oE "integrity: [-]*[0-9]*\.[0-9]+" | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "swe"$i".plt"
	    grep -E "r0.*Adaptions" $phase | grep -oE "load balancing: [-]*[0-9]*\.[0-9]+" | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "swe"$i".plt"
	    grep -E "r0.*Adaptions" $phase | grep -oE "update neighbors: [-]*[0-9]*\.[0-9]+" | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "swe"$i".plt"
	    grep -E "r0.*Time steps" $phase | grep -oE "(ET|time): [-]*[0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "swe"$i".plt"
	    grep -E "r0.*Displace" $phase | grep -oE "(ET|time): [-]*[0-9]*\.?[0-9]+" | grep -oE "[-]*[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "swe"$i".plt"
	    grep -E "r0.*Phase time" $phase | grep -oE "[-]*[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "swe"$i".plt"
        grep -E "r0.*Element throughput" $phase | grep -oE "[-]*[0-9]+\.[0-9]+" | tr "\n" " " | cat >> "swe"$i".plt"
        echo -n $layers >> "swe"$i".plt"
	    echo ""  >> "swe"$i".plt"

        i=$(( $i + 1 ))
    done

    rm -f xx*
done

#SWE: Cores Tics adap_time adap_ET integ_time lb_time updnb_time tstep_time tstep_ET displ_time displ_ET phase_time ET  layers
#     1     2    3         4       5          6       7          8          9        10         11       12         13  14

i=0
for phase in darcy*.plt ; do
    sort -t" " -n -k 1,1 $phase -o $phase
    i=$(( $i + 1 ))
done

i=0
for phase in swe*.plt ; do
    sort -t" " -n -k 1,1 $phase -o $phase
    i=$(( $i + 1 ))
done

gnuplot &>/dev/null << EOT

set terminal postscript enhanced color font ',20'
set xlabel "Processes, Threads, Layers"
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
set style line 9 lc rgb "yellow"

set style data histogram
set style histogram rowstacked
set style fill solid border 0
set boxwidth 0.75

max(x, y) = (x > y ? x : y)

#*******
# Darcy
#*******

do for [i=1:20] {
    infile = sprintf('darcy%i.plt', i)

    set title "Component breakdown"
    set ylabel "Sec. per Mio. Elements per core"

    unset output

plot infile u (10 / \$15 / max(\$20, 1)) ls 7 t "Pressure Solver", \
	    '' u (1 / \$11 / max(\$20, 1)) ls 5 t "Gradient", \
	    '' u (1 / \$9 / max(\$20, 1)) ls 4 t "Transport", \
	    '' u (1 / \$13 / max(\$20, 1)) ls 6 t "Permeability", \
	    '' u (1 / \$17 / max(\$20, 1)) ls 9 t "Error Estimate", \
	    '' u (\$5 / \$3 * 1 / \$4 / max(\$20, 1)) ls 2 t "Conformity", \
        '' u ((\$3 - \$7) / \$3 * 1 / \$4 / max(\$20, 1)):xtic(2) ls 1 t "Adaption", \
	    '' u (\$7 / \$3 * 1 / \$4 / max(\$20, 1)) ls 8 t "Neighbor search", \
	    '' u (\$6 / \$3 * 1 / \$4 / max(\$20, 1)) ls 3 t "Load Balancing"

    outfile = sprintf('| ps2pdf - darcy%i_components.pdf', i)
    set output outfile

    replot
}

#*****
# SWE
#*****

do for [i=1:20] {
    infile = sprintf('swe%i.plt', i)

    set title "SWE component breakdown"
    set xlabel "Processes, Threads"
    set ylabel "Sec. per Mio. Elements per core"

    unset output

    plot infile u (1 / \$9) ls 4 t "Time step", \
        '' u (\$5 / \$3 * 1 / \$4) ls 2 t "Conformity", \
        '' u ((\$3 - \$7) / \$3 * 1 / \$4):xtic(2) ls 1 t "Adaption", \
        '' u (\$7 / \$3 * 1 / \$4) ls 8 t "Neighbor search", \
	    '' u (\$6 / \$3 * 1 / \$4) ls 3 t "Load Balancing" #, \
	    #'' u (1 / \$11) ls 5 t "Displace"


    outfile = sprintf('| ps2pdf - swe%i_components.pdf', i)
    set output outfile

    replot
}

EOT
