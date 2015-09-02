#compile the code with different compilation flags
for mobility in brooks-corey linear quadratic
do
    scons config=gnu_release.py flux_solver=upwind asagi=noasagi mobility=$mobility -j4
    mv "bin/samoa_darcy_notasks_upwind_gnu" "bin/samoa_darcy_notasks_upwind_gnu_"$mobility    
done

#create directories (if they don't exist yet) for the output files and copy the reference data to subfolders
for scenario in density_rarefaction density_shock inflow_rarefaction inflow_shock
do
    mkdir -p "output/"$scenario

    for mobility in linear quadratic brooks-corey
    do
        mkdir -p "output/"$scenario"/"$mobility"/ref"

        cp "scripts/python/"$scenario"_"$mobility""* "output/"$scenario"/"$mobility"/ref"

        for dmax in {2..20..2}
        do
            mkdir -p "output/"$scenario"/"$mobility"/d"$dmax
        done
    done
done

#run the scenarios for the given parameter space
for dmax in {2..20..2}
do
    for mobility in brooks-corey linear quadratic
    do
        for scenario in density_rarefaction density_shock inflow_rarefaction inflow_shock
        do
            if [ $scenario = "inflow_rarefaction" ]
            then
                #run inflow rarefaction
                command="bin/samoa_darcy_notasks_upwind_gnu_"$mobility" -sections 1 -threads 4 -dmax "$dmax" -xmloutput -tout 1.0e4 -tmax 1.0e4 -epsilon 1.0e-12 -S_ref_th 0.2 -p_ref_th 1.0e20 -p_prod 1.0 -inflow 0.5434396505 -g_x 0.0 -courant 0.1"
            fi

            if [ $scenario = "inflow_shock" ]
            then
                #run inflow shock
                command="bin/samoa_darcy_notasks_upwind_gnu_"$mobility" -sections 1 -threads 4 -dmax "$dmax" -xmloutput -tout 1.0e4 -tmax 1.0e4 -epsilon 1.0e-12 -S_ref_th 0.2 -p_ref_th 1.0e20 -p_prod 1.0 -inflow -0.5434396505 -g_x 0.0 -courant 0.1"
            fi

            if [ $scenario = "density_rarefaction" ]
            then
                #run density rarefaction
                command="bin/samoa_darcy_notasks_upwind_gnu_"$mobility" -sections 1 -threads 4 -dmax "$dmax" -xmloutput -tout 1.0e4 -tmax 1.0e4 -epsilon 1.0e-12 -S_ref_th 0.2 -p_ref_th 1.0e20 -p_prod 1.0 -inflow 0.0 -g_x 9.80665 -courant 0.1"
            fi

            if [ $scenario = "density_shock" ]
            then
                #run density shock
                command="bin/samoa_darcy_notasks_upwind_gnu_"$mobility" -sections 1 -threads 4 -dmax "$dmax" -xmloutput -tout 1.0e4 -tmax 1.0e4 -epsilon 1.0e-12 -S_ref_th 0.2 -p_ref_th 1.0e20 -p_prod 1.0 -inflow 0.0 -g_x -9.80665 -courant 0.1"
            fi
            
            echo $command > "output/"$scenario"_"$mobility"_d"$dmax".log"
            $command >> "output/"$scenario"_"$mobility"_d"$dmax".log"

            mv output/*.*vtu output/*.log "output/"$scenario"/"$mobility"/d"$dmax
            rm output/*.csv
        done
    done
done
