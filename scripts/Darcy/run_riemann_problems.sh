#compile the code with different compilation flags

for mobility in brooks-corey linear quadratic
do
    exe="samoa_darcy_gnu_"$mobility
    scons config=gnu_release.py flux_solver=upwind asagi=noasagi mobility=$mobility exe=$exe -j4  
done

#create directories (if they don't exist yet) for the output files and copy the reference data to subfolders
for scenario in density_rarefaction density_shock inflow_rarefaction inflow_shock
do
    mkdir -p "output/"$scenario

    for mobility in linear quadratic brooks-corey
    do
        mkdir -p "output/"$scenario"/"$mobility"/ref"

        cp "scripts/python/"$scenario"_"$mobility""* "output/"$scenario"/"$mobility"/ref"
    done
done

#run the scenarios for the given parameter space
for dmax in {2..20..1}
do
    for mobility in brooks-corey linear quadratic
    do
        for scenario in density_rarefaction density_shock inflow_rarefaction inflow_shock
        do
            exe="samoa_darcy_gnu_"$mobility
            output_dir="output/"$scenario"/"$mobility"/d"$dmax
            mkdir -p $output_dir

            if [ $scenario = "inflow_rarefaction" ]
            then
                #run inflow rarefaction
                command="bin/$exe -sections 1 -threads 4 -dmax "$dmax" -xmloutput -tout 1.0e4 -tmax 1.0e4 -epsilon 1.0e-12 -S_ref_th 0.2 -p_ref_th 1.0e20 -p_prod 1.0 -inflow 0.5434396505 -g_x 0.0 -courant 1.0 -output_dir $output_dir"
            fi

            if [ $scenario = "inflow_shock" ]
            then
                #run inflow shock
                command="bin/$exe -sections 1 -threads 4 -dmax "$dmax" -xmloutput -tout 1.0e4 -tmax 1.0e4 -epsilon 1.0e-12 -S_ref_th 0.2 -p_ref_th 1.0e20 -p_prod 1.0 -inflow -0.5434396505 -g_x 0.0 -courant 1.0 -output_dir $output_dir"
            fi

            if [ $scenario = "density_rarefaction" ]
            then
                #run density rarefaction
                command="bin/$exe -sections 1 -threads 4 -dmax "$dmax" -xmloutput -tout 1.0e4 -tmax 1.0e4 -epsilon 1.0e-12 -S_ref_th 0.2 -p_ref_th 1.0e20 -p_prod 1.0 -inflow 0.0 -g_x 9.80665 -courant 1.0 -output_dir $output_dir"
            fi

            if [ $scenario = "density_shock" ]
            then
                #run density shock
                command="bin/$exe -sections 1 -threads 4 -dmax "$dmax" -xmloutput -tout 1.0e4 -tmax 1.0e4 -epsilon 1.0e-12 -S_ref_th 0.2 -p_ref_th 1.0e20 -p_prod 1.0 -inflow 0.0 -g_x -9.80665 -courant 1.0 -output_dir $output_dir"
            fi
            
            echo $command > $output_dir"/"$scenario"_"$mobility"_d"$dmax".log"
            $command | tee -a $output_dir"/"$scenario"_"$mobility"_d"$dmax".log"

            rm $output_dir/*.csv
        done
    done
done
