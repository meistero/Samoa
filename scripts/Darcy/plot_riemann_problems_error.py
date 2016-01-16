#!/home/meistero/Documents/ParaView-4.3.1-Linux-64bit/bin/pvpython

#### import the simple module from the paraview
from paraview.simple import *
from paraview import vtk
import os
import glob
import argparse
import re
from math import log 

def natural_sorted(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def main():
    parser = argparse.ArgumentParser(description='Plot error of Riemann problems in a csv file')
    
    parser.add_argument('-w', '--width', type=int, default=2048, help='width of the test plane')
    parser.add_argument('-t', '--height', type=int, default=64, help='height of the test plane')
    parser.add_argument('-o', '--once', type=bool, nargs='?', const=True, default=False, help='if true, the program is executed only on one data set')
    parser.add_argument('-d', '--data_dir', default='output/', help='data directory')
    
    args = parser.parse_args()

    print "scenario,model,dmax,L1 error,cell count,dof count,execution time, convergence rate (dx_min), convergence rate (DoFs)"

    for scenario_dir in natural_sorted(os.listdir(args.data_dir)):
        if os.path.isdir(os.path.join(args.data_dir, scenario_dir)):
            for model_dir in natural_sorted(os.listdir(os.path.join(args.data_dir, scenario_dir))):
                if os.path.isdir(os.path.join(args.data_dir, scenario_dir, model_dir)):
                    # load the second time step from the reference solution
                    ref_dir = 'ref'
                    ref = CSVReader(FileName = natural_sorted(glob.glob(os.path.join(args.data_dir, scenario_dir, model_dir, ref_dir, '*_1.csv'))))

                    # create a new 'Plane'
                    plane = Plane()
                    plane.Origin = [0.0, 0.0, 0.0]
                    plane.Point1 = [1.0, 0.0, 0.0]
                    plane.Point2 = [0.0, 1.0, 0.0]
                    plane.XResolution = args.width
                    plane.YResolution = args.height
            
                    # get active source
                    programmableFilter = ProgrammableFilter(ref)
                    programmableFilter.OutputDataSetType = 'vtkPolyData'

                    # set script
                    programmableFilter.Script = "    pdi = self.GetInput()\n    pdo =  self.GetPolyDataOutput()\n\n    numPoints = pdi.GetNumberOfRows()\n\n    coords = vtk.vtkPoints()\n    saturation = vtk.vtkFloatArray()\n    saturation.SetName('reference')\n    saturation.SetNumberOfComponents(1)\n    saturation.SetNumberOfTuples(2*numPoints)\n\n    x_array = pdi.GetRowData().GetArray('X')\n    s_array = pdi.GetRowData().GetArray(' saturation')\n\n    for i in xrange(numPoints):\n        coords.InsertPoint(2*i, x_array.GetValue(i), 0.0, 0.0)\n        coords.InsertPoint(2*i+1, x_array.GetValue(i), 1.0, 0.0)\n        saturation.SetValue(2*i, s_array.GetValue(i))\n        saturation.SetValue(2*i+1, s_array.GetValue(i))\n\n    pdo.SetPoints(coords)\n    pdo.GetPointData().AddArray(saturation)\n\n    pdo.Allocate(numPoints-1, 1)\n\n    for i in xrange(numPoints - 1):\n        pdo.InsertNextCell(9, 4, [2*i, 2*(i+1), 2*(i+1)+1, 2*i+1])"

                    resampled_ref = ResampleWithDataset(Input=programmableFilter, Source=plane)

                    sim_results = []

                    for sim_dir in natural_sorted(os.listdir(os.path.join(args.data_dir, scenario_dir, model_dir))):
                        if os.path.isdir(os.path.join(args.data_dir, scenario_dir, model_dir, sim_dir)) and sim_dir.find('d') >= 0:
                            if not glob.glob(os.path.join(args.data_dir, scenario_dir, model_dir, sim_dir, '*_1.pvtu')):
                                continue

                            # load the second time step from the simulation output
                            darcy = XMLPartitionedUnstructuredGridReader(FileName = natural_sorted(glob.glob(os.path.join(args.data_dir, scenario_dir, model_dir, sim_dir, '*_1.pvtu'))))

                            darcy.CellArrayStatus = []
                            darcy.PointArrayStatus = ['saturation']

                            # resample both the simulation and the reference data on the plane to compare them
                            resampled_sim = ResampleWithDataset(Input=darcy, Source=plane)

                            appendAttributes = AppendAttributes(Input=[resampled_sim, resampled_ref])
                            
                            # compute L1 norm: absolute contributions
                            calc_l1 = Calculator(Input=appendAttributes)
                            calc_l1.ResultArrayName = 'error_l1'
                            calc_l1.Function = 'abs(saturation - reference)'

                            integrateVariables = IntegrateVariables(Input=calc_l1)

                            phase_times = []
                            cell_counts = []
                            dof_counts = []
                            
                            f = open(glob.glob(os.path.join(args.data_dir, scenario_dir, model_dir, sim_dir, '*.log'))[0])
                            for line in f:
                                if re.match("^.*Phase time[  :]+[0-9]*\.[0-9]+ s$", line):
                                    phase_times += re.findall("[  :]+([0-9]*\.[0-9]+)", line)
                                
                                if re.match("^.*Cells[  :]+[0-9]+$", line):
                                    cell_counts += re.findall("[  :]+([0-9]+)", line)
                                
                                if re.match("^.*Inner nodes[  :]+[0-9]+$", line):
                                    dof_counts += re.findall("[  :]+([0-9]+)", line)

                                if re.match("^.*Boundary nodes.*[  :]+[0-9]+[  ]+[0-9]+$", line):
                                    b_nodes = re.findall("[  :]+([0-9]+)[  ]+([0-9]+)", line)
                                    dof_counts[-1] = repr(int(dof_counts[-1]) + int(b_nodes[0][0]) + int(b_nodes[0][1]))

                            # obtain error
                            rawData = servermanager.Fetch(integrateVariables)
                            sim_result = (int(sim_dir[1:]), rawData.GetPointData().GetArray('error_l1').GetValue(0), int(cell_counts[1]), int(dof_counts[1]), float(phase_times[1]))
                            #print scenario, model, depth, l1 error, cell count, DoF count, simulation time, convergence rate (dx_min), convergence rate (DoFs)"

                            sim_results.append(sim_result)
                            coarse_result = sim_results[0]
                            
                            try:
                                alpha = (log(coarse_result[1]) - log(sim_result[1])) / ((sim_result[0] - coarse_result[0])*log(2.0)/2.0)
                            except ZeroDivisionError:
                                alpha = 0.0                                

                            try:
                                beta = (log(coarse_result[1]) - log(sim_result[1])) / (log(sim_result[3]) - log(coarse_result[3]))
                            except ZeroDivisionError:
                                beta = 0.0

                            print "%s,%s,%i,%.4g,%i,%i,%.4g,%.4g,%.4g" % (scenario_dir, model_dir, sim_result[0], sim_result[1], sim_result[2], sim_result[3], sim_result[4],  alpha, beta)
                   
                    if (args.once):
                        break

                    # delete all objects in the pipeline
                    for f in GetSources().values():
                        Delete(f)

        if (args.once):
            break
main()
