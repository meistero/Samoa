#!/home/meistero/Documents/ParaView-4.3.1-Linux-64bit/bin/pvpython

#### import the simple module from the paraview
from paraview.simple import *
from paraview import vtk
import os
import sys
import glob
import argparse
import re
from math import log 

def main():
    parser = argparse.ArgumentParser(description='Plot error of Riemann problems in a csv file')
    parser.add_argument('-d', '--data_dir', default='output/', help='data directory')
    args = parser.parse_args()

    producers = xrange(5)
    
    #redirect output to stderr (we are not really interested in info messages from Paraview)
    default_stdout = sys.stdout    
    sys.stdout = sys.stderr
    
    pressure_dict = {}
    prod_dict = {}

    for result_dir in os.listdir(args.data_dir):
        if os.path.isdir(os.path.join(args.data_dir, result_dir)):
            (spe_tag, compiler, layers, averaging, depth) = result_dir.split("_")
            layers = int(layers[1:])
            depth = int(depth[1:])                

            pvtu = glob.glob(os.path.join(args.data_dir, result_dir, '*_0.pvtu'))

            if pvtu:
                darcy = XMLPartitionedUnstructuredGridReader(FileName = pvtu)
                darcy.CellArrayStatus = []
                darcy.PointArrayStatus = ['pressure']
                
                probe = ProbeLocation(Input=darcy, ProbeType='Fixed Radius Point Source')

                if layers > 0:
                    probe.ProbeType.Center = [0.0, 0.0, 51.81]
                else:
                    probe.ProbeType.Center = [0.0, 0.0, 0.0]
                
                probe.ProbeType.Radius = 0.0
                probe.Tolerance = 2.0e-16

                rawData = servermanager.Fetch(probe)

                pressure_dict[averaging, layers, depth] = rawData.GetPointData().GetArray('pressure').GetValue(0)

                Delete(darcy)
                Delete(probe)

            csv = glob.glob(os.path.join(args.data_dir, result_dir, '*_0.csv'))

            if csv:
                darcy_data = CSVReader(FileName = csv[0])
                rawData = servermanager.Fetch(darcy_data)
                points = TableToPoints(Input=darcy_data)

                rawData = servermanager.Fetch(points)

                prod_dict[averaging, layers, depth] = [rawData.GetPointData().GetArray(' oil rate').GetValue(i) for i in producers]

                Delete(darcy)
                Delete(points)
    
    #redirect output back to stdout    
    sys.stdout = default_stdout
    
    print "Pressure by layers, averaging, depth:"

    sys.stdout.write(" l aver ")
    
    for depth in xrange(30):
        if any([k[2] == depth for k in pressure_dict]):
            sys.stdout.write("%8d " % (depth))
    
    print ""

    for layers in xrange(86):
        if not any([k[1] == layers for k in pressure_dict]):
            continue

        for averaging in ['arithmetic', 'geometric', 'harmonic']:
            sys.stdout.write("%2d %s " % (layers, averaging[0:4]))

            for depth in xrange(30):
                try:
                    pressure = pressure_dict[averaging, layers, depth]
                    sys.stdout.write("%8.2f " % (pressure))
                except:
                    pass
            
            print ""

    for producer in producers:
        print "\nProduction rate "+ str(producer) + " by layers, averaging, depth:"

        sys.stdout.write(" l aver ")
        
        for depth in xrange(30):
            if any([k[2] == depth for k in prod_dict]):
                sys.stdout.write("%8d " % (depth))
        
        print ""

        for layers in xrange(86):
            if not any([k[1] == layers for k in prod_dict]):
                continue

            for averaging in ['arithmetic', 'geometric', 'harmonic']:
                sys.stdout.write("%2d %s " % (layers, averaging[0:4]))

                for depth in xrange(30):
                    try:
                        rates = prod_dict[averaging, layers, depth]
                        sys.stdout.write("%8.2f " % (rates[producer]))
                    except:
                        pass
                
                print ""
main()
