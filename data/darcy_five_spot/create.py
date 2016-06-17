#!/usr/bin/python

import argparse
import math
import numpy

def create_grid(width, height, depth, perm_file, por_file, nc_file, scaling):
    #init arrays
    nc_x = []
    nc_y = []
    nc_z = []
    nc_phi = {}
    nc_Kx = {}
    nc_Ky = {}
    nc_Kz = {}

    ### Try to import netcdf ##
    try:
        import netCDF4
        from netCDF4 import Dataset

        nc = Dataset(nc_file, 'w')
        nc.createDimension('x', scaling*width)
        nc.createDimension('y', scaling*height)
        nc.createDimension('z', scaling*depth)
        nc_x = nc.createVariable('x', 'f4', ('x'))
        nc_y = nc.createVariable('y', 'f4', ('y'))
        nc_z = nc.createVariable('z', 'f4', ('z'))

    except ImportError:
        # Install numpy-devel files
        ## Download http://code.google.com/p/netcdf4-python/downloads/detail?name=netCDF4-0.9.3.tar.gz
        # Extract
        # cd netCDF4-0.9.3
        # python setup-nc3.py build
        # sudo python setup-nc3.py install

        import warnings
        warnings.warn("NetCDF not found, file is not created")
        nc_file = False

    fwidth =  1200.0 * 12.0 * 0.0254    # 1200 ft = 1200 ft * 12 in/ft * 0.0254 m/in = 365.76 m
    fheight = 2200.0 * 12.0 * 0.0254    # 2200 ft = 1200 ft * 12 in/ft * 0.0254 m/in = 670.56 m
    fdepth =   170.0 * 12.0 * 0.0254    #  170 ft =  170 ft * 12 in/ft * 0.0254 m/in = 51.816 m

    #set x and y spacings
    nc_x[:] = numpy.linspace(-0.5 * fwidth + 0.5 * fwidth / float(scaling*width), 0.5 * fwidth - 0.5 * fwidth / float(scaling*width), scaling*width)
    nc_y[:] = numpy.linspace(-0.5 * fheight + 0.5 * fheight / float(scaling*height), 0.5 * fheight - 0.5 * fheight / float(scaling*height), scaling*height)
    nc_z[:] = numpy.linspace(0.5 * fdepth / float(scaling*depth), fdepth - 0.5 * fdepth / float(scaling*depth), scaling*depth)

    por_table_2D = numpy.loadtxt(perm_file)
    perm_table_2D = numpy.loadtxt(por_file)
    por_table = [item for sublist in por_table_2D for item in sublist]
    perm_table = [item for sublist in perm_table_2D for item in sublist]

    nc_Kx = nc.createVariable('Kx', 'f4', ('z', 'y', 'x'))
    nc_Ky = nc.createVariable('Ky', 'f4', ('z', 'y', 'x'))
    nc_Kz = nc.createVariable('Kz', 'f4', ('z', 'y', 'x'))
    nc_Phi = nc.createVariable('Phi', 'f4', ('z', 'y', 'x'))

    n = 0

    for k in xrange(0, depth):
        for j in xrange(0, height):
            for i in xrange(0, width):
                nc_Kx[ scaling * (depth - 1 - k) : scaling * (depth - k), \
                       scaling * (height - 1 - j) : scaling * (height - j), \
                       scaling * i : scaling * (i + 1)] = perm_table[n]

                n = n + 1
        print "Kx Slice: %i of %i" % (k + 1, depth)

    for k in xrange(0, depth):
        for j in xrange(0, height):
            for i in xrange(0, width):
                nc_Ky[ scaling * (depth - 1 - k) : scaling * (depth - k), \
                       scaling * (height - 1 - j) : scaling * (height - j), \
                       scaling * i : scaling * (i + 1)] = perm_table[n]

                n = n + 1
        print "Ky Slice: %i of %i" % (k + 1, depth)

    for k in xrange(0, depth):
        for j in xrange(0, height):
            for i in xrange(0, width):
                nc_Kz[ scaling * (depth - 1 - k) : scaling * (depth - k), \
                       scaling * (height - 1 - j) : scaling * (height - j), \
                       scaling * i : scaling * (i + 1)] = perm_table[n]

                n = n + 1
        print "Kz Slice: %i of %i" % (k + 1, depth)

    assert(n == len(perm_table))

    n = 0
    
    for k in xrange(0, depth):
        for j in xrange(0, height):
            for i in xrange(0, width):
                nc_Phi[ scaling * (depth - 1 - k) : scaling * (depth - k), \
                        scaling * (height - 1 - j) : scaling * (height - j), \
                        scaling * i : scaling * (i + 1)] = por_table[n]
                n = n + 1
        print "Phi Slice: %i of %i" % (k + 1, depth)

    assert(n == len(por_table))

    nc.close()

def main():
    #parse command line arguments
    
    parser = argparse.ArgumentParser(description='Creates a netcdf file from a dat file.')
    parser.add_argument('-k', '--perm_file', default=None, help='source permeability file name (default: none)')
    parser.add_argument('-p', '--por_file', default=None, help='source porosity file name (default: none)')
    parser.add_argument('-x', '--nx', type=int, default=60, help='width of the source grid (default: 60)')
    parser.add_argument('-y', '--ny', type=int, default=220, help='height of the source grid (default: 220)')
    parser.add_argument('-z', '--nz', type=int, default=85, help='depth of the source grid (default: 85)')

    parser.add_argument('-c', '--scaling', type=int, default=1, help='scaling of output grid (default: 1)')

    parser.add_argument('-n', '--nc_file', default=None, help='output netcdf file name (default: none)')

    args = parser.parse_args()
    
    create_grid(args.nx, args.ny, args.nz, args.perm_file, args.por_file, args.nc_file, args.scaling)

main()
