#!/usr/bin/python

import argparse
import math
import numpy
from pnoise import noise2D, noise3D

def create_grid(width = 128, height = 128, depth = 8, roughness = 0.2, nc_file = None):
    #init arrays
    nc_x = []
    nc_y = []
    nc_Kx = {}
    nc_Ky = {}
    nc_Kz = {}
    nc_Phi = {}

    ### Try to import netcdf ##
    if nc_file:
        try:
            import netCDF4
            from netCDF4 import Dataset

            nc = Dataset(nc_file, 'w')
            nc.createDimension('x', width)
            nc.createDimension('y', height)
            nc.createDimension('z', depth)
            nc_x = nc.createVariable('x', 'f4', ('x'))
            nc_y = nc.createVariable('y', 'f4', ('y'))
            nc_z = nc.createVariable('z', 'f4', ('z'))
            nc_Kx = nc.createVariable('Kx', 'f4', ('z', 'y', 'x'))
            nc_Ky = nc.createVariable('Ky', 'f4', ('z', 'y', 'x'))
            nc_Kz = nc.createVariable('Kz', 'f4', ('z', 'y', 'x'))
            nc_Phi = nc.createVariable('Phi', 'f4', ('z', 'y', 'x'))

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

    fwidth =  1.0
    fheight = 1.0
    fdepth =  1.0

    #set x and y spacings
    nc_x[:] = numpy.linspace(-0.5 * fwidth + 0.5 * fwidth / float(width), 0.5 * fwidth - 0.5 * fwidth / float(width), width)
    nc_y[:] = numpy.linspace(-0.5 * fheight + 0.5 * fheight / float(height), 0.5 * fheight - 0.5 * fheight / float(height), height)
    nc_z[:] = numpy.linspace(-0.5 * fdepth + 0.5 * fdepth / float(depth), 0.5 * fdepth - 0.5 * fdepth / float(depth), depth)

    #compute number of noise levels
    noise_levels = int(math.log(max(width, height, depth)))

    for k in xrange(depth):
        z = nc_z[k]

        for j in xrange(height):
            y = nc_y[j]

            for i in xrange(width):
                x = nc_x[i]

                out = 0.5 * noise3D(10.0 * x - 2.0, 10.0 * y, 10.0 * z, noise_levels, roughness) + 0.25 * (4.0 / (fheight * fheight) * (0.5 * fheight - y) * (0.5 * fheight + y) + 4.0 / (fdepth * fdepth) * (0.5 * fdepth - z) * (0.5 * fdepth + z))
                out = max(0.0, out)

                nc_Kx[k,j,i] = 10.0 ** (7.0 * out - 3.0)
                nc_Ky[k,j,i] = 10.0 ** (7.0 * out - 3.0)
                nc_Kz[k,j,i] = 10.0 ** (7.0 * out - 5.0)
                nc_Phi[k,j,i] = 0.2

    nc.close()

def main():
    #parse command line arguments
    
    parser = argparse.ArgumentParser(description='Creates a random noise field.')
    parser.add_argument('-x', '--nx', type=int, default=128, help='width of the resulting grid (default: 128)')
    parser.add_argument('-y', '--ny', type=int, default=128, help='height of the resulting grid (default: 128)')
    parser.add_argument('-z', '--nz', type=int, default=1, help='height of the resulting grid (default: 8)')
    parser.add_argument('-r', '--roughness', type=float, default=0.2, help='surface roughness, 0.0 = smooth, 1.0 = white noise (default: 0.2)')

    parser.add_argument('-n', '--nc_file', help='output netcdf file name (default: none)')

    args = parser.parse_args()

    create_grid(args.nx, args.ny, args.nz, args.roughness, args.nc_file)

main()
