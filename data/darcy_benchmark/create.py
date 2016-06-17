#!/usr/bin/python

import argparse
import math
import numpy
from pnoise import noise2D, noise3D

def create_grid(width = 128, height = 128, roughness = 0.2, nc_file = None, vtk_file = None):
	#init arrays
	nc_x = []
	nc_y = []
	nc_values = {}

	### Try to import netcdf ##
	if nc_file:
		try:
			import netCDF4
			from netCDF4 import Dataset

			nc = Dataset(nc_file, 'w')
			nc.createDimension('x', width)
			nc.createDimension('y', height)
			nc_x = nc.createVariable('x', 'f4', ('x'))
			nc_y = nc.createVariable('y', 'f4', ('y'))
			nc_values = nc.createVariable('z', 'f4', ('x', 'y'))

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

	### Create a vtk file ##	
	if vtk_file:
		vtk = open(vtk_file, 'w')
		vtk.write('# vtk DataFile Version 2.0\n')
		vtk.write('Perlin noise boundary\n')
		vtk.write('ASCII\n')
		vtk.write('DATASET STRUCTURED_POINTS\n')
		vtk.write('DIMENSIONS %i %i %i\n' %(width, height, 1))
		vtk.write('ORIGIN %g %g %g\n' %(-0.5, -0.5, 0.0))
		vtk.write('SPACING %g %g %g\n' %(1.0, 1.0, 1.0))	
		vtk.write('\n')
		vtk.write('POINT_DATA %i\n' %(width * height))
		vtk.write('\n')
		vtk.write('SCALARS noise float 1\n')
		vtk.write('LOOKUP_TABLE default\n')

	#set x and y spacings
	nc_x[:] = numpy.linspace(0.5 / float(height), (float(width) - 0.5) / float(height), width)
	nc_y[:] = numpy.linspace(0.5 / float(height) , 1.0 - 0.5 / float(height), height)

	#compute depth
	depth = int(math.log(height))

	for j in range(0, height):
		y = nc_y[j]

		for i in range(0, width):
			x = nc_x[i]
	
			out = -7.0e-8 * (noise2D(10.0 * x - 2.0, 10.0 * y, depth, roughness) + 0.7 - 4.0 * y * (1.0 - y)) + 0.5e-8
			out = max(0.0, min(1.0e-8, out))
	
			nc_values[j,i] = out

	if vtk_file:
		for j in range(0, height):
			for i in range(0, width):
				vtk.write("%g\n" %nc_values[i,j])
		vtk.close()

	if nc_file:
		nc.close()

def main():
	#parse command line arguments
	
	parser = argparse.ArgumentParser(description='Creates a random noise field.')
	parser.add_argument('-width', type=int, default=128, help='width of the resulting grid (default: 128)')
	parser.add_argument('-height', type=int, default=128, help='height of the resulting grid (default: 128)')
	parser.add_argument('-roughness', type=float, default=0.2, help='surface roughness, 0.0 = smooth, 1.0 = white noise (default: 0.2)')

	parser.add_argument('-nc_file', help='output netcdf file name (default: none)')
	parser.add_argument('-vtk_file', help='output vtk file name (default: none)')

	args = parser.parse_args()

	create_grid(args.width, args.height, args.roughness, args.nc_file, args.vtk_file)

main()
