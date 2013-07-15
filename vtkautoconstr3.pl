#!/usr/bin/perl
use strict;

my $numdepth = $ARGV[0];
my $numdepthMax = $ARGV[1];
my $numpartition = $ARGV[2];

print "vtk auto-construction tool, version 0.3\n";
print "USAGE: perl vtkautoconstr3.pl min_depth max_depth #procs $!\n";
print "Creating vtk Files..\n";

my $i=0;
while (1){
	my $test_file = sprintf("r_vtk__%02d_%02d_%04d_%04d.vtk", $numdepth, $numdepthMax, $i, 0);	
	
	if (not -e $test_file) {
		if ($i == 0) {
			print "File not found: $test_file\n";
			exit;
		} else {
			print "Done ($i file(s)).\n";
			exit;
		}
	}
	
	my $output_file = sprintf(">out_%02d_%02d_%04d.vtk", $numdepth, $numdepthMax, $i);	
	open(OUT, $output_file) or die "Can't open $output_file for writing: $!\n";

	my $j=0;
	while ($j < $numpartition){
		my $input_file = sprintf("r_vtk__%02d_%02d_%04d_%04d.vtk", $numdepth, $numdepthMax, $i, $j);

		open(IN, $input_file) 
		or die "Can't open " .$input_file. ": $!\n";
	      	print OUT <IN>;
	      	close(IN);
		$j++;
	}	

	my $j=0;
	while ($j < $numpartition){
		my $input_file = sprintf("r_cell_%02d_%02d_%04d_%04d.vtk", $numdepth, $numdepthMax, $i, $j);
			
		open(IN, $input_file) 
		or last;
	      	print OUT <IN>;
	      	close(IN);
		$j++;
	}

	my $j=0;
	while ($j < $numpartition){
		my $input_file = sprintf("r_dof__%02d_%02d_%04d_%04d.vtk", $numdepth, $numdepthMax, $i, $j);
		
		open(IN, $input_file) 
		or last;
	      	print OUT <IN>;
	      	close(IN);
		$j++;
	}

	my $j=0;
	while ($j < $numpartition){
		my $input_file = sprintf("r_u____%02d_%02d_%04d_%04d.vtk", $numdepth, $numdepthMax, $i, $j);
		
		open(IN, $input_file) 
		or last;
	      	print OUT <IN>;
	      	close(IN);
		$j++;
	}

	my $j=0;
	while ($j < $numpartition){
		my $input_file = sprintf("r_v____%02d_%02d_%04d_%04d.vtk", $numdepth, $numdepthMax, $i, $j);
		
		open(IN, $input_file) 
		or last;
	      	print OUT <IN>;
	      	close(IN);
		$j++;
	}

	my $j=0;
	while ($j < $numpartition){
		my $input_file = sprintf("r_proc_%02d_%02d_%04d_%04d.vtk", $numdepth, $numdepthMax, $i, $j);
		
		open(IN, $input_file) 
		or last;
	      	print OUT <IN>;
	      	close(IN);
		$j++;
	}
	
	close(OUT) or die "Can't close $output_file after writing: $!\n";  	
	$i++;
}
print "Output VTK files created!\n";
	
 		
  
