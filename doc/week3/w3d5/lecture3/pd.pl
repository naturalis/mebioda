#!/usr/bin/perl
use strict;
use warnings;
use Text::CSV;
use Getopt::Long;
use List::Util qw'shuffle sum';
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my $mapfile; # maps checklist columns to tree tips
my $dbfile; # location of the checklist-cleaned.csv
my $treefile; # location of the tree file
my $replicates = 100_000; # number of replicates
my @sizes = ( 2 .. 5 ); # subtree sizes (i.e. number of tips)
my $verbosity = WARN; # log level
my $outfile; # CSV file with resampled distances
GetOptions(
	'mapfile=s'    => \$mapfile,
	'dbfile=s'     => \$dbfile,
	'treefile=s'   => \$treefile,
	'replicates=s' => \$replicates,
	'sizes=s'      => \@sizes,
	'verbose+'     => \$verbosity,
	'outfile=s'    => \$outfile,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main'
);
my $csv = Text::CSV->new({
	'binary' => 1,
});
$log->info("Going to read newick tree from --treefile='$treefile'");
my $tree = parse_tree(
	'-format' => 'newick',
	'-file'   => $treefile
);
$log->info("Done reading newick tree");

my $mat = distmat( $tree );
make_distributions( $outfile, $mat, $replicates, @sizes );

# write two-column table with resampled distances
sub make_distributions {
	my ( $outfile, $mat, $replicates, @sizes ) = @_;
	open my $fh, '>', $outfile or die $!;	
	print $fh join( ',', map { "size$_" } @sizes ), "\n";	
	for my $rep ( 1 .. $replicates ) {
		my @res;
		for my $size ( @sizes ) {
			push @res, resample( $mat, $size );
		}
		print $fh join( ',', @res ), "\n";
	}
}

# compute square distance matrix
sub distmat {
	my $tree = shift;
	my @tips = @{ $tree->get_terminals };
	my %dist = map { $_->get_name => {} } @tips;
	for my $i ( 0 .. $#tips - 1 ) {
		for my $j ( $i + 1 .. $#tips ) {
			my $d = $tips[$i]->calc_patristic_distance( $tips[$j] );
			my $i_n = $tips[$i]->get_name;
			my $j_n = $tips[$j]->get_name;
			$dist{ $i_n }->{ $j_n } = $d;
			$dist{ $j_n }->{ $i_n } = $d;
			$log->debug("Distance $i_n => $j_n = $d");
		}
	}
	return \%dist;
}

# take a single sample of average distance between $size tips
sub resample {
	my ( $mat, $size ) = @_;
	my @tips = shuffle( keys %$mat );
	my @dist;
	for my $i ( 0 .. $size - 2 ) {
		my $i_n = $tips[$i];
		for my $j ( $i + 1 .. $size - 1 ) {
			my $j_n = $tips[$j];
			my $d = $mat->{$i_n}->{$j_n};
			push @dist, $d;
			$log->debug("Distance $i_n => $j_n = $d");
		}
	}
	return sum(@dist)/scalar(@dist);
}

# read mapfile
my %map;
if ( $mapfile ) {
	$log->info("Going to read mapping from --mapfile='$mapfile'");
	open my $fh, '<', $mapfile or die $!;
	while( my $row = $csv->getline($fh) ) {
		$map{ $row->[0] } = $row->[1];
	}
}

# iterate over instances 
if ( $dbfile ) {
	open my $fh, '<', $dbfile or die $!;
	my @header;
	while( my $row = $csv->getline($fh) ) {
		if ( not @header ) {
			@header = @$row;
		}
		else {
		
		}
	}	
}