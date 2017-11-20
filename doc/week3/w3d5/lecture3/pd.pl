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
my $replicates = 1000; # number of replicates
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

make_distributions( $outfile, $tree, $replicates, @sizes );

# write two-column table with resampled distances
sub make_distributions {
	my ( $outfile, $tree, $replicates, @sizes ) = @_;
	open my $fh, '>', $outfile or die $!;	
	print $fh join( ',', map { "size$_" } @sizes ), "\n";	
	for my $rep ( 1 .. $replicates ) {
		$log->debug("Computing replicate $rep");
		my @res;
		for my $size ( @sizes ) {
			push @res, resample( 
				parse_tree(
					'-format' => 'newick',
					'-file'   => $treefile,
				), 
				$size 
			);
		}
		print $fh join( ',', @res ), "\n";		
	}
}

# take a single sample of tree length
sub resample {
	my ( $tree, $size ) = @_;
	my @tips = shuffle( @{ $tree->get_terminals } );
	my @sample = @tips[ 0 .. $size - 1 ];
	my $subtree = $tree->keep_tips(\@sample);
	return $subtree->calc_tree_length;
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