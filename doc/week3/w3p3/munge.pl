#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO 'parse';
use Bio::Phylo::Util::CONSTANT ':objecttypes';

# process command line arguments
my $infile  = 'ungulates.nex';
my $outtree = 'ungulates.nwk';
my $data    = 'ungulates.tsv';
my $missing = 0.5;
GetOptions(
  'infile=s'  => \$infile,
  'outtree=s' => \$outtree,
  'data=s'    => \$data,
  'missing=f' => \$missing,
);

# read input file
my $proj = parse(
  '-format' => 'nexus',
  '-file'   => $infile,
  '-as_project' => 1,
);

# get data objects
my ($matrix) = grep { $_->get_type eq 'Continuous' } @{ $proj->get_matrices };
my ($tree)   = @{ $proj->get_items(_TREE_) };

# taxa in the tree but not in the matrix are going to be added as fully missing,
# taxa in the matrix but not in the tree are going to be pruned from the matrix
my @labels = @{ $matrix->get_charlabels };
my $nchar  = $matrix->get_nchar;
my @table;
for my $row ( @{ $matrix->get_entities } ) {
  my @char = $row->get_char;
  my $name = $row->get_name;
  if ( $tree->get_by_name($name) ) {
    @char = map { $_ eq '?' ? undef : $_ } @char;
    unshift @char, $name;
    push @table, \@char;
  }
}
for my $tip ( @{ $tree->get_terminals } ) {
  my $name = $tip->get_name;
  if ( not $matrix->get_by_name($name) ) {
    my @char = ($name);
    push(@char, undef) for 1 .. $nchar;
    push @table, \@char;
  }
}

# columns with more than $missing data are pruned
my @newlabels;
my @newmatrix = map { [ $_->[0] ] } @table;
for my $char ( 1 .. $nchar ) {
  my $hasdata = 0;
  for my $row ( 0 .. $#table ) {
    $hasdata++ if defined $table[$row]->[$char];
  }
  if ( ( $hasdata / scalar(@table) ) >= $missing ) {
    for my $row ( 0 .. $#table ) {
      push @{ $newmatrix[$row] }, $table[$row]->[$char];
    }
    my $label = $labels[ $char - 1 ];
    $label =~ s/'//g;
    push @newlabels, $label;
  }
}

# write table
{
  no warnings 'uninitialized';
  open my $outfh, '>', $data or die $!;
  print $outfh join("\t",undef,@newlabels), "\n";
  for my $row ( @newmatrix ) {
    print $outfh join("\t",@$row), "\n";
  }
}

# write tree
open my $outtreefh, '>', $outtree or die $!;
print $outtreefh $tree->to_newick;
