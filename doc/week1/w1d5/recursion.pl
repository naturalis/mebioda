#!/usr/bin/perl
use Bio::Phylo::IO 'parse_tree';

parse_tree(
	-format => 'nexus',
	-file   => 'tree.nex',
)->visit_depth_first(
	-pre    => sub { print 'PRE: ', shift->get_name, "\n" },
	-post   => sub { print 'POST: ', shift->get_name, "\n" }
);