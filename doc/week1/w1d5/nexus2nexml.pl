#!/usr/bin/perl
use Bio::Phylo::IO qw'parse';

print parse(
	-format     => 'nexus',
	-handle     => \*DATA,
	-as_project => 1,
)->to_xml; 

__DATA__
#NEXUS
begin taxa;
	dimensions ntax=5;
	taxlabels
		A
		B
		C
		D
		E
	;		
end;
begin trees;
	translate
		1 A,
		2 B,
		3 C,
		4 D,
		5 E;
	tree t1 = (((1,2),(3,4)),5);
end;