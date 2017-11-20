my %seen;
while(<>){
	my @fields = split /\t/, $_;
	print join( "\t", @fields ) if not $seen{$fields[0]}++; 
}

