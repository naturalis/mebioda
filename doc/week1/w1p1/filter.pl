my %map;
while(<>) {	
	my @fields = split /\t/;
	
	# the 4th column, index 3, holds the length
	next if $fields[3] < 150;
	
	# the first column has semicolon-separated fields
	my ( $otu, $size ) = split /;/, $fields[0];
	
	# we have not yet seen this OTU id
	if ( not $map{$otu} ) {
		$map{$otu} = [ @fields[1..3] ]; # store it
	}
	else {
	
		# we've already seen this OTU, but this match is higher
		if ( $fields[2] > $map{$otu}->[2] ) {
			warn $fields[2], '>', $map{$otu}->[2];
			$map{$otu} = [ @fields[1..3] ]; # replace it
		}
	}
}

# print output
for my $otu ( sort { $a cmp $b } keys %map ) {
	my $row = $map{$otu};
	print $otu, "\t", join("\t", @{ $map{$otu} } ), "\n";
}