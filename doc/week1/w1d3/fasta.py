from Bio import SeqIO # sudo pip install biopython
with open("Artiodactyla.fas", "rU") as handle:
	
	# retain COI-5P records for each species
	species = {}
	for record in SeqIO.parse(handle, "fasta"):
		fields = record.description.split('|')
		if fields[2] == 'COI-5P':
			sp = fields[1]
			if sp in species:
				species[sp].append(record)
			else:
				species[sp] = [ record ]

	# write longest record for each species
	for sp in species:
		length = 0
		longest = None
		for record in species[sp]:
			if len(record.seq) > length:
				length = len(record.seq)
				longest = record
		print '>' + longest.description
		print longest.seq