import sys
from Bio import SeqIO # sudo pip install biopython
with open(sys.argv[1], "rU") as handle:
    
    # Example: relabel sequences as Genus_species-ID
    for record in SeqIO.parse(handle, "fasta"):
        fields = record.description.split('|')
        name = fields[1].replace(' ', '_')
        print '>' + name + '-' + fields[0]
        print record.seq