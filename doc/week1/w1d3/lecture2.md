Barcoding workflow and MSA data
===============================

Barcoding
---------

![](barcode_pipeline.jpg)

Barcode Of Life Data Systems ([BOLDSYSTEMS](http://www.boldsystems.org/))
-------------------------------------------------------------------------

- Stores records about [specimens](http://www.boldsystems.org/index.php/Public_RecordView?processid=ABMC137-05)
- Includes marker sequence [data](fasta.fas), images, lat/lon coordinates, etc.
- Can query [taxonomically](http://www.boldsystems.org/index.php/Public_SearchTerms?query=Artiodactyla[tax])
  and download [all sequences](Artiodactyla.fas)
- Identification services:
  - COI for animals
  - ITS for fungi
  - rbcL and matK for plants
  
Fetching sequences through the [URL API](http://www.boldsystems.org/index.php/api_home)
---------------------------------------------------------------------------------------

Data from URLs can be downloaded on the command line using [curl](https://curl.haxx.se/):

```bash
# Fetch all sequences for Artiodactyla as FASTA
$ curl -o Artiodactyla.fas http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Artiodactyla
```

The BOLD sequence data service API returns a [FASTA file](Artiodactyla.fas), which holds 
multiple sequences, unaligned. The 
[definition line](https://en.wikipedia.org/wiki/FASTA_format#Description_line) is 
formatted as:

```
>ID|Scientific binomial|marker|...
``` 

With this we can use standard UNIX command line tools 
[grep, cut, sort and uniq](http://www.tldp.org/LDP/abs/html/textproc.html) to do some 
basic checks, e.g.:

```bash
$ grep '>' Artiodactyla.fas | cut -f 3 -d '|' | sort | uniq -c
  30 16S
2316 COI-5P
 184 COI-5P
 188 COII
 214 COII
 188 COXIII
 214 COXIII
 188 CYTB
 214 CYTB
  18 D-loop
 188 ND1
 188 ND2
 188 ND3
 188 ND4
 188 ND4L
 188 ND5-0
 188 ND6
 214 atp6
```

Filtering out markers
---------------------
It turns out there are multiple markers for this order. Unfortunately, because FASTA 
records are multiple lines (and the exact number is unpredictable), we can't easily use
command line tools (like `grep`). Instead, we might write a little script, e.g. in 
[biopython](http://biopython.org):

```python
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
	seen = {}
	for sp in species:
		length = 0
		longest = None
		for record in species[sp]:
			if len(record.seq) > length:
				length = len(record.seq)
				longest = record
		print '>' + longest.description
		print longest.seq
```

Run as:

```shell
$ python fasta.py > Artiodactyla.COI-5P.fas
```

Multiple sequence alignment
---------------------------

FASTA files can be aligned, for example, with [muscle](https://www.drive5.com/muscle/):

```shell
$ muscle -in Artiodactyla.COI-5P.fas -out Artiodactyla.COI-5P.muscle.fas
```

Resulting in an [alignment](Artiodactyla.COI-5P.muscle.fas), which is also a FASTA file. 

Alternatively, you might align with [mafft](https://mafft.cbrc.jp/alignment/software/), 
(or one of the many other multiple sequence alignment tools) which has additional 
functionality for more difficult markers (such as ITS):

```shell
$ mafft Artiodactyla.COI-5P.fas > Artiodactyla.COI-5P.mafft.fas
```

Resulting in this [file](Artiodactyla.COI-5P.mafft.fas). You can view both, for example, 
with this [web viewer](http://msa.biojs.net/app/). Are they different?

```shell
$ ls -la Artiodactyla.COI-5P.m*
-rw-r--r--  1 rutger.vos  NNM\Domain Users  400748 Jul 12 11:13 Artiodactyla.COI-5P.mafft.fas
-rw-r--r--@ 1 rutger.vos  NNM\Domain Users  400748 Jul 12 11:11 Artiodactyla.COI-5P.muscle.fas
```

Same number of bytes (so, same number of inserted gaps) but with different contents. Could
be capitalization, could be line folding, could be sequence order, or 
_actual differences in the alignment algorithms_. 

**How might we verify this further?**

Comparing different alignment results
-------------------------------------
Because the FASTA records in the files produced by `mafft` and `muscle` are in a different
order and the sequence data is capitalized differently, we can't easily compare the two
files. Here we sort the records, and capitalize the sequences, then write out to a
specified file format (e.g. `phylip`).

```python
import sys
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

with open(sys.argv[1], "rU") as handle:
    
    records = {}
    for seq in SeqIO.parse(handle, "fasta"):
        fields  = seq.description.split('|')
        seq.seq = seq.seq.upper()
        records[fields[0]] = seq
    
    aln = MultipleSeqAlignment([ records[key] for key in sorted(records.keys()) ])
    AlignIO.write(aln, sys.argv[2], sys.argv[3])
```

Usage:

```bash
$ python convert.py <infile> <outfile> <format>
```

Now we can compare the two alignments, e.g. using:

```bash
$ diff <mafft version> <muscle version>
```

It turns out there were no differences. Phew.


Bayesian evolutionary analysis by sampling trees (BEAST)
--------------------------------------------------------

[BEAST2](http://www.beast2.org/) is a modular system that can run many different types of 
analyses. The typical workflow usually goes:

1. Import data (e.g. a  FASTA alignment) into `beauti`, set up the analysis parameters, 
   possibly using a template
2. Start `beast filename.xml`, numerous output files (a.o. are log and tree files)
3. Inspect the log in `tracer` and run the analysis until the parameters all reach ESS>200
4. Summarize and interpret the results, e.g. build a consensus tree with `treeannotator`
   and visualize it with `figtree`

BEAST can read FASTA files, but it would be nice if the definition lines came out better
in trees, so we might relabel these:

```python
import sys
from Bio import SeqIO # sudo pip install biopython
with open(sys.argv[1], "rU") as handle:
    
    # Example: relabel sequences as Genus_species-ID
    for record in SeqIO.parse(handle, "fasta"):
        fields = record.description.split('|')
        name = fields[1].replace(' ', '_')
        print '>' + name + '-' + fields[0]
        print record.seq
```

Which gives us [this version](https://github.com/naturalis/mebioda/commit/76e9562db3f5ce1a8140f73f0b57d34e56e63b42)
to import in `beauti`, resulting in the [input file](BEAST/Danaus.mafft.xml).

Running a BEAST analysis
------------------------

![](BEAST/tracer.png)

If we run the [input file](BEAST/Danaus.mafft.xml) for 10*10^6 generations, the 
[log](BEAST/Danaus.log) file shows in tracer that all the parameters have been 
sufficiently sampled. If we compute a consensus, this is the result:

![](BEAST/Danaus.consensus.trees.png)
