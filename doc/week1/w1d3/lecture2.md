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

Run the script ([fasta.py](fasta.py)) as follows:

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

# sys.argv is the list of command line arguments
with open(sys.argv[1], "rU") as handle:
    
    records = {}
    for seq in SeqIO.parse(handle, "fasta"):
        fields  = seq.description.split('|')
        seq.seq = seq.seq.upper()
        records[fields[0]] = seq
    
    aln = MultipleSeqAlignment([ records[key] for key in sorted(records.keys()) ])
    AlignIO.write(aln, sys.argv[2], sys.argv[3])
```

Run the script ([convert.py](convert.py)) as follows:

```bash
$ python convert.py <infile> <outfile> <format>
```

Now we can compare the two alignments, e.g. using:

```bash
$ diff <mafft version> <muscle version>
```

It turns out there were no differences. Phew.


Exercise: generalizing the marker filter
----------------------------------------
The first python script we saw ([fasta.py](fasta.py)) extracts from the file 
`Artiodactyla.fas` for each species the longest sequence for the marker `COI-5P`. The 
script is useful, but only for a very specific case because the file name and the marker
name are "hard-coded" into the program. This is bad form: 

- What if the file has a different name? 
- What if we wanted to filter out all the markers and write them to separate files?

We can use some of the magic we saw in the second script ([convert.py](convert.py)) to
generalize the filter script, by using the special variable `sys.argv`. This variable
is a list. The first element in the list (`sys.argv[0]`) holds the name of the script.
The following elements hold the arguments that are passed to the script on the command
line. For example, if you run a script like this:

```bash
$ python script.py foo bar
```

Then `sys.argv` will hold:

0. `script.py`
1. `foo`
2. `bar`

*Assignment*: Open the script [fasta.py](fasta.py) in a text editor (like "editpad",
"notepad++", "bbedit", "gedit"). Modify the script so that the first argument is the
name of the FASTA file, and the second argument is the marker name. Save the script as
`fasta-filter.py`.