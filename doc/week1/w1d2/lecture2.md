Genomes in biodiversity research
================================

Adaptor sequences
-----------------

![](lecture2/fragmentation_and_ligation.png)

Read data may still have **adaptor sequences** that were ligated to the fragments during
library preparation.

Clipping adaptors
-----------------

![](lecture2/fragmentsize.png)

Depending on the sequencing platform, insert size, and additional services provided by 
the sequencing lab, the reads may already be sorted by adaptors (which are then clipped 
off) or you may have to do this yourself. 

Effect of adaptor clipping
--------------------------

**Sturm M, C Schroeder & P Bauer**, 2016. SeqPurge: highly-sensitive adapter trimming for 
paired-end NGS data. _BMC Bioinformatics._ **17**: 208.
doi:[10.1186/s12859-016-1069-7](http://doi.org/10.1186/s12859-016-1069-7)

![](lecture2/clipping.png)

There are many tools available for adaptor clipping. Some are faster than others, and they
all affect downstream analysis time differently.

Clipping primers
----------------

![](lecture2/libstructure.png)

**Sidenote about amplicon sequencing**

- In _amplicon_ sequencing, the fragment will have been ligated with a primer as well as 
  an adaptor sequence. 
- This allows for more samples to be multiplexed because the number of combinations then 
  becomes _n adaptors_ * _n primers_ 
- And you probably don't need the amount of coverage on a single marker that a 
  non-multiplexed run would give you anyway
- However, platform vendors cannot de-multiplex automatically (because they know their
  own adaptors, but not _your_ primers), and with degenerate primers you'd have to do
  [fuzzy matching](https://github.com/naturalis/fastq-simple-tools/blob/master/script/splitfastq.pl#L128) 
  against their sequences

Quality trimming
----------------

![](lecture2/fastqc.png)

The Phred quality of called bases is lower:

- At higher base positions in the read (i.e. near the "end")
- The first few bases
- Along homopolymers

Hence, read trimming by quality is a common procedure, for which, again, many tools are
available.

Additional filtering
--------------------

![](lecture2/chimera.gif)

- A (crude) proxy for a read possibly being chimeric is that it occurs only once, i.e. as
  a singleton, which you might therefore filter out. This is especially the case in
  genome, but not in amplicon sequencing (because the chimera might be PCR'ed).
- On platforms that have variable length reads, you might want to filter out all reads
  below a threshold length
- When paired-end sequencing, one of the two 'ends' might be filtered out, in which case
  you might filter out the opposite end as well

Overlap-layout-consensus assembly
---------------------------------

- In _Sanger_ sequencing assembly, a _graph_ is constructed where every _vertex_ is a 
  read, and _edges_ connect the reads that overlap
- Subsequently, the _Hamiltonian path_, which visits every _vertex_ (read) exactly once
  is searched for, and the consensus sequence along the path is computed
- However, this _overlap-layout-consensus_ approach is hard to solve (computer scientists
  call this _NP-complete_)

![](lecture2/overlap-layout-consensus.png)

An alternative way to traverse the graph
----------------------------------------

- Computing the _Hamiltonian path_ is too computationally intensive ("NP-complete") for 
  HTS data 
- Another approach to traverse graphs is along the _Eulerian path_, where every _edge_
  (instead of vertex) is visited exactly once
- The _Eulerian path_ can be traversed in _linear time_ (computer scientists notate this
  as "_O_(|_E_|)")
- However, for this _Eulerian path_ to exist, either zero or two _vertices_ may exist with 
  odd _degree_, so the Seven Bridges of Königsberg, which Euler studied, don't form a path
  (and the overlap graph might also not)

![](lecture2/koenigsberg.png)

Making the graph amenable to Eulerian traversal
-----------------------------------------------

- In a _De Bruijn_ graph, all vertices have even degree
- A _De Bruijn_ graph connects symbolic sequence data such that every vertex is a sequence
  string of length _k_ (a "_k-mer_") that is connected to other such vertices if the 
  the sequences are identical along the substring of length _k_-1, i.e. the sequences
  are shifted one step relative to one another, which creates _directedness_ in the
  graph
- The simplest cases, with binary sequence data, are shown for _k_ = 1..3 (imagine what
  this would look like for four symbols):

![](lecture2/DeBruijn-as-line-digraph.svg)

HTS sequence data and _k-mers_
------------------------------

- The _De Bruijn_ graph presupposes that for every _k-mer_ there are neighbours that 
  overlap with it along the substring _k-1_
- However, read data does not naturally come out meeting that assumption: there are biases
  causing reads to overlap more irregularly
- However, the reads can be re-processed to some smaller size _k_ that are shifted
  relative to one another (note that this collapses the duplicates that are then created):

![](lecture2/K-mer-example.png)

Optimal _k-mer_ size
--------------------

- Smaller _k_ makes for a smaller graph, but at the cost of more duplicate _k-mers_ 
  collapsed, causing information loss and an inability to cross over repeat regions
- Higher _k_ is more memory intensive, and may increase the risk of _k-mers_ having no
  neighbours, causing short contigs



Taxonomic diversity: tomato relatives
-------------------------------------

Functional diversity: snake venom
---------------------------------

Spatial diversity: phytophagous insects
---------------------------------------

