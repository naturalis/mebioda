Big data phylogenetics
======================

Phylogenetic data
-----------------
- Plays a role in a variety of different analyses
- Both intra- and interspecific (gene trees versus species trees)
- Need for metadata and annotation, but hard to represent as tables

The Newick / New Hampshire format
---------------------------------

![](phylogeny.png)

Newick representation:

```
(((A,B),(C,D)),E)Root;
```

- Optionally has branch lengths after each tip or node, as `:0.002321`
- ['Invented' on a napkin at Newick's Lobster House in Durham, New Hampshire, in 1986.](http://evolution.genetics.washington.edu/phylip/newicktree.html)
- Concise, but lacking all metadata

The Nexus format
----------------

![](phylogeny.png)

Nexus representation:

```
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
```

- Uses an extensible block structure that can also include character data (and other 
  things)
- Many different, mutually incompatible dialects that deviate from the original 
  [standard](NEXUS_final.pdf)
- More facilities for metadata (e.g. names of things, annotations for taxa)

New Hampshire eXtended
----------------------

PhyloXML
--------

NeXML
-----

Tabular representations
-----------------------

Phyloinformatic projects
------------------------
- Gene tree projects, e.g. TreeFam
- Many species tree projects, e.g. TreeBASE
- Big tree projects, e.g. OpenTree, ToLWeb, phylogenetic placement