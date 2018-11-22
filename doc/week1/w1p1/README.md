Week 1 Step-by-step analysis of micorrhizal diversity
-----------------------------------------------------

In this practical, we are going to analyze some HTS data sets (MiSeq) for soil fungi in samples taken at 
different locations along a mountain slope. The following are the steps we take:

### Downloading our data

1. Create a directory where to store your data, e.g. `mkdir w1p1`
2. Go into the directory and download the data, e.g. `wget https://github.com/naturalis/mebioda/blob/master/doc/week1/w1p1/C1P1.fastq.gz?raw=true`
3. Do this for all the \*.fastq.gz files
4. Extract the data, e.g. `gunzip *.gz`

### Trimming the sequences

Our sequencing reads may contain base calls with low support. These typically occur especially on the
ends of the reads. We are going to remove them. In addition, there still are primer sequences on both
ends (i.e. 3' and 5'), which have to be removed. For all of these, we use a tool called
[cutadapt](https://cutadapt.readthedocs.io/en/v1.10/installation.html). 

1. Install cutadapt: `pip install --user --upgrade cutadapt`
2. Trim the bad quality ends:

```bash
mkdir out_trimmed
for fq in *.fastq; do
  cutadapt -q 20,20 -o out_trimmed/"${fq%.fastq}_trimmed_ends.fastq" ${fq}
done
```

3. Trim the 3' primer:

```bash
for fq in out_trimmed/*.fastq; do 
  cutadapt -a TCCTCCGCTTATTGATAGC -o "${fq/_ends/_primers}" ${fq} 
done
```

4. Trim the 5' primer

```bash
for fq in out_trimmed/*_trimmed_primers.fastq; do 
  cutadapt -g GTGARTCATCGAATCTTTG -o "${fq/_primers/_primers2}" ${fq}
done
```

### Analyzing the data

We are going to use the tool `vsearch`, first described in a paper by
[Rognes et al., 2016.](https://dx.doi.org/10.7717%2Fpeerj.2584). Part of the
analysis will be a continuation of data pre-processing steps, but in ways that
`cutadapt` can't do. There is additional quality filtering, and deduplication of
the reads. Then, the read IDs are given additional labels to identify the 
different soil samples. This is because we are going to merge the data and
cluster the sequences across all samples - but we still want to be able to
recognize where the sequences came from originally, and that information would
otherwise be lost.

Subsequently, the sequences are going to be clustered: the ones that are at least
97% similar are all going to end up in the same putative taxon (OTU). Once the
clustering is done, we can also figure out which of the sequences are weird 
singletons. Those are probably artefacts of the sequencing process (_chimeras_),
which we will remove. We then compare the remaining, plausible, OTUs with a 
reference database to see if we can recognize and classify the OTUs. This step
will produce the biologically relevant information that we need (though we will
process the output a bit more so that we can easily interpret it in a spreadsheet
program).

1. Install `vsearch` from its source code. The steps to take on the command line
   are [here](https://github.com/torognes/vsearch#download-and-install)
2. Do the additional quality filtering:

```bash   
for fq in out_trimmed/*_trimmed_primers2.fastq; do
  parameters="--fastq_qmax 46 --fastq_maxee_rate 1 --fastq_trunclen 200"
  vsearch $parameters --fastq_filter ${fq} --fastaout "${fq%.fastq}.fa"
done
```

   The parameters mean the following:
   
   - `--fastq_qmax`: specify the maximum quality score accepted when reading FASTQ files. 
   - `--fastq_maxee_rate 1`: we specify the expected error per base: 1 
   - `--fastq_trunclen 200`, all the sequences will be truncated at 200 bp; 

3. Deduplicate the reads at sample level and relabel with the sample number:

```bash
cd out_trimmed
for fq in *_trimmed_primers2.fa; do
  vsearch --derep_fulllength $fq --output ${fq/trimmed_primers2/uniques} --relabel ${fq/trimmed_primers2/seq} --sizeout --minsize 2
done 
```

   The parameters mean the following:
   
   - `--derep_fulllength`: Merge strictly identical sequences contained in filename.   
   - `--sizeout`: the number of occurrences (i.e. abundance) of each sequence is indicated 
     at the end of their fasta header.
   - `--relabel`: sequence renaming is done to add sample prefixes to the sequence identifiers, 
     so that we can later merge the sequences
        
4.  Merge the data across all samples:

```bash
rm *_ends.fastq *_primers.fastq *_primers2.fastq *_primers2.fa
cat *uniques.fa* > CPuniques.fasta
```

5. Cluster the reads at 97% similarity

```bash
vsearch \
--cluster_size CPuniques.fasta \
--id 0.97 \
--strand plus \
--sizein \
--sizeout \
--fasta_width 0 \
--uc CP_otus.uc \
--relabel OTU_ \
--centroids CP.otus.fasta \
--otutabout CP.otutab.txt
```

   The parameters mean the following:
   
   - `--cluster_size`: sorts sequences by decreasing abundance before
     clustering.
   - `--id`: reject the sequence match if the pairwise identity is lower than 0.97
   - `--fasta_width`: Fasta files produced by vsearch are wrapped (sequences are written on lines of integer
     nucleotides, 80 by default). To eliminate the wrapping, the value is set to zero.
   - `--centroids`: The centroid is the sequence that seeded the cluster (i.e. the first sequence of the cluster).     
   - `--otutabout`: Output an OTU table in the classic tab-separated plain text format as a matrix containing
     the abundances of the OTUs in the different samples.
 
6. Detect likely chimeras _de novo_

```bash
vsearch \
--uchime_denovo CP.otus.fasta \
--sizein \
--sizeout \
--fasta_width 0 \
--nonchimeras CP.otus.nonchim.fasta 
```

   - `--uchime_denovo`: Detect chimeras present in the fasta-formatted filename, without external references (i.e. de novo)       

7. Compare clusters with reference database

```bash
vsearch \
--usearch_global CP.otus.nonchim.fasta \
-db ../utax_reference_dataset_10.10.2017.fasta \
-id 0.7 \
-blast6out CPotus.m8 \
-strand both \
-maxaccepts 1 \
-maxrejects 256
```

   - `--usearch_global`: Compare target sequences (--db) to the fasta-formatted query sequences contained in
     filename, using global pairwise alignment.

### Post-processing

We now have the all the output files that we need. However, the `CPotus.m8` needs to be
processed in the following way:

- The fourth column contains the sequence length. We only want sequences longer than
  150 bp
- The first column contains both the OTU id and the cluster size, like this: 
  `OTU_1;size=9159;`. We just want to keep the `OTU_1` part, so before the first 
  semicolon.
- For each distinct OTU id, we want to keep the match with the highest percentage (there
  might be multiple matches for the same OTU).

This [script](filter.pl) performs these operations. It is executed like this:

```bash
perl ../filter.pl CPotus.m8 > CPotus.uniq.tsv
```

We now have two tables that each have OTU identifiers in the first column. One table
contains the taxonomic identification of these OTUs (i.e. their nearest matches in 
the Linnean taxonomy), the other table contains the abundances of these OTUs in the
different samples. We want to combine these by matching up the OTU identifiers. This
is a database operation called a [join](https://en.wikipedia.org/wiki/Relational_algebra#Joins_and_join-like_operators).
The Linux environment contains a standard tool that allows you to do this with 
tab-separated text files, producing output that consists of the joined records. 
Run it like this:

```bash
join CPotus.uniq.tsv CP.otutab.txt > CPotus.samples.tsv
 ```

### Interpretation

We now need to get the file `CPotus.samples.tsv` on our own computers. Depending on
the platform (Windows, Mac, Linux) there are different ways to do that. On
Windows you would use [putty](https://www.putty.org/). On Mac and Linux you use the
standard `scp` command. Once you have the file, import it in a spreadsheet program
such as Excel, as a tab-separated file. In the spreadsheet, we are going to summarize
our results:

- After the OTUs' abundance values you can calculate the sum, e.g `=SUM(a1:g1)`
- Next to the sum we are going to write an IF condition that will return the value 1 if the 
  abundance value is greater than 0, otherwise it will return 0. `=IF(e1>0,1,0)`
- Finally, we want the sum of each absence-presence column, so at the end of each column 
  calculate the sum.

