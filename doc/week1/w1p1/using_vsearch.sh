#cutadapt installation
pip install --user --upgrade cutadapt

#First, we trim the bad quality ends using cutadapt with the command -q 20,20 which trim both the 3' ends and the 5' ends

mkdir out_trimmed

for fq in *.fastq 
do
	cutadapt -q 20,20 -o out_trimmed/"${fq%.fastq}_trimmed_ends.fastq" ${fq}
done
#-----------------------------------------------------------------------------------

#we will now trim the 3' primer:

for fq in out_trimmed/*.fastq
do 
cutadapt -a TCCTCCGCTTATTGATAGC -o "${fq/_ends/_primers}" ${fq} 
done 

#-----------------------------------------------------------------------------------
  #Then, the 5' primer:
  for fq in out_trimmed/*_trimmed_primers.fastq
do 
cutadapt -g GTGARTCATCGAATCTTTG -o "${fq/_primers/_primers2}" ${fq}
done

#-----------------------------------------------------------------------------------

#Using vsearch, we will do the the quality filtering of the sequences, the commands used are: fastq_qmax:specify the maximum quality score accepted when reading FASTQ files. the default values is 41. fastq_maxee_rate 1: we specify the expected error per base; --fastq_minlen 200 discards the sequences smaller than 200 bp; --fastq_maxns 0 we specify that we do no want any N base.
for fq in out_trimmed/*_trimmed_primers2.fastq
do
vsearch --fastq_qmax 46 --fastq_filter ${fq} --fastq_maxee_rate 1 --fastq_minlen 200 --fastq_maxns 0 --fastaout "${fq%.fastq}.fa"
done

#-----------------------------------------------------------------------------------
# derep_prefix:merge sequences with identical prefixes contained in filename. A short sequence identical to an initial segment (prefix) of another sequence is considered a replicate of the longer sequence; --sizeout: the number of occurrences (i.e. abundance) of each sequence is indicated at the end of their fasta header using the pattern. 
for fq in out_trimmed/*_trimmed_primers2.fa
do
vsearch -derep_prefix ${fq} --output "${fq/trimmed_primers2/uniques}" --sizeout 
done 

#-----------------------------------------------------------------------------------

#Sequence renaming is done to add sample prefixes to the sequence identifiers, so that we can later merge the sequences.
cd out_trimmed
for fq in *_uniques.fa
do 
vsearch --rereplicate ${fq} --relabel ${fq}seq --sizeout --output ${fq/uniques/relabeled}.fasta
done
#Now, we can combine all files in one
cat *relabelled.fa* > CPuniques.txt

#-----------------------------------------------------------------------------------
#--sortbysize; Sort by decreasing abundance the sequences contained in filename (missing abundance values are assumed to be ’;size=1’)
vsearch --sortbysize CPuniques.fasta --output CPuniquessort.fasta

#-----------------------------------------------------------------------------------
# OTU clustering
#--cluster_size: sorts sequences by decreasing abundance before clustering 
vsearch --cluster_size CPuniquessort.fasta --id 0.97 --sizein --sizeout --fasta_width 0 --uc all_clustered.uc --relabel OTU_ --centroids all_otu.fasta --otutabout alla_otu.txt

vsearch -usearch_global all_otu.fasta -db utax_reference_dataset_10.10.2017.fasta -id 0.7 -blast6out CPotus2017.m8 -strand both -maxaccepts 1 -maxrejects 256