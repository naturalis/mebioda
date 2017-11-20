#cutadapt installation
pip install --user --upgrade cutadapt

#First, we trim the bad quality ends using cutadapt with the command -q 20,20 which trim both the 3' ends and the 5' ends

mkdir out_trimmed

for fq in *.fastq 
do
	cutadapt -q 15,15 -o out_trimmed/"${fq%.fastq}_trimmed_ends.fastq" ${fq}
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

#Using vsearch, we are going to do quality filtering of the sequences, the commands used are: fastq_qmax:specify the maximum quality score accepted when reading FASTQ files. the default values is 41. fastq_maxee_rate 1: we specify the expected error per base; --fastq_trunclen 200, all the sequences will have the same length: 200 bp; 
for fq in out_trimmed/*_trimmed_primers2.fastq
do
vsearch --fastq_qmax 46 --fastq_filter ${fq} --fastq_maxee_rate 1 --fastq_trunclen 200 --fastaout "${fq%.fastq}.fa"
done

#-----------------------------------------------------------------------------------
# derep_fulllength:merge sequences with identical prefixes contained in filename. A short sequence identical to an initial segment (prefix) of another sequence is considered a replicate of the longer sequence; --sizeout: the number of occurrences (i.e. abundance) of each sequence is indicated at the end of their fasta header using the pattern. Sequence renaming is done to add sample prefixes to the sequence identifiers, so that we can later merge the sequences.
cd out_trimmed
for fq in *_trimmed_primers2.fa
do
vsearch --derep_fulllength  ${fq} --output "${fq/trimmed_primers2/uniques}" --relabel "${fq/trimmed_primers2/seq}" --sizeout 
done 

#merge all the samples
cat *uniques.fa* > CPuniques.fasta

#Dereplicate across samples and remove singletons
vsearch --derep_fulllength CPuniques.fasta --minuniquesize 2 --sizein --sizeout --fasta_width 0 --uc all.derep.uc --output all.derep.fasta
#-----------------------------------------------------------------------------------

#Precluster at 98% before chimera detection
 vsearch --cluster_size all.derep.fasta --id 0.98 --strand plus --sizein --sizeout --fasta_width 0 --uc all.preclustered.uc --centroids all.preclustered.fasta

#De novo chimera detection
 vsearch --uchime_denovo all.preclustered.fasta --sizein --sizeout --fasta_width 0 --nonchimeras all.denovo.nonchimeras.fasta 

#Reference chimera detection
 vsearch --uchime_ref all.denovo.nonchimeras.fasta --db utax_reference_dataset_10.10.2017.fasta --sizein --sizeout --fasta_width 0 --nonchimeras all.ref.nonchimeras.fasta

#-----------------------------------------------------------------------------------

#Extract all non-chimeric, non-singleton sequences, dereplicated
perl map.pl all.derep.fasta all.preclustered.uc all.ref.nonchimeras.fasta > all.nonchimeras.derep.fasta

#Extract all non-chimeric, non-singleton sequences in each sample
perl map.pl CPuniques.fasta all.derep.uc all.nonchimeras.derep.fasta > all.nonchimeras.fasta

#Cluster at 97% and relabel with OTU_n, generate OTU table
 vsearch --cluster_size all.nonchimeras.fasta --id 0.97  --strand plus --sizein --sizeout --fasta_width 0 --uc all.clustered.uc --relabel OTU_ --centroids all.otus.fasta --otutabout all.otutab.txt
vsearch --usearch_global all.otus.fasta -db utax_reference_dataset_10.10.2017.fasta -id 0.7 -blast6out CPotus.m8 -strand both -maxaccepts 1 -maxrejects 256
