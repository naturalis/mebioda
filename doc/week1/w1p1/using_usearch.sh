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
#usearch9 installation, usearch10 gives problems when installing, not only to me I asked for confirmation
ln -s usearch9.2.64_ilinux32 usearch9
sudo ln -s /home/irene/Documents/usearch9 /bin/usearch9
chmod +x usearch9
#-----------------------------------------------------------------------------------
for fq in *_trimmed_primers2.fastq
do
	usearch9 -fastq_filter ${fq} -fastq_maxee 1.0 -fastq_minlen 200 -fastaout "${fq%.fastq}.fa"
done #fatal error to check, wrote to Robert Edgar, waiting for the reply. 
#-----------------------------------------------------------------------------------

for fq in *_trimmed_primers2.fa
do
usearch9 -fastx_uniques ${fq} -sizeout -fastaout "${fq/trimmed_primers2/uniques}"
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

usearch9 -sortbysize CPuniques.fasta -fastaout CPuniquessort.fasta

#-----------------------------------------------------------------------------------

#clustering
usearch9 -cluster_otus CPuniquessort.fasta -otu_radius_pct 3 -otus CPotus.fa -uparseout CPout.up -relabel OTU_ 

usearch9 -usearch_global all_otu.fasta -db utax_reference_dataset_10.10.2017.fasta -id 0.7 -blast6out CPotus2017.m8 -strand both -maxaccepts 1 -maxrejects 256