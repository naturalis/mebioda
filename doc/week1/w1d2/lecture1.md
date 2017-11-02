(High-throughput) DNA sequencing
================================

Sequencing techniques
---------------------
- (Meta)barcoding of specific markers
- Probing for ultraconserved regions
- RNA sequencing
- etc.

Sequencing platforms
--------------------
- Old school: Sanger sequencing
- Short read standard: Illumina HiSeq
- Long reads: PacBio
- Longer reads: Nanopore MinION

File formats
------------
- FASTA
- FASTQ
- SAM/BAM
- VCF
- BED
- GFF

The FASTQ format
----------------
- Sequential format returned by most HTS platforms
- The initial mountain of data to deal with
- Includes base calling quality scores

Record layout:

1. `@`+identifier (note paired-end sequencing)
2. Sequence data, IUPAC single character nucleotides
3. +
4. Quality lines map phred scores to ASCII characters

Example:

    @FAKE0005
    ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
    +
    @ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~

Quality (phred) scores
----------------------

Different platforms map phred scores in different ways to ASCII:
- sanger: 33..126
- solexa: 59..126
- illumina: 64..126

|code|char|code|char|code |char|code |char|
|----|----|----|----|-----|----|-----|----|
| 33 | !  | 57 | 9  | 81  | Q  | 105 | i  |
| 34 | "  | 58 | :  | 82  | R  | 106 | j  |
| 35 | #  | 59 | ;  | 83  | S  | 107 | k  |
| 36 | $  | 60 | <  | 84  | T  | 108 | l  |
| 37 | %  | 61 | =  | 85  | U  | 109 | m  |
| 38 | &  | 62 | >  | 86  | V  | 110 | n  |
| 39 | '  | 63 | ?  | 87  | W  | 111 | o  |
| 40 | (  | 64 | @  | 88  | X  | 112 | p  |
| 41 | )  | 65 | A  | 89  | Y  | 113 | q  |
| 42 | *  | 66 | B  | 90  | Z  | 114 | r  |
| 43 | +  | 67 | C  | 91  | [  | 115 | s  |
| 44 | ,  | 68 | D  | 92  | \  | 116 | t  |
| 45 | -  | 69 | E  | 93  | ]  | 117 | u  |
| 46 | .  | 70 | F  | 94  | ^  | 118 | v  |
| 47 | /  | 71 | G  | 95  | _  | 119 | w  |
| 48 | 0  | 72 | H  | 96  | \` | 120 | x  |
| 49 | 1  | 73 | I  | 97  | a  | 121 | y  |
| 50 | 2  | 74 | J  | 98  | b  | 122 | z  |
| 51 | 3  | 75 | K  | 99  | c  | 123 | {  |
| 52 | 4  | 76 | L  | 100 | d  | 124 | \| |
| 53 | 5  | 77 | M  | 101 | e  | 125 | }  |
| 54 | 6  | 78 | N  | 102 | f  | 126 | ~  |
| 55 | 7  | 79 | O  | 103 | g  |     |    |
| 56 | 8  | 80 | P  | 104 | h  |     |    |

The SAM/BAM/CRAM format
-----------------------
- Format to represent (FASTQ) reads aligned to a reference sequence
- Textual (SAM) and binary representations (BAM)
- Accessed using tools such as samtools, picard, EMBOSS, (Bio::SamTools, Galaxy)

The VCF/BCF format
------------------
- Format for variants (SNPs, indels, microsats) computed from a SAM/BAM/CRAM file
- Concise, good for resequencing projects, but lossy

```
##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM  POS     ID          REF ALT     QUAL    FILTER  INFO                                FORMAT      NA00001         NA00002         NA00003
20      14370   rs6054257   G   A       29      PASS    NS=3;DP=14;AF=0.5;DB;H2             GT:GQ:DP:HQ 0|0:48:1:51,51  1|0:48:8:51,51  1/1:43:5:.,.
20      17330   .           T   A       3       q10     NS=3;DP=11;AF=0.017                 GT:GQ:DP:HQ 0|0:49:3:58,50  0|1:3:5:65,3    0/0:41:3
20      1110696 rs6040355   A   G,T     67      PASS    NS=2;DP=10;AF=0.333,0.667;AA=T;DB   GT:GQ:DP:HQ 1|2:21:6:23,27  2|1:2:0:18,2    2/2:35:4
20      1230237 .           T   .       47      PASS    NS=3;DP=13;AA=T                     GT:GQ:DP:HQ 0|0:54:7:56,60  0|0:48:4:51,51  0/0:61:2
20      1234567 microsat1   GTC G,GTCT  50      PASS    NS=3;DP=9;AA=G                      GT:GQ:DP    0/1:35:4        0/2:17:2        1/1:40:3
```

1. a good simple SNP 
2. a possible SNP that has been filtered out because its quality is below 10
3. a site at which two alternate alleles are called, with one of them (T) being ancestral 
   (possibly a reference sequencing error)
4. a site that is called monomorphic reference (i.e. with no alternate alleles)
5. a microsatellite with two alternative alleles, one a deletion of 2 bases (TC), and the 
   other an insertion of one base (T). 
   
Genotype data are given for three samples, two of which are phased and the third unphased, 
with per sample genotype quality, depth and haplotype qualities (the latter only for the 
phased samples) given as well as the genotypes. The microsatellite calls are unphased.

BED files
---------

The GFF format
--------------