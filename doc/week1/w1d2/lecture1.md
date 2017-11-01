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

    @identifier (note paired-end sequencing)
    Sequence
    +(Repeat of title line)
    Quality lines map phred scores to ASCII characters

Example:

    @FAKE0005
    ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
    +
    @ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~

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

The SAM/BAM format
------------------
- Format to represent (FASTQ) reads aligned to a reference sequence
- Textual (SAM) and binary representations (BAM)
- Accessed using tools such as samtools, picard, EMBOSS, (Bio::SamTools, Galaxy)

The VCF format
--------------
- Format for variants (SNPs, indels) computed from a SAM/BAM file
- Concise, good for resequencing projects, but lossy

BED files
---------

The GFF format
--------------