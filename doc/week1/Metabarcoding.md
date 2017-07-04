- multiplex sequencing MiSeq at Macrogen
- demultiplexed on adaptors at Macrogen
- assemble the contigs for each read (by ID, in Mothur)
- trim the primers and low quality ends, locally in Geneious
- export the FASTQ from Geneious
- cleanup with USearch, throw out low Phred score reads, truncate to 200nt (also throw out <200)
- (option: resample to normalize)
- group identical reads into distinct, unique types per sample, retain counts per type
- rename sequences in Geneious, give prefix of the sample (name)
- merge samples into a single file, cluster OTUs, remove chimeric
- make OTU matrix, identify OTUs taxonomically (perhaps normalize by rarerifaction)


