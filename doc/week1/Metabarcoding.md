- multiplex sequencing MiSeq at Macrogen (_preliminary, done by lab_)
- demultiplexed on adaptors at Macrogen (_preliminary, done by lab_)
- assemble the contigs for each read (by ID, in Mothur, _installed on VM_)
- trim the primers and low quality ends, locally in Geneious (_see if FASTX-toolkit can be used for this_)
- export the FASTQ from Geneious (_see if FASTX-toolkit can be used for this_)
- cleanup with USearch, throw out low Phred score reads, truncate to 200nt (also throw out <200)
- (option: resample to normalize)
- group identical reads into distinct, unique types per sample, retain counts per type
- rename sequences in Geneious, give prefix of the sample (name)
- merge samples into a single file, cluster OTUs, remove chimeric
- make OTU matrix, identify OTUs taxonomically (perhaps normalize by rarifaction)


