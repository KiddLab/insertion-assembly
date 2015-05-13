# insertion-assembly


This repository is contains python scripts for a pipeline for assembling
transposable element (i.e., Alu) calls using Illumina short read data.  The pipeline
utilizes the CAP3 assembler, and is based on initial detection of candidate insertions 
based on the program RetroSeq.

## Required Components
The following programs must be in the *path* and be *executable*:
bedtools, CAP3, RepeatMasker, miropeats, ps2pdf, samtools, blat,
stretcher (from the EMBOSS tool set),
fastalength (from the exonerate utilities),
fastq_to_fasta.py (from http://nebc.nerc.ac.uk/downloads/scripts/parse/fastq_to_fasta.py)

The pipeline is designed to work with results output from RetroSeq for discovery of mobile
element insertions, but can easily be modified for local assembly of other variants or using
other programs for initial candidate discovery.
 
## Example usage description

To illustrate program usage, we consider analysis of Alu calls identified using RetroSeq
based on 2x100 bp Illumina WGS data from the CHM1 mole sample ( [Chaisson et al,](http://www.ncbi.nlm.nih.gov/pubmed/25383537) 
SRA: [SRX652547](http://www.ncbi.nlm.nih.gov/sra/SRX652547[accn])  )

