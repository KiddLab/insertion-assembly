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

# Initial preperation

We begin with a BAM file of Illumina reads mapped and processed using standard approaches.
RetroSeq is then used discover candidate Alu insertions, as described.  Given the problems
associated with calling mobile element insertions in regions next to existing elements of the
same type in the reference, we first filter out calls located within 500 bp of reference elements.
This is done using bed tools.

```
**Inputs**
CHM1_lib1.SINE.calls.out.PE.vcf  --> output call file from RetroSeq
hg19.RM.Alu.sites.sorted.slop500.bed  --> location of Alu in refernece, expanded by
500bp using bedtools slop
**Outputs**
CHM1_lib1.SINE.calls.out.PE.notRef500.vcf --> filtered VCF, without candidates near reference
elements
```

Command:
intersectBed -v -a ../CHM1_lib1.SINE.calls.out.PE.vcf \
-b hg19.RM.Alu.sites.sorted.slop500.bed \
> CHM1_lib1.SINE.calls.out.PE.notRef500.vcf

