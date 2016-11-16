# insertion-assembly


This repository is contains python scripts for a pipeline for assembling
transposable element (i.e., Alu) calls using Illumina short read data.  The pipeline
utilizes the CAP3 assembler, and is based on initial detection of candidate insertions 
based on the program RetroSeq.

For more information, please see: 

Wildschutte JH, Baron A, Diroff NM, Kidd JM. Discovery and characterization of Alu repeat sequences via precise local read
assembly. Nucleic Acids Res. 2015 Dec 2;43(21):10292-307. doi: 10.1093/nar/gkv1089. Epub 2015 Oct 25.

https://www.ncbi.nlm.nih.gov/pubmed/26503250

## Required Components
The following programs must be in the *path* and be *executable*:
bedtools, CAP3, RepeatMasker, miropeats, ps2pdf, samtools, blat,
stretcher (from the EMBOSS tool set),
fastalength (from the exonerate utilities),
fastq_to_fasta.py (from http://nebc.nerc.ac.uk/downloads/scripts/parse/fastq_to_fasta.py)

The pipeline is designed to work with results output from RetroSeq for discovery of mobile
element insertions, but can easily be modified for local assembly of other variants or using
other programs for initial candidate discovery.

**Note**
Manual changes must be made to RepeatMasker options (for example, for different species)
and to output parsing (for example, LINE vs SINE element matching) at appropriate places
in some scripts.  Also, for proper execution on some systems may require changes to the "scriptDir"
parameter in some files.  See code for details and for where changes are required.

 
## Example usage description

To illustrate program usage, we consider analysis of Alu calls identified using RetroSeq
based on 2x100 bp Illumina WGS data from the CHM1 mole sample ( [Chaisson et al,](http://www.ncbi.nlm.nih.gov/pubmed/25383537) 
SRA: [SRX652547](http://www.ncbi.nlm.nih.gov/sra/SRX652547[accn])  )

## Initial Preparation

We begin with a BAM file of Illumina reads mapped and processed using standard approaches.
RetroSeq is then used discover candidate Alu insertions, as described.  Given the problems
associated with calling mobile element insertions in regions next to existing elements of the
same type in the reference, we first filter out calls located within 500 bp of reference elements.
This is done using bed tools.

```
Inputs:
CHM1_lib1.SINE.calls.out.PE.vcf  --> output call file from RetroSeq
hg19.RM.Alu.sites.sorted.slop500.bed  --> location of Alu in refernece, expanded by
500bp using bedtools slop

Outputs:
CHM1_lib1.SINE.calls.out.PE.notRef500.vcf --> filtered VCF, without candidates near reference
elements

Command:
intersectBed -v -a ../CHM1_lib1.SINE.calls.out.PE.vcf \
-b hg19.RM.Alu.sites.sorted.slop500.bed \
> CHM1_lib1.SINE.calls.out.PE.notRef500.vcf
```

Next, select sites with RetroSeq support level >=6 for assembly

```
Inputs:
CHM1_lib1.SINE.calls.out.PE.notRef500.vcf --> filtered VCF, without candidates near reference
elements
Outputs:
CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678 --> only calls with support level >=6

Command:
python select-retroseq-level678.py \
--in CHM1_lib1.SINE.calls.out.PE.notRef500.vcf \
--out CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678
```

## Assembly Pipeline

In step 1, soft clipped reads in region of each call are gathered

```
Inputs:
CHM1_lib1.markdup.bam --> BAM file used for RetroSeq 
Outputs:
CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678.softclip --> file of reads with softclipping
that are in each call region

Command:
python 01_gather-soft-clipped-persample.py \
--vcf CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678 \
--bam CHM1_lib1.markdup.bam \
--out CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678.softclip

```

In step 2, we gather all the read pairs that support candidate calls.  To get both
ends of a read the BAM file must be traversed entirely, which may take some time.

```
Inputs:
CHM1_lib1.SINE.discovery.tab --> file from the *discovery* step of RetroSeq
Outputs:
CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678.w500.sam --> SAM file of read pairs
from the soft clip and original support files from RetroSeq.

Command:
python 02_gather-seq-sam-persample.py \
--vcf CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678 \
--bam  CHM1_lib1.markdup.bam \
--softclip CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678.softclip \
--discover  CHM1_lib1.SINE.discovery.tab  \
--out CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678.w500.sam
```

In step 3, reads are matched with the call they support using bedtools windowBed.  
Additional formating also occurs prior to assembly.

```
Inputs:
Files from previous steps
Outputs:
CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678.reads.seq --> file of sequence information
listing name, read number, quality, and sample name
CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678.reads.names --> read names and coordinates

Command:
python 03_combine-sam-and-evidence.py \
--vcf CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678 \
--softclip CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678.softclip \
--discover  CHM1_lib1.SINE.discovery.tab   \
--sam CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678.w500.sam \
--outseq CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678.reads.seq \
--outnames CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678.reads.names \
--sample CHM1
```

In step 4, reads supporting each candidate call are assembled using CAP3
```
Inputs:
CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678.reads.names.w500 --> output of windowbed 
as executed in step 3
CHM1/initial_assembly --> directory for initial assembly
Outputs:
initial assemblies are reported under the directory indicated.  Information is also written
to the indicated log file.

Command:
python 04_run-cap3-assembly.py \
--vcf CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678 \
--seq CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678.reads.seq \
--namesintersect CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678.reads.names.w500 \
--outdirbase CHM1/initial_assembly \
--log logs/CHM1_lib1.678.assem.log
```

Contigs from CAP3 are linked together based on CAP3 reported linking informaiton to form
scaffolds.  Contigs are separated by a run of 'N'.  RepeatMasker is then ran on the assembled
sequences.

```
Command:
python 05_combine-contigs-based-on-links.py \
--vcf CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678 \
--outdirbase CHM1/initial_assembly 
```

## Parsing of Assembly Results

Results of the assembly are analyzed to identify insertions with assembled sequence supporting
the givne insertion type.  These candidate are then further analyzed in the breakpoint 
analysis step.

In step 6, scaffolds are analyzed for matches to the element of interest.

```
Outputs:
CHM1_lib1.assembly_SINE.txt --> summary table of repeat hits for element type of interest
for each assembled sequence

Command:
python 06_parse-sine-cap3-only-assembly.py \
--vcf CHM1_lib1.SINE.calls.out.PE.notRef500.vcf.sel678 \
--outdirbase CHM1/initial_assembly \
--out CHM1_lib1.assembly_SINE.txt
```

In step 7, candidates insertion sties with 30 bp of sequence on each end are selected

```
Command:
python 07_get-has-30-both.py \
--in CHM1_lib1.assembly_SINE.txt
```

In step 8, individual assembled contigs passing criteria are selected

```
Command:
python 08_select-fragment-that-is-full.py \
--in CHM1_lib1.assembly_SINE.txt.30both
```

In step 9, the longest hit is chosen for each site

```
Command:
python 09_select-longest.py \
--in CHM1_lib1.assembly_SINE.txt.30both.sel
```

In step 10, an additional check is performed to ensure that at least one end of the 
RepeatMasked position does not end in a gap in the assembled sequence

```
Command:
python 10_check-not-both-ends-gap.py \
--outdirbase CHM1/initial_assembly \
--in CHM1_lib1.assembly_SINE.txt.30both.sel.longest
```

## Alignment to Reference Genome and Breakpoint Determination

Next, candidate insertion will be aligned to the reference and characterized. This
includes making of annotated images using miropeats for comparison, then ultimately
extraction of sequence and 3-way alignment to determine breakpoints

In step 11, the directory and files required for subsequent steps are setup and the 
orientaiton of the contig (i.e, direct or reverse complement) is determined.

```
Input:
CHM1_lib1.assembly_SINE.txt.30both.sel.longest.notbothgap --> selected contig for each site from
step 10
hg19.fa --> reference genome, indixed using samtools
initial_assembly/ --> directory of initial assembly results from above
parse/initial/ --> directory for initial files for reference comparison
Output:
required files are written in directories created under the indicated location.  Also,
new file summary of contig orientaiton and information is written to 
CHM1_lib1.assembly_SINE.txt.30both.sel.longest.notbothgap.getcontig

Command:
python 11_get-breakpoints-setup.py \
--in CHM1_lib1.assembly_SINE.txt.30both.sel.longest.notbothgap \
--ref bwa-0.5.9-index/hg19.fa \
--assembase initial_assembly/  \
--parsebase parse/initial/
```

In step 12, miropeats is used to get an initial breakpoint guess for subsequent steps.

```
Input:
parse/miropeats --> directory where miropeats analysis files will be created and used
Output:
parse/CHM1_lib1.assembly_SINE.txt.initial_miropeats --> table of initial miropeats results
and parse information.
For each candidate insertion, and file *.initialMR.annotated.pdf is also created giving
annotated version of miropeats comparison.

Command:
python 12_get-breakpoints-initialMR.py \
--ref bwa-0.5.9-index/hg19.fa \
--tmp tmp/ \
--assembase initial_assembly/  \
--miropeats_base parse/miropeats \
--in CHM1_lib1.assembly_SINE.txt.30both.sel.longest.notbothgap.getcontig \
--out parse/CHM1_lib1.assembly_SINE.txt.initial_miropeats
```

In step 13, the initial attempt at extracting sequencing and aligning the three breakpoint
regions to determine exact breakpoints is performed.  Results are summarized in new file and
PDF images for each site are constructed.


```
Output:
CHM1_lib1.assembly_SINE.txt.initial_miropeats.3wayalign and  --> table summarizing parse 
results
For each candidate insertion, and file *.3way.combined.annotated.pdf is also created giving
annotated version of miropeats comparison.

python 13_get-breakpoints-3way-align.py \
--ref bwa-0.5.9-index/hg19.fa \
--miropeats_base parse/miropeats \
--in parse/CHM1_lib1.assembly_SINE.txt.initial_miropeats
```

At this point, candidate breakpoints have been determined for each site.  They should then
be assessed for accuracy, sites dropped without support (for example, if variant is actually
flanked by gaps without evidence of continuity), or repeated with new sequences extracted
for breakpoint alignment (for example, if the automated miropeats analysis selected the wrong 
candidate breakpoint region, or if aligned sequence is too short to show clear breakpoint
pattern.

This can be facilitated using combine-pdf-images.py which makes a single (sorted) PDF of 
the individual images for review.

```
python combine-pdf-images.py \
--miropeats_base parse_2015-04_26/miropeats \
--in parse/CHM1_lib1.assembly_SINE.txt.initial_miropeats.3wayalign \
--out parse/combined-pdfs/CHM1_lib1.SINE.3wayalign.pdf
```

Sites can then be dropped, and changes updated iteratively using 14_update-3way-align.py,
which based on coded comments from manual review repeats breakpoint identification process
based on new sequences that are extracted for alignment.


```
Input:
parse/CHM1_lib1.assembly_SINE.txt.initial_miropeats.3wayalign.changes --> table giving
changes in sequence to extract for breakpoint alignment for sites requiring changes
Output:
parse/CHM1_lib1.assembly_SINE.txt.initial_miropeats.3wayalign.changes.a1 --> breakpoint table
for updated sites

Command:
python 14_update-3way-align.py \
--ref bwa-0.5.9-index/hg19.fa \
--miropeats_base parse_2015-04_26/miropeats \
--in parse/CHM1_lib1.assembly_SINE.txt.initial_miropeats.3wayalign.changes \
--out parse/CHM1_lib1.assembly_SINE.txt.initial_miropeats.3wayalign.changes.a1
```




