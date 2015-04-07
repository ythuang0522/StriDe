# Introduction
The StriDe Assembler integrates string and de Bruijn graph by decomposing reads within error-prone regions, while extending paire-end read into long reads for assembly through repetitive regions. The entire implementation is done by revising Simpson's SGA components for our own purpose, porting Li's ropebwt2 for FM-index construction, and adding key components of this assembler. 

# Compile
To compile StriDe assembler, simply type 

      1. ./autogen.sh 
      2. ./configure
      3. make

An executable program called stride will be found under the StriDe folder.

# Execution
The program can be executed in all-in-one mode or step-by-step mode, depending on the commands specified.

All-in-one Commands:

      all	  Perform error correction, long-read generation, overlap computation, and assembly in one run

Step-by-step Commands:

      preprocess  filter and quality-trim reads
      index       build FM-index for a set of reads
      correct     correct sequencing errors in reads 
      fmwalk      merge paired reads into long reads via FM-index walk
      filter      remove redundant reads from a data set
      overlap     compute overlaps between reads
      assemble    generate contigs from an assembly graph


For instance, given two pair-end reads data sets in fastq format (A_R1.fq, A_R2.fq, B_R1.fq, B_R2.fq), simply type

      stride all A_R1.fq A_R2.fq B_R1.fq B_R2.fq

The entire preprocess, index, correction, fmwalk/decomposition, ... will be performed.

