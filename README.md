# Introduction
The StriDe Assembler integrates string and de Bruijn graph by decomposing reads within error-prone regions, while extending paire-end read into long reads for assembly through repetitive regions. The entire implementation is done by revising Simpson's [SGA][1] components for our own purpose, porting Li's [ropebwt2][2] for FM-index construction, and adding key components of this assembler. 

# Compile by yourself
To compile StriDe assembler in your specific environment, type 

      1. ./autogen.sh 
      2. ./configure
      3. make

An executable program called stride will be found under the StriDe folder.

# Execution
The program can be executed in all-in-one mode or step-by-step mode, depending on the commands specified.

All-in-one Commands (experimental):

      all	  Perform error correction, long-read generation, overlap computation, and assembly in one run

Step-by-step Commands:

      preprocess  filter and quality-trim reads
      index       build FM-index for a set of reads
      correct     correct sequencing errors in reads 
      fmwalk      merge paired reads into long reads via FM-index walk
      filter      remove redundant reads from a data set
      overlap     compute overlaps between reads
      assemble    generate contigs from an assembly graph

For instance, given two pair-end reads data sets in fastq format (A_R1.fq, A_R2.fq), simply type

      stride all A_R1.fq A_R2.fq

The entire preprocess, index, correction, fmwalk/decomposition, ... will be performed.

Example of a Step-by-step script

      stride preprocess --discard-quality -p 1 A_R1.fq A_R2.fq -o reads.fa
      stride index -a ropebwt2 -t 30 reads.fa
      stride correct -a overlap -t 30 -k 31 -x 3 reads.fa -o READ.ECOLr.fasta
      stride index  -t 30 READ.ECOLr.fasta
      stride fmwalk -m 80 -M 95 -t 30 -L 32 -I 400 -k 31 -p READ.ECOLr READ.ECOLr.fasta
      cat READ.ECOLr.merge.fa READ.ECOLr.kmerized.fa >merged.fa
      stride index -t 30 merged.fa
      stride filter -t 30 --no-kmer-check merged.fa
      stride overlap -m 30 -t 30 merged.filter.pass.fa
      stride assemble -k 31 -t 3 -p READ.ECOLr -r 100 -i 200 merged.filter.pass.asqg.gz

[1]: https://github.com/jts/sga
[2]: https://github.com/lh3/ropebwt2
