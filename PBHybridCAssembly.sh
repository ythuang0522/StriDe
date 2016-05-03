R1=$1
R2=$2
PB=$3
SD=../StriDe/stride
k=750
ovl=749

startTime=`date "+%s"`

# preprocess short reads from StriDe GitHub:
# https://github.com/ythuang0522/StriDe
$SD preprocess --discard-quality -p 1 $1 $2 -o reads.fa
$SD index -a ropebwt2 -t 30 reads.fa
$SD correct -a overlap -t 30 -k 31 -x 3 reads.fa -o READ.ECOLr.fasta
$SD index -t 30 READ.ECOLr.fasta

# pacbio hybrid correction
$SD pbhc -p READ.ECOLr -t 30 -s 21 -k 31 -M 91 -L 256 $3
$SD index -a ropebwt2 -t 30 PB.PBHybridCor.fa

# decompose
$SD fmwalk -a validate -t 30 PB.PBHybridCor.fa -m $ovl -k $k -L 128
cat PB.PBHybridCor.origin.fa PB.PBHybridCor.kmerized.fa > merged.fa

# filter redundant reads
$SD index -a ropebwt2 -t 30 merged.fa
$SD filter -t 30 merged.fa

# LSSF overlap
$SD overlap -m $ovl -e 0.05 -l 50 merged.filter.pass.fa -t 30

# assembly
# -i is the median or N50 PacBio read length
$SD asmlong -i 13000 -p PB.PBHybridCor merged.filter.pass.asqg.gz

endTime=`date "+%s"`
elapsed=$((endTime-startTime))
printf "Total time : %i:%02i:%02i\n" $((elapsed/3600)) $(((elapsed/60)%60)) $(($elapsed%60))
