#seg=NS
miseq=c1
minion=flu-6-4-rebasecalled.pass
ref=all-fludb-H1N1-H3N2-FluB.fasta
subtype="H3N2"

for seg in NS NA NP M PB1 PB2 PA HA
do
# sort the segments by their hits
sort -nrk4 $minion.all-fludb.blast.seg-counts.tsv | grep -P $subtype'\t'$seg | head | cut -f1 > $minion-top-$subtype-$seg.tsv

# get the fasta sequences
grep --no-group-separator -A1 -F -f $minion-top-$subtype-$seg.tsv $ref > $minion-top-$subtype-$seg.fasta

# map dat fool
bowtie2-build $minion-top-$subtype-$seg.fasta $minion-top-$subtype-$seg.fasta
bowtie2 -p 40 --very-sensitive-local -x $minion-top-$subtype-$seg.fasta  -1 $miseq.r1.fastq.gz -2 $miseq.r2.fastq.gz | samtools view -bS - | samtools sort - miseq-$minion-top-$subtype-$seg && samtools index miseq-$minion-top-$subtype-$seg.bam

# make consensus
for i in $(samtools idxstats miseq-$minion-top-$subtype-$seg.bam | cut -f1 | grep -v '*'); do python ~/custom-scripts/bamfile_consensus_generate/bamfile_consensus_generate.py miseq-$minion-top-$subtype-$seg.bam $i ; done > miseq-$minion-top-$subtype-$seg.consensus.fasta

# map minion data
/opt/bioinformatics-software/bwa-0.7.12/bwa index $minion-top-$subtype-$seg.fasta
/opt/bioinformatics-software/bwa-0.7.12/bwa mem -x ont2d -t 40 $minion-top-$subtype-$seg.fasta $minion.fastq | samtools view -bS - | samtools sort - $minion-top-$subtype-$seg && samtools index $minion-top-$subtype-$seg.bam

# make minion consensus
for i in $(samtools idxstats miseq-flu-11-09-top-$subtype-$seg.bam | cut -f1 | grep -v '*'); do python ~/custom-scripts/bamfile_consensus_generate/bamfile_consensus_generate.py $minion-top-$subtype-$seg.bam $i; done > $minion-top-$subtype-$seg.consensus.fasta
done