#!/bin/bash
# creatation: 2023-08-03

# Stop on error
set -e

### MapForPaired.sh - This script is used to map paied-end ChIP-seq data to the reference genome.


### parameters for mapping
threads=$1
assembly=$2
spike_in=$3
basename=$4
file1=$5
file2=$6

mkdir -p $basename
cd $basename
mv ../$file1 ./
mv ../$file2 ./

###1.adapter trimming
echo "trimming adaptor"
time=`date`
echo $time
cutadapt --times 1 -e 0.1 -O 3  -m 30 -q 25,25 -u 8 -a AGATCGGAAGAGC -A AGATCGGAAGAGC  -o ${basename}_R1.trim.fq.gz -p ${basename}_R2.trim.fq.gz $file1 $file2
echo "trimming adaptor finished "
time=`date`
echo $time

###2.mapping to genome by bowtie2
echo "Mapping by bowtie2"
time=`date`
echo $time
bowtie2 -X 1000 -p $threads -x /share/Genomes/${assembly}/bowtie2/${assembly} -1 ${basename}_R1.trim.fq.gz -2 ${basename}_R2.trim.fq.gz -S ${basename}.sam >> ${basename}_ChIP-seq_mapping_summary.txt 2>&1
sed -i '/^@PG/d' ${basename}.sam
samtools sort -@ 3 ${basename}.sam -o ${basename}_sorted.bam
samtools index ${basename}_sorted.bam
echo "" >> ${basename}_ChIP-seq_mapping_summary.txt
echo "Reads duplication statistics" >> ${basename}_ChIP-seq_mapping_summary.txt
#remove duplicates
picard MarkDuplicates INPUT=${basename}_sorted.bam OUTPUT=${basename}_rmdup.bam METRICS_FILE=${basename}.PCR_duplicates REMOVE_DUPLICATES=true 
grep "Warning" -v ${basename}.PCR_duplicates >> ${basename}_KAS-seq_mapping_summary.txt
samtools index ${basename}_rmdup.bam

###3.create bed file and bedGraph file
echo "" >> ${basename}_ChIP-seq_mapping_summary.txt
echo "Average length of DNA fragments" >> ${basename}_ChIP-seq_mapping_summary.txt
samtools view -h  ${basename}_rmdup.bam | python /share/home/zhongyiting/pipeline/ChIP-seq/src/SAMtoBED.py  -i - -o  ${basename}.bed -x -v >> ${basename}_ChIP-seq_mapping_summary.txt 2>&1
bedSort ${basename}.bed ${basename}.sort.bed
genomeCoverageBed -bg -i ${basename}.sort.bed -g /share/Genomes/${assembly}/Sequence/${assembly}.chrom.sizes > ${basename}.bg


###4.spike-in
#i.map to spike-in genome
bowtie2 -X 1000 -p $threads -x /share/Genomes/${spike_in}/bowtie2/${spike_in} -1 ${basename}_R1.trim.fq.gz -2 ${basename}_R2.trim.fq.gz -S ${basename}.spike.sam >> ${basename}_ChIP-seq_mapping_summary.txt 2>&1
sed -i '/^@PG/d' ${basename}.spike.sam
samtools sort -@ 3 ${basename}.spike.sam -o ${basename}_sorted.spike.bam
samtools index ${basename}_sorted.spike.bam
echo "" >> ${basename}_ChIP-seq_mapping_summary.txt
echo "Reads duplication statistics" >> ${basename}_ChIP-seq_mapping_summary.txt
picard MarkDuplicates INPUT=${basename}_sorted.spike.bam OUTPUT=${basename}_rmdup.spike.bam METRICS_FILE=${basename}.spike.PCR_duplicates REMOVE_DUPLICATES=true 
grep "Warning" -v ${basename}.spike.PCR_duplicates >> ${basename}_KAS-seq_mapping_summary.txt
samtools index ${basename}_rmdup.spike.bam

#ii.create bed file
echo "" >> ${basename}_ChIP-seq_mapping_summary.txt
echo "Average length of DNA fragments" >> ${basename}_ChIP-seq_mapping_summary.txt
samtools view -h  ${basename}_rmdup.spike.bam | python /share/home/zhongyiting/pipeline/ChIP-seq/src/SAMtoBED.py  -i - -o  ${basename}.spike.bed -x -v >> ${basename}_ChIP-seq_mapping_summary.txt 2>&1

###5.normalize and make bigwig file
nor=$(wc -l ${basename}.spike.bed | cut -d ' ' -f 1)
awk '{printf("%s\t%d\t%d\t%.2f\n",$1,$2,$3,$4*100000000/'${nor}')}' ${basename}.bg > ${basename}.norm.bg
bedSort ${basename}.norm.bg ${basename}.norm.sort.bg
bedGraphToBigWig ${basename}.norm.sort.bg /share/Genomes/${assembly}/Sequence/${assembly}.chrom.sizes ${basename}.norm.bw

###6.organize files and folders
rm -f *sam
rm -f ${basename}_sorted.bam.bai
rm -f ${basename}_sorted.bam
rm -f ${basename}_sorted.spike.bam.bai
rm -f ${basename}_sorted.spike.bam
rm -f ${basename}_rmdup.spike.bam.bai
rm -f ${basename}_rmdup.spike.bam
rm -f ${basename}.PCR_duplicates
rm -f ${basename}.spike.PCR_duplicates
rm -f ${basename}.bed 
rm -f ${basename}.norm.bg
mv $file1 ../
mv $file2 ../

cd ..
mkdir -p Bam_files
cd Bam_files
mv ../${basename}/${basename}_rmdup.bam ./
mv ../${basename}/${basename}_rmdup.bam.bai ./
cd ..

mkdir -p BedGraph_files
cd BedGraph_files
mv ../${basename}/${basename}.bg ./
mv ../${basename}/${basename}.norm.sort.bg ./
cd ..

mkdir -p Bed_files
cd Bed_files
mv ../${basename}/${basename}.sort.bed ./
mv ../${basename}/${basename}.spike.bed ./
cd ..

mkdir -p BigWig_files
cd BigWig_files
mv ../${basename}/${basename}.norm.bw ./
cd ..

mkdir -p Mapping_summary
cd Mapping_summary
mv ../${basename}/${basename}_ChIP-seq_mapping_summary.txt ./
cd ..
rm -r $basename

###7.prepare the track for UCSC
cd BedGraph_files
lines=$(wc -l ${basename}.norm.sort.bg | cut -d ' ' -f 1)
if [ $lines -le 50000000 ]; then
sed -i '1i\track type=bedGraph name='${basename}' description='${basename}' color=0,0,0' ${basename}.norm.sort.bg
gzip ${basename}.norm.sort.bg

elif [ $lines -gt 50000000 ]; then
l=$(expr $line + 1)
a=$(expr $l / 2)
head ${basename}.norm.sort.bg -n $a > ${basename}.1.norm.sort.bg
tail ${basename}.norm.sort.bg -n $a > ${basename}.2.norm.sort.bg
sed -i '1i\track type=bedGraph name='${basename}.1' description='${basename}.1' color=0,0,0' ${basename}.1.norm.sort.bg
sed -i '1i\track type=bedGraph name='${basename}.2' description='${basename}.2' color=0,0,0' ${basename}.2.norm.sort.bg
gzip ${basename}.1.norm.sort.bg
gzip ${basename}.2.norm.sort.bg
fi

echo "=== All done successfully."
