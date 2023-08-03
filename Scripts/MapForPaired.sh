#!/bin/bash
# creatation: 2023-08-03

# Stop on error
set -e

#mapping with bowtie2
file1 = $1
file2 = $2
basename = $3
threads = $4
bowtie2_index_path = $5
assembly = $6


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
bowtie2 -X 1000 -p $threads -x $bowtie2_index_path$assembly -1 ${basename}_R1.trim.fq.gz -2 ${basename}_R2.trim.fq.gz -S ${basename}.sam >> ${basename}_ChIP-seq_mapping_summary.txt 2>&1
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

###4.organize files and folders
rm -f *sam
rm -f ${basename}_sorted.bam.bai
rm -f ${basename}_sorted.bam
rm -f ${basename}.PCR_duplicates
rm -f ${basename}.bed 
mv $raw_fastq_read1 ../
mv $raw_fastq_read2 ../

cd ..
mkdir -p Bam_files
cd Bam_files
mv ../${basename}/${basename}_rmdup.bam ./
mv ../${basename}/${basename}_rmdup.bam.bai ./
cd ..

mkdir -p BedGraph_files
cd BedGraph_files
mv ../${basename}/${basename}.bg ./
cd ..

mkdir -p Bed_files
cd Bed_files
mv ../${basename}/${basename}.sort.bed ./
cd ..

mkdir -p Mapping_summary
cd Mapping_summary
mv ../${basename}/${basename}_ChIP-seq_mapping_summary.txt ./
cd ..
rm -r $basename

echo "=== All done successfully."
