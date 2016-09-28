# split fastq files
#!/bin/bash

# 97-99

for CELL in {97..99} 
do
gunzip S00${CELL}.fastq.gz
grep -A 3 \ 1:N S00${CELL}.fastq | sed '/^--$/d' > /data/yosef2/users/xiuwei/fastq/S00${CELL}_1_R1_combined.fastq
grep -A 3 \ 2:N S00${CELL}.fastq | sed '/^--$/d' > /data/yosef2/users/xiuwei/fastq/S00${CELL}_2_R1_combined.fastq
done         

for CELL in {100..384} 
do
gunzip S0${CELL}.fastq.gz
grep -A 3 \ 1:N S0${CELL}.fastq | sed '/^--$/d' > /data/yosef2/users/xiuwei/fastq/S0${CELL}_1_R1_combined.fastq
grep -A 3 \ 2:N S0${CELL}.fastq | sed '/^--$/d' > /data/yosef2/users/xiuwei/fastq/S0${CELL}_2_R1_combined.fastq
done         

for CELL in {97..99} 
do
mkdir S00${CELL}
mv S00${CELL}_1.fastq S00${CELL}/S00${CELL}_1_R1_combined.fastq
mv S00${CELL}_2.fastq S00${CELL}/S00${CELL}_2_R1_combined.fastq
done

for CELL in {100..384} 
do
mkdir S0${CELL}
mv S0${CELL}_1.fastq S0${CELL}/S0${CELL}_1_R1_combined.fastq
mv S0${CELL}_2.fastq S0${CELL}/S0${CELL}_2_R1_combined.fastq
done

for CELL in {98..99} 
do
cd S00${CELL}
mv S00${CELL}_1_R1_combined.fastq S00${CELL}_R1_combined.fastq
gzip S00${CELL}_R1_combined.fastq
mv S00${CELL}_2_R1_combined.fastq S00${CELL}_R2_combined.fastq
gzip S00${CELL}_R2_combined.fastq
cd ..
done

for CELL in {100..384} 
do
cd S0${CELL}
mv S0${CELL}_1_R1_combined.fastq S0${CELL}_R1_combined.fastq
gzip S0${CELL}_R1_combined.fastq
mv S0${CELL}_2_R1_combined.fastq S0${CELL}_R2_combined.fastq
gzip S0${CELL}_R2_combined.fastq
cd ..
done

