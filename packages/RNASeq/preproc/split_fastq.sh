# split fastq files
#!/bin/bash

# 97-99

for CELL in {97..99} do
gunzip S00$CELL.fastq.gz
grep -A 3 \ 1:N S00$CELL.fastq | sed '/^--$/d' > /data/yosef2/users/xiuwei/fastq/S00$CELL_1.fastq
grep -A 3 \ 2:N S00$CELL.fastq | sed '/^--$/d' > /data/yosef2/users/xiuwei/fastq/S00$CELL_2.fastq
let CELL=CELL+1 
done

for CELL in {100..384} do
gunzip S0$CELL.fastq.gz
grep -A 3 \ 1:N S0$CELL.fastq | sed '/^--$/d' > /data/yosef2/users/xiuwei/fastq/S0$CELL_1.fastq
grep -A 3 \ 2:N S0$CELL.fastq | sed '/^--$/d' > /data/yosef2/users/xiuwei/fastq/S0$CELL_2.fastq
let CELL=CELL+1 
done
         


