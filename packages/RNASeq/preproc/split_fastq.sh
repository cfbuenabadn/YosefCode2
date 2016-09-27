# split fastq files
#!/bin/bash

# 97-99
# 100-385
CELL=100
while [  $CELL -lt 385 ]; do
gunzip S0$CELL.fastq.gz grep -A 3 \ 1:N | sed '/^--$/d' > /data/yosef2/users/xiuwei/fastq/S0$CELL_1.fastq
gunzip S0$CELL.fastq.gz grep -A 3 \ 2:N | sed '/^--$/d' > /data/yosef2/users/xiuwei/fastq/S0$CELL_2.fastq
let CELL=CELL+1 
done
         


