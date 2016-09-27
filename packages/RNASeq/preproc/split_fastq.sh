# split fastq files
grep -A 3 \ 1:N S0097.fastq | sed '/^--$/d' > S0097_1.fastq
grep -A 3 \ 2:N S0097.fastq | sed '/^--$/d' > S0097_2.fastq
