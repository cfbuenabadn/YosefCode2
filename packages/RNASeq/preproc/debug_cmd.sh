#!/bin/sh

Brain June 2015 cleanslate:
python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/processFolder.py -N JuneBrain --rsem_bowtie_maxins 1000 --kallisto_fragment_length 540 -p 1 -r mm10 -o /data/yosef/BRAIN/processed_June2015/150521_HS3A/ /data/yosef/BRAIN/sources/150521_HS3A/1_mismatch_better/Project_Ngai
python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/collectRsem.py --skip_collecting_dup_genes -r mm10 -o /home/eecs/allonwag/data/BRAIN/processed_June2015/rsem "/home/eecs/allonwag/data/BRAIN/processed_15_05/150515_HS3A;/home/eecs/allonwag/data/BRAIN/processed_15_05/150521_HS3A/;"


perl /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/collect_dat_cufflinks.pl "/home/eecs/allonwag/data/BRAIN/processed_15_05/150515_HS3A;/home/eecs/allonwag/data/BRAIN/processed_15_05/150521_HS3A/;" /home/eecs/allonwag/data/BRAIN/processed_June2015/cuff /data/yosef/index_files/mm10_4brain/index/rsem_index/rsemDictionary/mm10_4brain_rsemGeneMapping.txt 0

https://www.broadinstitute.org/~picard/picard_metric_definitions.html

python ~/project/singleCell/allon_script/preproc/normalizeBrainFastaInput.py /home/eecs/allonwag/data/BRAIN/sources/150309_HS1A
python ~/project/singleCell/allon_script/preproc/processFolder.py -N tfh_1st --paired_end --kallisto_fragment_length 300 -p 2 -r mm10 -o /data/yosef/TFH/processed/FC_01481 /data/yosef/TFH/sources/FC_01481
python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/collectRsem.py  --skip_collecting_dup_genes -r mm10 -o /home/eecs/allonwag/data/TFH/processed/FC_01481/rsem /home/eecs/allonwag/data/TFH/processed/FC_01481/
perl /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/collect_dat_cufflinks.pl /home/eecs/allonwag/data/TFH/processed/FC_01481/ /data/yosef/index_files/mm10_4brain/index/rsem_index/rsemDictionary/mm10_4brain_rsemGeneMapping.txt 0

Brain May 2015 runs:
python ~/project/singleCell/allon_script/preproc/processFolder.py -N May21Brain --rsem_bowtie_maxins 1000 --kallisto_fragment_length 540 -p 1 -r mm10 -o /data/yosef/BRAIN/processed_15_05/150521_HS3A/ /data/yosef/BRAIN/sources/150521_HS3A/1_mismatch_better/Project_Ngai
python ~/project/singleCell/allon_script/preproc/processFolder.py -N May21Olfa --rsem_bowtie_maxins 1000 --kallisto_fragment_length 540 -p 1 -r mm10 -o /data/yosef/BRAIN/processed_15_05/150515_HS3A /data/yosef/BRAIN/sources/150515_HS3A/Project_Ngai8by8

python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/collectRsem.py  --skip_collecting_dup_genes -r mm10 /data/yosef/BRAIN/processed_15_05/150521_HS3A/
python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/collectRsem.py  --skip_collecting_dup_genes -r mm10 /data/yosef/BRAIN/processed_15_05/150515_HS3A/
perl /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/collect_dat_cufflinks.pl /data/yosef/BRAIN/processed_15_05/150521_HS3A/ /data/yosef/index_files/mm10_4brain/index/rsem_index/rsemDictionary/mm10_4brain_rsemGeneMapping.txt 0
perl /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/collect_dat_cufflinks.pl /data/yosef/BRAIN/processed_15_05/150515_HS3A/ /data/yosef/index_files/mm10_4brain/index/rsem_index/rsemDictionary/mm10_4brain_rsemGeneMapping.txt 0


debug dup:
python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/processSingleSample.py --skip_tophat --skip_rsem --skip_kallisto --skip_trimmomatic -p 30 -r mm10 -o ~/archive/TFH/script/output_single ~/archive/TFH/script/gbc.fastq.gz

python ~/project/singleCell/allon_script/preproc/processSingleSample.py --kallisto_bootstrap_samples 100 --skip_tophat --skip_tophat_qc --skip_rsem --skip_rsem_qc --skip_qc --skip_trimmomatic --do_not_rely_on_previous_trimmomatic -p 16 -r mm10 -o ~/archive/TFH/script/output_single  ~/archive/TFH/script/gbc.fastq.gz

#python ~/project/singleCell/allon_script/preproc/processSingleSample.py --paired_end --kallisto_bootstrap_samples 100 --skip_tophat --skip_tophat_qc --skip_rsem --skip_rsem_qc --skip_qc --skip_trimmomatic --do_not_rely_on_previous_trimmomatic -p 16 -r mm10 -o ~/archive/TFH/script/output_paired  ~/archive/TFH/script/tfh1.fastq.gz ~/archive/TFH/script/tfh2.fastq.gz

echo "" > /var/mail/allonwag; rm -f *zen*.OU; rm -f *zen*.ER; python ~/project/singleCell/allon_script/preproc/processFolder.py -N tryRnaPipe -r mm10 -o ~/archive/TFH/pipeout_single ~/data/BRAIN/sources/olfactory/truncated_GBC_L01


python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/processSingleSample.py --do_not_clean_intermediary_files --paired_end --skip_kallisto --skip_rsem --skip_qc -p 20 -r mm10 -o ~/archive/TFH/script/output_paired /data/yosef/BRAIN/sources/150202_HS2A/Project_Ngai/Sample_OEP02_N712_S504/OEP02_N712_S504_GTAGAGGA-AGAGTAGA_L002_R1_combined.fastq.gz /data/yosef/BRAIN/sources/150202_HS2A/Project_Ngai/Sample_OEP02_N712_S504/OEP02_N712_S504_GTAGAGGA-AGAGTAGA_L002_R2_combined.fastq.gz
python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/processSingleSample.py --do_not_clean_intermediary_files --skip_kallisto --skip_rsem --skip_qc -p 20 -r mm10 -o ~/archive/TFH/script/output_single /data/yosef/BRAIN/sources/150202_HS2A/Project_Ngai/Sample_OEP02_N711_S502/OEP02_N711_S502_AAGAGGCA-CTCTCTAT_L002_R1_combined.fastq.gz


on queue yosef:
python ~/project/singleCell/allon_script/preproc/processFolder.py -N tryGBC_L01 -r mm10 -o ~/archive/TFH/pipeout_GBC_L01 ~/data/BRAIN/sources/olfactory/GBC_L01

on queue genomics:
python ~/project/singleCell/allon_script/preproc/processFolder.py -N tryGBC_L01 -r mm10 -o ~/archive/TFH/pipeout_P02-P03 ~/data/BRAIN/sources/olfactory/GBC_P02-P03

for jim (only rsem, on yosef_test) 
first batch:
python ~/project/singleCell/allon_script/preproc/processFolder.py --paired_end -N try4Jim -r mm10 -p 2 -o /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/MiSeq_PE25/LCMV_Day7_Plate1/processed /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/MiSeq_PE25/LCMV_Day7_Plate1/raw_data
python ~/project/singleCell/allon_script/preproc/processFolder.py --paired_end -N try4Jim -r mm10 -p 2 -o /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/MiSeq_PE25/LCMV_Day30_Plate1/processed /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/MiSeq_PE25/LCMV_Day30_Plate1/raw_data
second batch:
python ~/project/singleCell/allon_script/preproc/processFolder.py -N try4JimB --skip_tophat --skip_tophat_qc -r mm10 -p 2 -o /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/LCMV_Pilot/LCMV_Plates1_d7_Arm/processed2 /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/LCMV_Pilot/LCMV_Plates1_d7_Arm/raw_data
python ~/project/singleCell/allon_script/preproc/processFolder.py -N try4JimB --skip_tophat --skip_tophat_qc -r mm10 -p 2 -o /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/LCMV_Pilot/LCMV_Plates2_d7_Arm/processed2 /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/LCMV_Pilot/LCMV_Plates2_d7_Arm/raw_data
python ~/project/singleCell/allon_script/preproc/processFolder.py -N try4JimB --skip_tophat --skip_tophat_qc -r mm10 -p 2 -o /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/LCMV_Pilot/LCMV_Plates3_d30_Clone/processed2 /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/LCMV_Pilot/LCMV_Plates3_d30_Clone/raw_data
python ~/project/singleCell/allon_script/preproc/processFolder.py -N try4JimB --skip_tophat --skip_tophat_qc -r mm10 -p 2 -o /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/LCMV_Pilot/LCMV_Plates4_d7_Arm/processed2 /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/LCMV_Pilot/LCMV_Plates4_d7_Arm/raw_data

python ~/project/singleCell/allon_script/preproc/processSingleSample.py --paired_end --skip_trimmomatic --skip_tophat --skip_rsem --rsem_bowtie_maxins 3000 -r mm10 -p 4 -o /home/eecs/allonwag/data/BRAIN/processed_debug/OEP02_N712_S508_GTAGAGGA-CTAAGCCT_L002_R1_combined /home/eecs/allonwag/data/BRAIN/sources/150202_HS2A/Project_Ngai/Sample_OEP02_N712_S508/OEP02_N712_S508_GTAGAGGA-CTAAGCCT_L002_R1_combined.fastq.gz /home/eecs/allonwag/data/BRAIN/sources/150202_HS2A/Project_Ngai/Sample_OEP02_N712_S508/OEP02_N712_S508_GTAGAGGA-CTAAGCCT_L002_R2_combined.fastq.gz
python ~/project/singleCell/allon_script/preproc/processSingleSample.py --skip_trimmomatic --skip_tophat --skip_rsem --rsem_bowtie_maxins 3000 -r mm10 -p 4 -o /home/eecs/allonwag/data/BRAIN/processed_debug_single/OEP02_N712_S508_GTAGAGGA-CTAAGCCT_L002_R1_combined /home/eecs/allonwag/data/BRAIN/sources/150202_HS2A/Project_Ngai/Sample_OEP02_N712_S508/OEP02_N712_S508_GTAGAGGA-CTAAGCCT_L002_R1_combined.fastq.gz


python ~/project/singleCell/allon_script/preproc/processFolder.py --skip_tophat --skip_tophat_qc -N brainP02P03 -r mm10 -p 2 -o /home/eecs/allonwag/data/BRAIN/processed2/GBC_P02-P03 /home/eecs/allonwag/data/BRAIN/sources/olfactory/GBC_P02-P03
python ~/project/singleCell/allon_script/preproc/processFolder.py --skip_tophat --skip_tophat_qc -N brainL01 -r mm10 -p 2 -o /home/eecs/allonwag/data/BRAIN/processed2/GBC_L01 /home/eecs/allonwag/data/BRAIN/sources/olfactory/GBC_L01
python ~/project/singleCell/allon_script/preproc/processFolder.py --skip_tophat --skip_tophat_qc -N brainL02 -r mm10 -p 2 -o /home/eecs/allonwag/data/BRAIN/processed2/GBC_L02 /home/eecs/allonwag/data/BRAIN/sources/olfactory/GBC_L02
python ~/project/singleCell/allon_script/preproc/processFolder.py --skip_trimmomatic --skip_tophat --skip_rsem --rsem_bowtie_maxins 3000 --paired_end -N brainPaired -r mm10 -p 1 -o /home/eecs/allonwag/data/BRAIN/processed3/150202_HS2A/Project_Ngai /home/eecs/allonwag/data/BRAIN/sources/150202_HS2A/Project_Ngai
python ~/project/singleCell/allon_script/preproc/processFolder.py --skip_trimmomatic --skip_tophat --skip_rsem --rsem_bowtie_maxins 3000 -N brainPairedAsSingle -r mm10 -p 1 -o /home/eecs/allonwag/data/BRAIN/processed3/150202_HS2A/Project_Ngai_AsSingle /home/eecs/allonwag/data/BRAIN/sources/150202_HS2A/Project_Ngai
python ~/project/singleCell/allon_script/preproc/processFolder.py --skip_tophat --skip_tophat_qc -N brainPairedAsSingleOtherEnd -r mm10 -p 2 -o /home/eecs/allonwag/data/BRAIN/processed2/150202_HS2A/Project_Ngai_AsSingleOtherEnd /home/eecs/allonwag/data/BRAIN/sources/150202_HS2A/Project_Ngai
python ~/project/singleCell/allon_script/preproc/processFolder.py --skip_trimmomatic --skip_tophat --skip_rsem --rsem_bowtie_maxins 3000 --paired_end -N brain96_192 -r mm10 -p 4 -o /home/eecs/allonwag/data/BRAIN/processed3/2015_3_13/150309_HS3A/Project_Ngai /home/eecs/allonwag/data/BRAIN/sources/2015_3_13/150309_HS3A/Project_Ngai 
python ~/project/singleCell/allon_script/preproc/processFolder.py --rsem_bowtie_maxins 3000 -N firstDavidBatch -r mm10 -p 2 -o /data/yosef/BRAIN/processed2/150309_HS1A /data/yosef/BRAIN/sources/150309_HS1A


python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/processFolder.py --skip_trimmomatic --skip_tophat --skip_rsem --trimmomatic_window 4:20 --rsem_bowtie_maxins 3000 --paired_end -N brainPaired -r mm10 -p 1 -o /home/eecs/allonwag/data/BRAIN/processed4/150202_HS2A/Project_Ngai /home/eecs/allonwag/data/BRAIN/sources/150202_HS2A/Project_Ngai
python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/processFolder.py --skip_trimmomatic --skip_tophat --skip_rsem --trimmomatic_window 4:20 --rsem_bowtie_maxins 3000 -N brainPairedAsSingle -r mm10 -p 1 -o /home/eecs/allonwag/data/BRAIN/processed4/150202_HS2A/Project_Ngai_AsSingle /home/eecs/allonwag/data/BRAIN/sources/150202_HS2A/Project_Ngai
python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/collectRsem.py -r mm10 -o /home/eecs/allonwag/data/BRAIN/processed4/150202_HS2A/Project_Ngai/rsem  /home/eecs/allonwag/data/BRAIN/processed4/150202_HS2A/Project_Ngai
python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/collectRsem.py  -r mm10 -o /home/eecs/allonwag/data/BRAIN/processed4/150202_HS2A/Project_Ngai_AsSingle/rsem /home/eecs/allonwag/data/BRAIN/processed4/150202_HS2A/Project_Ngai_AsSingle
perl /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/collect_dat_cufflinks.pl /home/eecs/allonwag/data/BRAIN/processed4/150202_HS2A/Project_Ngai /data/yosef/index_files/mm10_4brain/index/rsem_index/rsemDictionary/mm10_4brain_rsemGeneMapping.txt 0
perl /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/collect_dat_cufflinks.pl /home/eecs/allonwag/data/BRAIN/processed4/150202_HS2A/Project_Ngai_AsSingle /data/yosef/index_files/mm10_4brain/index/rsem_index/rsemDictionary/mm10_4brain_rsemGeneMapping.txt 0
python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/processFolder.py --trimmomatic_window 4:30 --rsem_bowtie_maxins 3000 --paired_end -N brainPaired -r mm10 -p 1 -o /home/eecs/allonwag/data/BRAIN/processed5/150202_HS2A/Project_Ngai /home/eecs/allonwag/data/BRAIN/sources/150202_HS2A/Project_Ngai
python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/processFolder.py --trimmomatic_window 4:30 --rsem_bowtie_maxins 3000 -N brainPairedAsSingle -r mm10 -p 1 -o /home/eecs/allonwag/data/BRAIN/processed5/150202_HS2A/Project_Ngai_AsSingle /home/eecs/allonwag/data/BRAIN/sources/150202_HS2A/Project_Ngai



python ~/project/singleCell/allon_script/preproc/collectRsem.py  -r mm10 -o /home/eecs/allonwag/data/BRAIN/processed2/GBC_L01/rsem /home/eecs/allonwag/data/BRAIN/processed2/GBC_L01
python ~/project/singleCell/allon_script/preproc/collectRsem.py  -r mm10 -o /home/eecs/allonwag/data/BRAIN/processed2/GBC_L02/rsem /home/eecs/allonwag/data/BRAIN/processed2/GBC_L02
python ~/project/singleCell/allon_script/preproc/collectRsem.py  -r mm10 -o /home/eecs/allonwag/data/BRAIN/processed2/GBC_P02-P03/rsem /home/eecs/allonwag/data/BRAIN/processed2/GBC_P02-P03
python ~/project/singleCell/allon_script/preproc/collectRsem.py  -r mm10 -o /home/eecs/allonwag/data/BRAIN/processed3/150202_HS2A/Project_Ngai/rsem /home/eecs/allonwag/data/BRAIN/processed3/150202_HS2A/Project_Ngai
python ~/project/singleCell/allon_script/preproc/collectRsem.py  -r mm10 -o /home/eecs/allonwag/data/BRAIN/processed2/150202_HS2A/Project_Ngai_AsSingleOtherEnd/rsem /home/eecs/allonwag/data/BRAIN/processed2/150202_HS2A/Project_Ngai_AsSingleOtherEnd
python ~/project/singleCell/allon_script/preproc/collectRsem.py  -r mm10 -o /home/eecs/allonwag/data/BRAIN/processed3/150202_HS2A/Project_Ngai_AsSingle/rsem /home/eecs/allonwag/data/BRAIN/processed3/150202_HS2A/Project_Ngai_AsSingle
python ~/project/singleCell/allon_script/preproc/collectRsem.py -r mm10 -o /home/eecs/allonwag/data/BRAIN/processed3/2015_3_13/150309_HS3A/Project_Ngai/rsem  /home/eecs/allonwag/data/BRAIN/processed3/2015_3_13/150309_HS3A/Project_Ngai
python ~/project/singleCell/allon_script/preproc/collectRsem.py -r mm10 -o /data/yosef/BRAIN/processed2/150309_HS1A/rsem  /data/yosef/BRAIN/processed2/150309_HS1A


%shows that the parameter that controls the insert size works, but it simply does not increase above 200
python ~/project/singleCell/allon_script/preproc/processSingleSample.py --paired_end --rsem_bowtie_maxins 100 -p 8 -r mm10 --skip_tophat --skip_tophat_qc -o /home/eecs/allonwag/data/BRAIN/processed_debug/rsemMaxIns100 ~/data/BRAIN/sources_debug/100/R1.fastq.gz   ~/data/BRAIN/sources_debug/100/R2.fastq.gz
python ~/project/singleCell/allon_script/preproc/processSingleSample.py --paired_end --rsem_bowtie_maxins 200 -p 8 -r mm10 --skip_tophat --skip_tophat_qc -o /home/eecs/allonwag/data/BRAIN/processed_debug/rsemMaxIns200 ~/data/BRAIN/sources_debug/200/R1.fastq.gz   ~/data/BRAIN/sources_debug/200/R2.fastq.gz
python ~/project/singleCell/allon_script/preproc/processSingleSample.py --paired_end --rsem_bowtie_maxins 400 -p 8 -r mm10 --skip_tophat --skip_tophat_qc -o /home/eecs/allonwag/data/BRAIN/processed_debug/rsemMaxIns400 ~/data/BRAIN/sources_debug/400/R1.fastq.gz   ~/data/BRAIN/sources_debug/400/R2.fastq.gz
python ~/project/singleCell/allon_script/preproc/processSingleSample.py --paired_end --rsem_bowtie_maxins 600 -p 8 -r mm10 --skip_tophat --skip_tophat_qc -o /home/eecs/allonwag/data/BRAIN/processed_debug/rsemMaxIns600 ~/data/BRAIN/sources_debug/600/R1.fastq.gz   ~/data/BRAIN/sources_debug/600/R2.fastq.gz
python ~/project/singleCell/allon_script/preproc/processSingleSample.py --paired_end --rsem_bowtie_maxins 750 -p 8 -r mm10 --skip_tophat --skip_tophat_qc -o /home/eecs/allonwag/data/BRAIN/processed_debug/rsemMaxIns750 ~/data/BRAIN/sources_debug/750/R1.fastq.gz   ~/data/BRAIN/sources_debug/750/R2.fastq.gz
python ~/project/singleCell/allon_script/preproc/processSingleSample.py --paired_end --rsem_bowtie_maxins 1000 -p 8 -r mm10 --skip_tophat --skip_tophat_qc -o /home/eecs/allonwag/data/BRAIN/processed_debug/rsemMaxIns1000 ~/data/BRAIN/sources_debug/1000/R1.fastq.gz ~/data/BRAIN/sources_debug/1000/R2.fastq.gz
python ~/project/singleCell/allon_script/preproc/processSingleSample.py --paired_end --rsem_bowtie_maxins 1500 -p 8 -r mm10 --skip_tophat --skip_tophat_qc -o /home/eecs/allonwag/data/BRAIN/processed_debug/rsemMaxIns1500 ~/data/BRAIN/sources_debug/1500/R1.fastq.gz ~/data/BRAIN/sources_debug/1500/R2.fastq.gz


/150202_HS2A/Project_Ngai

python ~/project/singleCell/allon_script/preproc/tempCollect.py /home/eecs/allonwag/data/BRAIN/processed3/150202_HS2A/Project_Ngai_AsSingle


DEBUG RSEM OUTPUT - REMOVE WEIGHTS
/opt/pkg/rsem-1.2.19/bin/rsem-calculate-expression --sampling-for-bam --num-threads 16 --bowtie2 --estimate-rspd --output-genome-bam --paired-end /home/eecs/allonwag/data/BRAIN/processed/150202_HS2A/Project_Ngai/OEP01_N710_S501_CGAGGCTG-TAGATCGC_L001_R1_combined/trimmomatic_output/1.Ptrim.fastq /home/eecs/allonwag/data/BRAIN/processed/150202_HS2A/Project_Ngai/OEP01_N710_S501_CGAGGCTG-TAGATCGC_L001_R1_combined/trimmomatic_output/2.Ptrim.fastq /data/yosef/index_files/mm10_4brain/index/rsem_index/GRCm38.p3_4brain_rsem /home/eecs/allonwag/data/BRAIN/processed/150202_HS2A/Project_Ngai/OEP01_N710_S501_CGAGGCTG-TAGATCGC_L001_R1_combined/bla/rsem_output
/opt/pkg/rsem-1.2.19/bin/rsem-calculate-expression --sampling-for-bam --num-threads 16 --bowtie2 --estimate-rspd --output-genome-bam /home/eecs/allonwag/data/BRAIN/processed/150202_HS2A/Project_Ngai_AsSingle/OEP01_N710_S501_CGAGGCTG-TAGATCGC_L001_R1_combined/trimmomatic_output/1.Utrim.fastq /data/yosef/index_files/mm10_4brain/index/rsem_index/GRCm38.p3_4brain_rsem /home/eecs/allonwag/data/BRAIN/processed/150202_HS2A/Project_Ngai_AsSingle/OEP01_N710_S501_CGAGGCTG-TAGATCGC_L001_R1_combined/bla/rsem_output
python ~/project/singleCell/allon_script/preproc/processFolder.py \
--paired_end -r mm10 -p 16 \
--skip_tophat_qc --skip_tophat \
-o ~/data/BRAIN/sources/150202_HS2A/forDebugProcessed \
~/data/BRAIN/sources/150202_HS2A/forDebug

python ~/project/singleCell/allon_script/preproc/processFolder.py \
-r mm10 -p 16 \
--skip_tophat_qc --skip_tophat \
-o ~/data/BRAIN/sources/150202_HS2A/forDebugProcessedAsSingle \
~/data/BRAIN/sources/150202_HS2A/forDebug

rsem keeps multiple reads for all possibilities unless you tell him to sample the file
How to count mapped reads: https://www.biostars.org/p/69109/
samtools view -S -c -q 100 accepted_hits.sam

FIX DAVIDE PROBLEM - RSEM abuse of number of reads - multiple alignments by sampling
Nir - updated count reads script
mark duplicates - in single end? in paired end?
Read - single end, fragment - both ends of the pair


For Michael:
python ~/project/singleCell/allon_script/preproc/processSingleSample.py --skip_trimmomatic --paired_end --skip_tophat --skip_tophat_qc -r hg38 -o ~/archive/TFH/script/output_michael /data/yosef/HIV/dat/T_Cell_Runs10-11/test_batch_for_Allon/singlecell_test_sample/small_R1_combined.fastq.gz /data/yosef/HIV/dat/T_Cell_Runs10-11/test_batch_for_Allon/singlecell_test_sample/small_R2_combined.fastq.gz

python ~/project/singleCell/allon_script/preproc/processFolder.py -N tryMichael -r hg38 -p 2 --paired_end --skip_tophat --skip_tophat_qc -o ~/archive/TFH/pipeout_4Michael /data/yosef/HIV/dat/T_Cell_Runs10-11/data



brain concatenation - done
brain output rsem - multiple files
brain - rsem collect
ribosomal file

example for split sample:
~/data/BRAIN/sources/olfactory/GBC_P02-P03/Sample_GBCP02_N707_S506

DEBUG error:
GBCP02_N709_S506_GCTACGCT-ACTGCATA_L007
GBCP02_N707_S506_CTCTCTAC-ACTGCATA_L007





echo "" > /var/mail/allonwag; rm -rf ~/archive/TFH/pipeout_single/RFSN_N708; rm -rf ~/archive/TFH/pipeout_single/RFSN_N709; rm -f *zen*.OU; rm -f *zen*.ER; python ~/project/singleCell/allon_script/preproc/processFolder.py -N tryRnaPipe -r mm10 -o ~/archive/TFH/pipeout_single ~/data/BRAIN/sources/olfactory/truncated_GBC_L01


#In the PicardOutput folder:
#BamIndexStats INPUT=sorted.bam
#samtools view -h sorted.bam > sorted.sam


perl check_contam.pl /home/eecs/allonwag/archive/TFH/script/gbc.fastq -o  /home/eecs/allonwag/archive/TFH/script/output_single/fastqc_output/primer.1.txt






/opt/genomics/bin/bowtie2  ~/data/index_files/mm10_4brain/index/rsem_index/GRCm38.p3_4brain_rsem \
--sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 \
--no-mixed --no-discordant \
-p 6 -1 gbc.fastq -S ./rsem_output/aligned_by_bowtie2.sam


OLD:
#rsemComand = Template("/opt/pkg/rsem-1.2.19/bin/rsem-calculate-expression --sam-header-info samHeader.txt --bowtie2 --estimate-rspd --paired-end --output-genome-bam $SAMPLE_FILE1 $SAMPLE_FILE2 $RSEM_INDEX $OUTPUT_FOLDER/rsem_output/rsem_output").substitute(OUTPUT_FOLDER=args.output_folder, RSEM_INDEX=RSEM_INDEX, SAMPLE_FILE1=args.sampleFile1, SAMPLE_FILE2=args.sampleFile2);

#	print("Running alignment for rsem");
#	#I discovered that I have to run the rsem alignment separately: if I let rsem align, and then use the option that outputs a genomic  bam (in addition to the usual transcript bam) - then it does not include the ERCC spikeins and other transcripts that I added - it includes only genes that appear in the gtf file...
#	#so for the purpose of QC I run the alignment separately
#	#The bowtie2 command line args are the ones that are used by rsem by default when it uses bowtie as its internal aligner
#	if(args.paired_end):
#		rsemAlignmentComand = Template("/opt/genomics/bin/bowtie2 $RSEM_INDEX \
#						--sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 \
#						--no-mixed --no-discordant -p 6 \
#						-1 $SAMPLE_FILE1 -2 $SAMPLE_FILE2 -S $OUTPUT_FOLDER/rsem_output/aligned_4rsem.sam").substitute( \
#						OUTPUT_FOLDER=args.output_folder, RSEM_INDEX=RSEM_INDEX, SAMPLE_FILE1=args.sampleFile1, SAMPLE_FILE2=args.sampleFile2);
#	else:
#		rsemAlignmentComand = Template("/opt/genomics/bin/bowtie2 $RSEM_INDEX \
#						--sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 \
#						--no-mixed --no-discordant -p 6 \
#						-1 $SAMPLE_FILE1 -2 $SAMPLE_FILE2 -S $OUTPUT_FOLDER/rsem_output/aligned_4rsem.sam").substitute( \
#						OUTPUT_FOLDER=args.output_folder, RSEM_INDEX=RSEM_INDEX, SAMPLE_FILE1=args.sampleFile1, SAMPLE_FILE2=args.sampleFile2);
#
#		rsemAlignmentComand = Template("/opt/pkg/rsem-1.2.19/bin/rsem-calculate-expression --bowtie2 --estimate-rspd --output-genome-bam $SAMPLE_FILE1 $RSEM_INDEX $OUTPUT_FOLDER/rsem_output/rsem_output").substitute(OUTPUT_FOLDER=args.output_folder, RSEM_INDEX=RSEM_INDEX, SAMPLE_FILE1=args.sampleFile1);
#		
#	print(rsemAlignmentComand)
#	returnCode = subprocess.call(rsemAlignmentComand, shell=True);
#	if(returnCode != 0):
#		raise Exception("rsem failed");