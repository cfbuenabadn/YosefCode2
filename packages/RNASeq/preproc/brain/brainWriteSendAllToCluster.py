from string import Template
import os
import stat
from brainSources import cortical_sourceFolders, olfactory_sourceFolders, bateup_sourceFolders, samIsrael_sourceFolders

#TODO: change process folder to detect automatically if this is single end or paired-end?

SEND_SCRIPT_FILE_NAME = "/project/eecs/yosef-archive/users/allonwag/temp/sendAllBrain.sh"
COLLECT_SCRIPT_FILE_NAME = "/project/eecs/yosef-archive/users/allonwag/temp/collectAllBrain.sh"

#OUTPUT_FOLDER = "/data/yosef/BRAIN/processed_June2015_b/"
#OUTPUT_FOLDER = "/data/yosef/BRAIN/processed_Sep2015/"
#OUTPUT_FOLDER = "/data/yosef/BRAIN/processed_July2015/"
#OUTPUT_FOLDER = "/data/yosef/BRAIN/processed_Bateup_Aug2015/"
#OUTPUT_FOLDER = "/data/yosef/BRAIN/processed_Zebrafish_Oct2015/"
#OUTPUT_FOLDER = "/data/yosef2/BRAIN/processed_olfactory_Jun2016/"
OUTPUT_FOLDER = "/data/yosef2/BRAIN/processed_cortical_Oct2016/"

#Cortical / Olfactory / Bateup / SamIsrael
PROJECT = "Cortical"
if PROJECT == "Cortical":
    sourceFolders = cortical_sourceFolders
    reference_genome = "mm10"
elif PROJECT == "Olfactory":
    sourceFolders = olfactory_sourceFolders
    reference_genome = "mm10"
elif PROJECT == "Bateup":
    sourceFolders = bateup_sourceFolders
    reference_genome = "mm10"
elif PROJECT == "SamIsrael":
    sourceFolders = samIsrael_sourceFolders
    reference_genome = "z10"
else:
    raise Exception("unrecognized project!")



allBatchAliases = []
with open(SEND_SCRIPT_FILE_NAME, "wt") as fout:
    fout.write("#!/bin/sh" + '\n\n\n')

    for (batchDir, isPairedEnd) in sourceFolders:
        batchAlias = batchDir.replace("/", "-") #make a hierarchy of one batch folder and beneath it all the cell folders
        send_cmd = Template("python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/processFolder.py -N JuneBrain_b" +
                            (" --paired_end" if isPairedEnd else "") +
                            " --skip_rsem --skip_kallisto --skip_rsem_qc --skip_kallisto_qc --rsem_bowtie_maxins 1000 --mean_fragment_length 540 --std_fragment_length 250 -p 1 -r $REFERENCE_GENOME" +
                            " -o $OUT_FOLDER $INPUT_FOLDER").substitute(
                            OUT_FOLDER=os.path.join(OUTPUT_FOLDER, batchAlias),
                            INPUT_FOLDER=os.path.join("/data/yosef/BRAIN/sources/", batchDir),
                            REFERENCE_GENOME = reference_genome
                            )
        fout.write(send_cmd + '\n\n')

        allBatchAliases.append(batchAlias)

with open(COLLECT_SCRIPT_FILE_NAME, "wt") as fout:
    fout.write("#!/bin/sh" + '\n\n\n')

    fout.write("#Collect all batches together:\n\n")

    collect_cmd = Template("python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/collectPreprocResults.py --skip_collecting_dup_genes" +
                           " -r $REFERENCE_GENOME -o $OUTPUT_FOLDER" +
                            " \"$DIRECTORORIES_TO_PROCESS\"").substitute(OUTPUT_FOLDER=os.path.join(OUTPUT_FOLDER, "collect"), DIRECTORORIES_TO_PROCESS=';'.join([OUTPUT_FOLDER + alias for alias in allBatchAliases]), REFERENCE_GENOME = reference_genome)

    fout.write(collect_cmd + '\n\n')

    #no longer needed - collect cufflinks was merged into the collect rsem
    # collect_cuff_cmd = Template("perl /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/collect_dat_cufflinks.pl" +
    # " \"$DIRECTORORIES_TO_PROCESS\" $OUTPUT_FOLDER" +
    # " /data/yosef/index_files/mm10_4brain/index/rsem_index/rsemDictionary/mm10_4brain_rsemGeneMapping.txt 0").substitute(OUTPUT_FOLDER=OUTPUT_FOLDER, DIRECTORORIES_TO_PROCESS=';'.join([OUTPUT_FOLDER + alias for alias in allBatchAliases]))
    #
    # fout.write(collect_cuff_cmd + '\n\n')
    fout.write('###############################################################' + '\n')
    fout.write('###############################################################' + '\n\n')
    fout.write("#Collect each batch separately:\n\n")

    #write individual collect commands for each of the libraries
    for batchAlias in allBatchAliases:

            currentFolder = OUTPUT_FOLDER + batchAlias
            collect_cmd = Template("# python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/collectPreprocResults.py --skip_collecting_dup_genes" +
                           " -r $REFERENCE_GENOME -o $OUTPUT_FOLDER" +
                            " $INPUT_FOLDER").substitute(OUTPUT_FOLDER=os.path.join(currentFolder, "collect"), INPUT_FOLDER=currentFolder, REFERENCE_GENOME = reference_genome)

            fout.write(collect_cmd + '\n\n')

            #no longer needed - collect cufflinks was merged into the collect rsem
            # collect_cuff_cmd = Template("perl /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/collect_dat_cufflinks.pl" +
            #                             " $INPUT_FOLDER $OUTPUT_FOLDER" +
            #                             " /data/yosef/index_files/mm10_4brain/index/rsem_index/rsemDictionary/mm10_4brain_rsemGeneMapping.txt 0").\
            #                             substitute(OUTPUT_FOLDER=currentFolder, INPUT_FOLDER=currentFolder)
            #
            # fout.write(collect_cuff_cmd + '\n\n')



for filename in [SEND_SCRIPT_FILE_NAME, COLLECT_SCRIPT_FILE_NAME]:
    #add chmod u+x (add user run permissions) to the files
    st = os.stat(filename)
    os.chmod(filename, st.st_mode | stat.S_IEXEC)

