from string import Template
import os
import stat

#TODO: change process folder to detect automatically if this is single end or paired-end?

sourceFolders = [("olfactory/GBC_L01", False),
                 ("olfactory/GBC_L02", False),
                 ("olfactory/GBC_P02-P03", False),
                 ("olfactory/HBC", False),
                 ("150202_HS2A/Project_Ngai", True),
                 ("150309_HS1A", False),
                 ("150309_HS3A/Project_Ngai_new", True),
                 ("150515_HS3A/Project_Ngai8by8", False),
                 ("150521_HS3A/1_mismatch_better/Project_Ngai", False)
]

SEND_SCRIPT_FILE_NAME = "/project/eecs/yosef-archive/users/allonwag/temp/sendAllBrain.sh"
COLLECT_SCRIPT_FILE_NAME = "/project/eecs/yosef-archive/users/allonwag/temp/collectAllBrain.sh"

allBatchAliases = []
with open(SEND_SCRIPT_FILE_NAME, "wt") as fout:
    fout.write("#!/bin/sh" + '\n\n\n')

    for (batchDir, isPairedEnd) in sourceFolders:
        batchAlias = batchDir.replace("/", "-") #make a hierarchy of one batch folder and beneath it all the cell folders
        send_cmd = Template("python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/processFolder.py -N JuneBrain_b" +
                            (" --paired_end" if isPairedEnd else "") +
                            " --rsem_bowtie_maxins 1000 --kallisto_fragment_length 540 -p 1 -r mm10" +
                            " -o $OUT_FOLDER $INPUT_FOLDER").substitute(
                            OUT_FOLDER=os.path.join("/data/yosef/BRAIN/processed_June2015_b/", batchAlias),
                            INPUT_FOLDER=os.path.join("/data/yosef/BRAIN/sources/", batchDir)
                            )
        fout.write(send_cmd + '\n\n')

        allBatchAliases.append(batchAlias)

with open(COLLECT_SCRIPT_FILE_NAME, "wt") as fout:
    fout.write("#!/bin/sh" + '\n\n\n')

    collect_cmd = Template("python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/collectRsem.py --skip_collecting_dup_genes" +
                           " -r mm10 -o /home/eecs/allonwag/data/BRAIN/processed_June2015/rsem" +
                            " $DIRECTORORIES_TO_PROCESS").substitute(DIRECTORORIES_TO_PROCESS=';'.join(allBatchAliases))

    fout.write(collect_cmd + '\n\n')



for filename in [SEND_SCRIPT_FILE_NAME, COLLECT_SCRIPT_FILE_NAME]:
    #add chmod u+x (add user run permissions) to the files
    st = os.stat(filename)
    os.chmod(filename, st.st_mode | stat.S_IEXEC)

