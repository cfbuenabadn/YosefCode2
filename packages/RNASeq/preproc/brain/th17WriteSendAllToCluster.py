from string import Template
import os
import stat
from brainSources import cortical_sourceFolders, olfactory_sourceFolders, bateup_sourceFolders, samIsrael_sourceFolders

#TODO: change process folder to detect automatically if this is single end or paired-end?

SEND_SCRIPT_FILE_NAME = "/project/eecs/yosef-archive/users/allonwag/temp/sendAllTh17.sh"
COLLECT_SCRIPT_FILE_NAME = "/project/eecs/yosef-archive/users/allonwag/temp/collectAllTh17.sh"

OUTPUT_FOLDER = "/data/yosef2/Published_Data/TH17/processed_aw20160630"

th17_sourceFolders = [
    ("fqs", True)
]


allBatchAliases = []
with open(SEND_SCRIPT_FILE_NAME, "wt") as fout:
    fout.write("#!/bin/sh" + '\n\n\n')

    for (batchDir, isPairedEnd) in th17_sourceFolders:
        batchAlias = batchDir.replace("/", "-") #make a hierarchy of one batch folder and beneath it all the cell folders
        send_cmd = Template("python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/processFolder.py -N Th17" +
                            (" --paired_end" if isPairedEnd else "") +
                            " --single_cell_prior -p 1 -r mm10" +
                            " -o $OUT_FOLDER $INPUT_FOLDER").substitute(
                            OUT_FOLDER=os.path.join(OUTPUT_FOLDER, batchAlias),
                            INPUT_FOLDER=os.path.join("/data/yosef2/Published_Data/TH17/", batchDir),
                            )
        fout.write(send_cmd + '\n\n')

        allBatchAliases.append(batchAlias)

with open(COLLECT_SCRIPT_FILE_NAME, "wt") as fout:
    fout.write("#!/bin/sh" + '\n\n\n')

    fout.write("#Collect all batches together:\n\n")

    collect_cmd = Template("python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/collectPreprocResults.py --skip_collecting_dup_genes" +
                           " -r mm10 -o $OUTPUT_FOLDER" +
                            " \"$DIRECTORORIES_TO_PROCESS\"").substitute(OUTPUT_FOLDER=os.path.join(OUTPUT_FOLDER, "collect"), DIRECTORORIES_TO_PROCESS=';'.join([OUTPUT_FOLDER + alias for alias in allBatchAliases]))

    fout.write(collect_cmd + '\n\n')


    fout.write('###############################################################' + '\n')
    fout.write('###############################################################' + '\n\n')
    fout.write("#Collect each batch separately:\n\n")

    #write individual collect commands for each of the libraries
    for batchAlias in allBatchAliases:

            currentFolder = OUTPUT_FOLDER + batchAlias
            collect_cmd = Template("# python /data/yosef/users/allonwag/YosefCode/packages/RNASeq/preproc/collectPreprocResults.py --skip_collecting_dup_genes" +
                           " -r mm10 -o $OUTPUT_FOLDER" +
                            " $INPUT_FOLDER").substitute(OUTPUT_FOLDER=os.path.join(currentFolder, "collect"), INPUT_FOLDER=currentFolder)




for filename in [SEND_SCRIPT_FILE_NAME, COLLECT_SCRIPT_FILE_NAME]:
    #add chmod u+x (add user run permissions) to the files
    st = os.stat(filename)
    os.chmod(filename, st.st_mode | stat.S_IEXEC)

