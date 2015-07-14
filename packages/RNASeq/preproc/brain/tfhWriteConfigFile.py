from brainSources import sourceFolders
import os
import os.path
import glob
import openpyxl
import re
from collections import namedtuple

#folder name, isPairedEnd
sourceFolders = [("FC_01481", True)
]

#The sourceFolders is a tuple of dir and bool (isPairedEnd) - get rid of the latter field which is not needed here
sourceFolders = [x[0] for x in sourceFolders]

SOURCES_DIR = "/data/yosef/TFH/sources"
PROCESSED_DIR = "/data/yosef/TFH/processed"

CONFIG_OUTPUT_FILE = os.path.join(PROCESSED_DIR, 'collect/config_tfh.xlsx')

#observe that the cuff cell list and the rsem cell list (and in the future kallisto cell list?!) should all be the same
outputCellListFileName = os.path.join(PROCESSED_DIR, 'collect/cuff/cell_list.txt')


def ProcessOutputCellList(outputCellListFileName):
    cellID2outputName_dict = {}

    outputCellList = [row.strip() for row in open(outputCellListFileName).readlines()]

    #make the output cell list into a dictionary where the cells unique id points to its dir structure
    extractUniqueIDFromCellName = re.compile(r"^(?P<path>[\d\w_\-]+)/(?P<uniqueID>[\w_]+)_[ACGT]+\-[ACGT]+_L\d\d\d$")
    for cell in outputCellList:
        #print cell
        m = re.match(extractUniqueIDFromCellName, cell)

        if not(m.groupdict().has_key("path")) or not(m.groupdict().has_key("uniqueID")):
            raise Exception("Reading collect output list: Cannot parse cell name from the output list! (cell name: %s)" % cell)

        cell_uniqueID = m.group("uniqueID")
        if(cellID2outputName_dict.has_key(cell_uniqueID)):
            raise Exception("Reading collect output list: Same cell encountered twice in the output?! (cell name: %s)" % cell)

        cellID2outputName_dict[cell_uniqueID] = cell


    return cellID2outputName_dict


def WriteConfigFile(elementsToWrite):
    #the values in the element list are: (cellName, metadata_for_cell, outputName_for_cell)

    #use only the values, sort them by the cell's unique id which is the first element in the tuple
    elementsToWrite = sorted(elementsToWrite.values(), key=lambda x: x[0])

    wb = openpyxl.Workbook(write_only=True)
    ws = wb.create_sheet()

    #compose and print the header line - it is made of the header fields, I just want to place the unique ID first before all the rest of the fields
    #and prefix every MD on the cell with "MD_"
    headerLineColumns = ['unique_id', 'output_name']
    ws.append(headerLineColumns)

    for (cellName, outputName_for_cell) in elementsToWrite:
        cellColumns = [cellName, outputName_for_cell]

        ws.append(cellColumns)

    wb.save(CONFIG_OUTPUT_FILE)




output_cell_dict = ProcessOutputCellList(outputCellListFileName)


elementsToWrite = {}
allCellsInAllFolders = [] #elementsToWrite.keys() is not equivalent to the list of all cells for which there is a source folder because I do not add a uniqueID to that list if there's no metadata file for it, for example
for sourceFolder in [os.path.join(SOURCES_DIR, folder) for folder in sourceFolders]:
    print "reading sources folder: %s" % sourceFolder
    cellsInFolder = sorted([d for d in os.listdir(sourceFolder) \
                     if os.path.isdir(os.path.join(sourceFolder, d)) and d.lower().startswith("sample_lib") and d.lower() != "Undetermined_indices".lower() and d.lower() != "home".lower()])

    #remove the prefix "Sample_" from the cell name if it is there
    cellsInFolder = [cellName[len("Sample_"):] if cellName.startswith("Sample_") else cellName for cellName in cellsInFolder]

    allCellsInAllFolders += cellsInFolder

    for cellName in cellsInFolder:

        outputName_for_cell = output_cell_dict.get(cellName)
        if not(outputName_for_cell):
            print "Reading source data directories: Warning: cell %s has no collected output associated with it - skipping it..." % cellName
            continue

        if(elementsToWrite.has_key(cellName)):
            raise Exception("Reading source data directories: unique cell ID (%s) encountered twice?!")

        #the value associated with the key is a tuple: the uniqueID, its metadata and its output name
        elementsToWrite[cellName] = (cellName, outputName_for_cell)

    print "done reading sources folder: %s" % sourceFolder
    print "*********************************************\n\n"


numCellDirectoriesRead = len(allCellsInAllFolders)
allCellsInAllFolders = frozenset(allCellsInAllFolders)
if(numCellDirectoriesRead > len(allCellsInAllFolders)):
    raise Exception("The list of unique IDs read from the source folders is not unique...")


#now, see if there are any collect output entries for which no file was found
unclaimed_output_entries = sorted(frozenset(output_cell_dict.keys()).difference(allCellsInAllFolders))
for entry in unclaimed_output_entries:
    print "Warning: cell %s appeared in collect output list but had no source directory associated with it - skipping it..." % entry

#you collected all the elements to be written, now write the config file
WriteConfigFile(elementsToWrite)

print "all done"
