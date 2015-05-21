
ncRNA_dict = {}
MM10_GFF = '/data/yosef/index_files/mm10_4brain/index/GCF_000001635.23_GRCm38.p3_genomic.gff'
with open(MM10_GFF) as fin:
    line = fin.readline()

    while(line):
        if("ncrna_class" in line):
            #from stack overflow: initialize a dictionary from a string (http://stackoverflow.com/questions/186857/splitting-a-semicolon-separated-string-to-a-dictionary-in-python)
            # s= "Name1=Value1;Name2=Value2;Name3=Value3"
            # dict(item.split("=") for item in s.split(";"))

            tags = line.split('\t')[-1].split(';')
            tags = dict(item.split('=') for item in tags)

            #dictionary gene name --> ncrna_class
            geneName = tags['gene']
            ncRnaClass = tags['ncrna_class']

            if(not(geneName) in ncRNA_dict):
                ncRNA_dict[geneName] = ncRnaClass
            else:
                if(ncRNA_dict[geneName] != ncRnaClass):
                    print "error: mismatch in ncrna_class"

        line = fin.readline()


MM10_GFF_NCRNA_TABLE = '/data/yosef/index_files/mm10_4brain/derived_tables/ncRNA_table.txt'
with open(MM10_GFF_NCRNA_TABLE, 'wt') as outFile:
    outFile.writelines("Gene Name\tncRNA class\n")
    outFile.writelines("%s\t%s\n" % (item[0], item[1]) for item in ncRNA_dict.items())



    for item in ncRNA_dict.items():
        outFile.writelines("%s\t%s\n" % (item[0], item[1]))

