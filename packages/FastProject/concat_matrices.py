from __future__ import division, print_function;
import pandas as pd
import numpy as np
import logging
import sys
import os

out_file = sys.argv[1];
files = sys.argv[2:];

for filename in files:
    if(not os.path.isfile(filename)):
        raise ValueError("File does not exist:", filename);

if(os.path.isfile(out_file)):
    raise ValueError("Output file,", out_file, "already exists.  Please remove before running the script.");

#Drop duplicate indices
def drop_dup_indices(df):
    df['index_to_drop'] = df.index;
    df.drop_duplicates(subset='index_to_drop', take_last=False, inplace=True);
    del df['index_to_drop']

#First file
print("Reading", files[0], "...");
output = pd.read_table(files[0], header=0, index_col = 0);
drop_dup_indices(output);

#Merge in other files
for filename in files[1:]:
    print("Reading", filename, "...");
    new_data = pd.read_table(filename, header=0, index_col = 0);
    print("Merging", filename, "...");
    drop_dup_indices(new_data);
    output = pd.concat((output, new_data), join='outer', axis=1);

output[output == 0] = np.nan;

print("Making column names unique");
#Make columns unique
col_names = output.columns.values;
unique_col_names = set(col_names.tolist());
for i, name in enumerate(col_names):
    if(name in unique_col_names):
        new_name = name + "x";
        while(new_name in unique_col_names):
            new_name = new_name + "x";
        col_names[i] = new_name;

output.columns = col_names;

print("Writing result to", out_file);
#Output to one file
output.to_csv(out_file, sep='\t', na_rep='0');
