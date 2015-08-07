# -*- coding: utf-8 -*-
"""
Created on Sat Jan 03 21:29:39 2015

@author: daved_000
"""
import os;
import numpy as np;
import matplotlib;
import shutil;
matplotlib.use("svg")
import matplotlib.pyplot as plt;
from . import HtmlViewer;

def make_dirs(root_directory):
    """
    Creates relevant directories that will be needed for output files.

    :param root_directory: Root directory for FastProject output.
    :return: None
    """

    try:
        os.makedirs(root_directory);
    except OSError:
        pass;

    html_sub_directory = os.path.join(root_directory, HtmlViewer.OUTPUT_RESOURCE_DIRECTORY);
    try:
        os.makedirs(html_sub_directory);
    except OSError:
        pass;

def read_matrix(filename='', delimiter = '\t'):
    """Reads data in matrix format from <filename>  Returns data_matrix, 
    row_labels, and column_labels"""
    
    if(filename == ''):
        from Tkinter import Tk
        from tkFileDialog import askopenfilename
        Tk().withdraw();
        filename = askopenfilename();
    
    print("Loading data from " + filename + " ...");
    #load matrix data
    #First read in list of cell names
    
    with open(filename,'r') as ff:
        col_labels = ff.readline().strip().split(delimiter);
        col_labels = [col_label.strip() for col_label in col_labels];

        #col_labels in first line may or may not be shifted by one to line up with
        #data in following lines.  Determine # of datapoints from second line
        line2 = ff.readline().strip().split(delimiter);
        if(len(line2) == len(col_labels)):
            col_labels = col_labels[1:];


    #read all columns, start at 1 to skip first column of row labels
    cols_to_read = np.arange(1,len(col_labels)+1);
    
    data = np.loadtxt(filename,delimiter=delimiter,skiprows=1,usecols=cols_to_read);
    
    if(data.std() > 10 or data.mean() > 10):
        print "Log-Transforming data..."
        data = np.log(data + 1);
    
    #read in gene names
    with open(filename,'r') as ff:
        xx = ff.readline();

        row_labels = list();
        for line in ff:
          entries = line.split(delimiter);
          row_labels.append(entries[0].upper());
    

    #Test for unique row and column labels
    if(len(set(row_labels)) != len(row_labels)):
        #Find duplicates
        row_labels_copy = [r for r in row_labels];
        for row_label in set(row_labels_copy):
            row_labels_copy.remove(row_label);

        #Warning for duplicate genes
        print("WARNING: Row labels (gene identifiers) are not unique.\n" + ", ".join(row_labels_copy));

    if(len(set(col_labels)) != len(col_labels)):
        #Find duplicates
        for col_label in set(col_labels):
            col_labels.remove(col_label);

        raise ValueError("Column labels (sample names) are not unique.\n" + ", ".join(col_labels) + "\nExiting.");

    return (data, row_labels, col_labels);

def write_matrix(filename, data, row_labels, col_labels):
    ff = open(filename, 'w');
    
    ff.write('\t' + '\t'.join(col_labels) + '\n'); #Extra tab so col labels line up with data
    for i in range(data.shape[0]):
        ff.write(row_labels[i] + '\t');
        ff.write('\t'.join(["{:.5f}".format(num).rstrip('0').rstrip('.') for num in data[i,:]]));
        ff.write('\n');
    
    ff.close();  

def write_data(filename, data):
    write_matrix(filename, data, data.row_labels, data.col_labels); 
    
def read_matrix_nolabels(filename = '', delimiter = '\t'):
    
    if(filename == ''):
        from Tkinter import Tk
        from tkFileDialog import askopenfilename
        Tk().withdraw();
        filename = askopenfilename();

    #read data    
    data = np.loadtxt(filename, delimiter=delimiter);
    
    #generate labels
    num_rows = data.shape[0];
    num_cols = data.shape[1];
    
    row_labels = [str(i) for i in xrange(num_rows)];
    col_labels = [str(i) for i in xrange(num_cols)];
    
    return (data, row_labels, col_labels);
    
def write_signature_scores(filename, sig_scores_dict, col_labels):
    with open(filename, 'w') as ff:
        #Write column headers to first row
        ff.write('\t' + '\t'.join(col_labels) + '\n');

        #Iterate over signatures and write scores
        for name, sig_scores in sig_scores_dict.items():
            if(sig_scores.isFactor):
                row = [name] + sig_scores.scores;
            else:
                row = [name] + ["{:.5f}".format(num).rstrip('0').rstrip('.') for num in sig_scores.scores];

            ff.write('\t'.join(row) + '\n');


#def write_scatter_plot(filename, x_coords, y_coords, colors=[], xlabel='', ylabel='', title=''):
#    #Add .png extension if filename lacks it
#    if(not filename.endswith(".png")):
#        filename = filename + ".png";
#
#    plt.ioff();
#    ff = plt.figure();
#    if(len(colors) == len(x_coords)):
#        plt.scatter(x_coords, y_coords, c=colors);
#        plt.colorbar();
#    else:
#        plt.scatter(x_coords, y_coords);
#
#    plt.xlabel(xlabel);
#    plt.ylabel(ylabel);
#    plt.title(title);
#
#    plt.savefig(filename, bbox_inches='tight');
#
#    plt.close(ff);

def write_projection_file(filename, sample_labels, projections):
    """
    Outputs the coordinates for each projection to a file.

    Parameters
    ----------
    filename : String
        File to write
    sample_labels : list(String) len=Num_Samples
        Labels for each sample
    projections : dict(string, (2 x Num_Samples) numpy.ndarray)
        dictionary mapping the projection type (e.g. "tSNE") to an array containing
        the two-dimensional coordinates for each sample in the projection.

    """

    ff = open(filename,'w');

    for proj in projections.keys():
        coordinates = projections[proj];
        for i in range(coordinates.shape[1]):
            ff.write(proj + '\t');
            ff.write(sample_labels[i] + '\t');
            ff.write('{:.5f}'.format(coordinates[0,i]) + '\t');
            ff.write('{:.5f}'.format(coordinates[1,i]) + '\t');
            ff.write('\n');

    ff.close();

def write_cluster_file(filename, sample_labels, clusters):
    """
    Outputs cluster assignments to a file

    :param filename: String
        File to be written
    :param sample_labels: list(String)
        Labels for each sample
    :param clusters: Dict(projection_name (String) -> Projection Clusters)
        Projection Clusters: Dict(cluster_name (String) -> Cluster Assignments)
        Cluster Assignments: numpy.ndarray (len = num_samples)
            Cluster assignments for each sample

    """

    with open(filename, 'w') as fout:
        fout.write('\t'.join(sample_labels) + '\n');

        for proj in clusters.keys():
            proj_clusters = clusters[proj];
            for cluster in proj_clusters.keys():
                assignments = proj_clusters[cluster];
                fout.write(proj + " " + cluster + "\t");
                fout.write('\t'.join([str(x) for x in assignments]));
                fout.write('\n');

def write_filter_file(filename, filter):
    """
    Outputs cluster assignments to a file

    :param filename: String
        File to be written
    :param filter: set(String)
        Set of gene names included in the filter

    """

    with open(filename, 'w') as fout:
        fout.write('\t'.join(filter) + '\n');

def write_qc_file(directory, sample_passes, sample_scores, sample_labels):
    """
    Outputs a file with Quality Report information for the samples

    :param directory: String
        Output directory to use
    :param sample_passes: np.ndarray, dtype=bool, size=NUM_SAMPLES
        Whether or not the sample passes the quality check
    :param sample_scores: np.ndarray, dtype=float, size=NUM_SAMPLES
        A samples quality score
    :param sample_labels: list(String)
        Labels for each sample
    :return: None
    """

    #Build SVG
    svg_file_name = os.path.join(directory, HtmlViewer.OUTPUT_RESOURCE_DIRECTORY, "qc.svg");

    plt.style.use("fivethirtyeight");
    quantity, boundaries, patches = plt.hist(sample_scores, 30);

    ##Compute MAD
    med = np.median(sample_scores);
    mad = np.median(np.abs(sample_scores - med));
    cutoff = med - 1.6*mad;
    plt.plot([cutoff, cutoff], [0, np.max(quantity)*1.1], "red");


    plt.savefig(svg_file_name);

    #Build html Table
    ii = np.argsort(sample_scores);

    sample_scores = sample_scores[ii];
    sample_labels = [sample_labels[i] for i in ii];

    out_table = [];
    out_table.append("<table>");
    out_table.append("<tr>");
    out_table.append("<th>Sample</th>");
    out_table.append("<th>Score</th>");
    out_table.append("</tr>");

    for score, label in zip(sample_scores, sample_labels):
        if(score < cutoff):
            out_table.append('<tr class="failed_row">');
        else:
            out_table.append("<tr>");
        out_table.append("<td>" + label + "</td>");
        out_table.append("<td>" + "{:.3f}".format(score) + "</td>");
        out_table.append("</tr>");

    out_table.append("</table>");


    #Output Html File
    out_html = os.path.join(directory,"QC_Report.html");

    shutil.copy(HtmlViewer.RESOURCE_DIR + os.sep + "QC_Report.html",out_html);

    with open(out_html, "r") as fin:
        contents = fin.read();

    contents = contents.replace("<% table %>", "\n".join(out_table));

    with open(out_html, "w") as fout:
        fout.write(contents);





