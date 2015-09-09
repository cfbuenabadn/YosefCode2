# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 16:34:49 2015

@author: David
"""
from __future__ import division, print_function;

import numpy as np;
from sklearn.metrics.pairwise import pairwise_distances;
from scipy.spatial.distance import cdist;
from scipy.stats import norm, rankdata;
from .Utils import ProgressBar;
import HtmlViewer;

#This is used to cache the background distribution used when evaluating
#Signatures vs projections.  No need to regenerate the random indices
#when the size has not changed.  Saves significant time for large N_SAMPLES
_bg_dist = np.zeros((0,0));
def get_bg_dist(N_SAMPLES, NUM_REPLICATES):
    global _bg_dist;
    if(_bg_dist.shape[0] != N_SAMPLES or _bg_dist.shape[1] != NUM_REPLICATES):
        _bg_dist = np.random.rand(N_SAMPLES, NUM_REPLICATES);
        _bg_dist = np.argsort(_bg_dist, axis=0);

    return _bg_dist;

def read_signatures(filename='', match_terms=[]):
    """Calls either read_signatures_txt or read_signatures_gmt when appropriate
    """

    if(filename == ''):
        from Tkinter import Tk
        from tkFileDialog import askopenfilename
        Tk().withdraw();
        filename = askopenfilename();

    if(filename.lower().endswith('.gmt')):
        return read_signatures_gmt(filename, match_terms);
    else:
        return read_signatures_txt(filename, match_terms);


def read_signatures_txt(filename='', match_terms=[]):
    """Reads in a signature file.  Returns a list of Signature objects
    
    Parameters
    ----------
    filename : string
        Name (and path if not in working directory) of signature file to read
        If no file is entered, opens a file dialog
    match_terms : list(String)
        List of terms to be matched against signature names.
        If empty or omitted, reads all signatures
        Signature is retained if any of the terms match.
    
    Returns
    -------
    signatures : list(FastProject.Signature)
        The filtered signature list
    """
    
    if(filename == ''):
        from Tkinter import Tk
        from tkFileDialog import askopenfilename
        Tk().withdraw();
        filename = askopenfilename();
    
    ff = open(filename, 'r');
    found_signatures = dict();    ## dictionary of signatures
    
    if(type(match_terms) is str):
        match_terms = [match_terms];
    
    match_terms = [term.lower() for term in match_terms];
    
    try:    

        for i, line in enumerate(ff):
            sline = line.strip();
            if(sline[0] == '#' or sline[0:2] == '//'):
                continue;

            row_data = sline.split('\t');
            if(len(row_data) == 2):
                signed = False;
            elif(len(row_data) == 3):
                signed = True;
            else:
                raise ValueError("Line " + str(i) + " Signature file should contain 2 (unsigned) or 3 (signed) columns");

            name = row_data[0];
            if(not found_signatures.has_key(name)):  ## Add the signature if we haven't seen it yet

                ## Only add signatures if they match one of the supplied terms
                matched = False if len(match_terms) > 0 else True
                lname = name.lower();                
                for term in match_terms:
                    if(lname.find(term) >= 0):
                        matched = True;
                        break;
                
                if(matched):
                    sig = Signature(dict(), signed, filename, name);
                    found_signatures.update({name: sig});
                else:
                    continue;
            
            sig_dict = found_signatures[name].sig_dict;
            
            
            if signed:
                sig_sign = row_data[1].lower();
                gene_name = row_data[2].upper();
                if(sig_sign == 'plus'):
                    sig_val = 1;
                elif(sig_sign == 'minus'):
                    sig_val = -1;
                elif(sig_sign == 'both'):
                    sig_val = 0;
                elif(sig_sign == 'mius'):  #LOL, close enough
                    sig_val = -1;
                ## Other spelling variants go here
                else:
                    print("Error on line ", str(i), " Couldn't read signature value.");
                    print("   :", line);
                    continue;

            else:  #unsigned case
                sig_val = 0;
                gene_name = row_data[1].upper();

            sig_dict[gene_name] = sig_val;
    except:
        raise;
    finally:
        ff.close();
    
    return [found_signatures[key] for key in found_signatures.keys()]  #dict to list

def read_signatures_gmt(filename='', match_terms=[]):
    """Reads in a signature file in GMT format.

    Each row is tab delimited with gene set name, description, and list of genes
    Signed signatures should be on adjacent rows with identical set names except for
    a _plus or _minus suffix.

    Example:
    Memory_plus description gene1 gene2 gene3
    Memory_minus description gene1 gene2 gene3

    Returns a list of Signature objects

    Parameters
    ----------
    filename : string
        Name (and path if not in working directory) of signature file to read
        If no file is entered, opens a file dialog
    match_terms : list(String)
        List of terms to be matched against signature names.
        If empty or omitted, reads all signatures
        Signature is retained if any of the terms match.

    Returns
    -------
    signatures : list(FastProject.Signature)
        The filtered signature list
    """

    plus_terms = ['plus', 'up'];
    minus_terms = ['minus', 'dn', 'down'];

    def sig_parts(name):
        """
        Returns: The name of the signature with the plus or minus suffix removed.
                and whether or not the signature is positive or negative
                    : as a string, either 'plus' or 'minus'
        """
        lname = name.lower();
        sign = 'unsigned';
        ii = len(lname);
        for term in plus_terms:
            if(lname.endswith('_' + term)):
                ii = lname.rfind('_' + term);
                sign = 'plus';

        if(sign != 'plus'):
            for term in minus_terms:
                if(lname.endswith('_' + term)):
                    ii = lname.rfind('_' + term);
                    sign = 'minus';

        return name[0:ii], sign;

    if(filename == ''):
        from Tkinter import Tk
        from tkFileDialog import askopenfilename
        Tk().withdraw();
        filename = askopenfilename();

    if(type(match_terms) is str):
        match_terms = [match_terms];

    match_terms = [term.lower() for term in match_terms];

    with open(filename, 'r') as ff:
        found_signatures = dict();    ## dictionary of signatures
        end_of_file = False;

        while(not end_of_file):
            line = ff.readline().strip();
            if(line == ""):
                end_of_file = True;
                continue;
            if(line[0] == '#' or line[0:2] == '//'):
                continue;

            row_data = line.split('\t');
            name = row_data[0];
            root_name, sign = sig_parts(name);

            ## Only add signatures if they match one of the supplied terms
            matched = False if len(match_terms) > 0 else True
            lname = root_name.lower();
            for term in match_terms:
                if(term in lname):
                    matched = True;
                    break;
            if(not matched): continue;

            isSigned = sign == 'plus' or sign == 'minus';

            #Check if sig exists, create/insert if not
            if(found_signatures.has_key(root_name)):
                sig = found_signatures[root_name];
            else:
                sig = Signature(dict(), isSigned, filename, root_name);
                found_signatures.update({root_name: sig});

            if(sign == 'plus'):
                sig_val = 1;
            elif(sign == 'minus'):
                sig_val = -1;
            elif(sign == 'unsigned'):
                sig_val = 0;
            else:
                raise Exception("This should not happen");

            for gene_name in row_data[2:]:
                sig.sig_dict.update({gene_name.upper(): sig_val});

    return found_signatures.values();  #dict to list

def filter_sig_list(signatures, match_terms):
    """
    Filters the list of signature objects to retain only those that contain
    one of the specified terms
    
    Parameters
    ----------
    signatures : list(Signature)
        List of signatures to be filtered
    match_terms : list(String)
        List of terms to be matched against signature names.
        Signature is retained if any of the terms match.
    
    Returns
    -------
    filtered_signatures : list(FastProject.Signature)
        The filtered signature list
        
    
    """
    
    if(type(match_terms) is str):       #So you dont need to wrap a single string in []
        match_terms = [match_terms];    
    
    match_terms = map(lambda term: term.lower(), match_terms);    
    
    filtered_signatures = list();
    for sig in signatures:
        name = sig.name.lower();
        
        for term in match_terms:            
            if(name.find(term) >= 0):
                filtered_signatures.append(sig);
                break;
    
    
    return filtered_signatures;

def sigs_vs_projections(projections, sig_scores_dict, NEIGHBORHOOD_SIZE = 0.33, subsample_size = None):
    """
    Evaluates the significance of each signature vs each projection

    :param projections: dict of (string) => (numpy.ndarray of shape 2xNum_Samples)
        Maps projections to their spatial coordinates for each sample
    :param sig_scores: dict of (string) => (numpy.ndarray of shape Num_Samples)
        Maps signature names to their value at each coordinate
    :return:
    """
    sp_row_labels = [];
    sp_row_labels_factors = [];

    #Remove signatures that are factor signatures
    for name, sig_scores in sig_scores_dict.items():
        if(sig_scores.isFactor):
            sp_row_labels_factors.append(name);
        else:
            sp_row_labels.append(name);

    sp_col_labels = projections.keys();
    sp_col_labels.sort();

    N_SAMPLES = len(sig_scores_dict[sp_row_labels[0]].sample_labels);
    N_SIGNATURES = len(sp_row_labels);
    N_SIGNATURES_FACTORS = len(sp_row_labels_factors);
    N_PROJECTIONS = len(sp_col_labels);


    sig_proj_matrix   = np.zeros((N_SIGNATURES,N_PROJECTIONS));
    sig_proj_matrix_p = np.zeros((N_SIGNATURES,N_PROJECTIONS));

    factor_sig_proj_matrix   = np.zeros((N_SIGNATURES_FACTORS,N_PROJECTIONS));
    factor_sig_proj_matrix_p = np.zeros((N_SIGNATURES_FACTORS,N_PROJECTIONS));

    #Build a matrix of all signatures
    sig_score_matrix = np.zeros((N_SAMPLES, N_SIGNATURES));

    for j, sig in enumerate(sp_row_labels):
        sig_score_matrix[:,j] = sig_scores_dict[sig].ranks;

    #Build one-hot matrices for each factor
    factor_dict = dict();
    for sig in sp_row_labels_factors:
        factor_values = sig_scores_dict[sig].scores;
        factor_levels = list(set(factor_values)); #Makes unique
        factor_frequencies = np.zeros(len(factor_levels));
        factor_matrix = np.zeros((N_SAMPLES, 0));
        for j, fval in enumerate(factor_levels):
            factor_matrix_row = np.zeros((N_SAMPLES, 1));
            equal_ii = [i for i,val in enumerate(factor_values) if val == fval];
            factor_matrix_row[equal_ii] = 1;
            factor_frequencies[j] = len(equal_ii) / len(factor_values);
            factor_matrix = np.concatenate((factor_matrix, factor_matrix_row), axis=1);
        factor_dict[sig] = (factor_levels, factor_frequencies, factor_matrix);




    if(subsample_size):
        ii_sub = np.random.choice(N_SAMPLES, subsample_size, replace=False);
    else:
        ii_sub = np.arange(N_SAMPLES);

    print();
    print("Evaluating Signatures against Projections");
    pp = ProgressBar(N_PROJECTIONS);
    for i, proj in enumerate(sp_col_labels):
        data_loc = projections[proj];

        distance_matrix = cdist(data_loc[:,ii_sub].T, data_loc.T, metric='euclidean');

        weights = np.exp(-1 * distance_matrix**2 / NEIGHBORHOOD_SIZE**2);
        weights[np.arange(ii_sub.size), ii_sub] = 0; #Don't count self
        weights /= np.sum(weights, axis=1, keepdims=True);

        neighborhood_prediction = np.dot(weights, sig_score_matrix);


        ##Neighborhood dissimilarity score = |actual - predicted|
        dissimilarity = np.abs(sig_score_matrix[ii_sub,:] - neighborhood_prediction);
        med_dissimilarity = np.median(dissimilarity, axis=0);

        NUM_REPLICATES = 10000;

        #random_sig_values of a given size is cached in this module
        random_sig_values = get_bg_dist(N_SAMPLES, NUM_REPLICATES);

        random_predictions = np.dot(weights, random_sig_values);
        random_scores = np.median(np.abs(random_sig_values[ii_sub,:] - random_predictions), axis=0);

        mu = np.mean(random_scores);
        sigma = np.std(random_scores);
        if(sigma == 0):
            sigma = 1e-3;

        p_values = norm.cdf((med_dissimilarity - mu)/sigma);


        sig_proj_matrix[:,i] = med_dissimilarity;
        sig_proj_matrix_p[:,i] = p_values;

        #Calculate significance for Factor signatures
        for j,sig in enumerate(sp_row_labels_factors):
            factor_levels, factor_frequencies, factor_matrix = factor_dict[sig];
            N_LEVELS = len(factor_levels);
            factor_predictions = np.dot(weights, factor_matrix);

            if(subsample_size):
                factor_values = factor_matrix[ii_sub,:];
            else:
                factor_values = factor_matrix;

            dissimilarity = 1 - np.sum(factor_values * factor_predictions, axis=1);
            med_dissimilarity = np.median(dissimilarity);

            #Now...compute a background?
            NUM_REPLICATES = 1000;
            rand_factors = np.random.rand(N_SAMPLES, NUM_REPLICATES);
            column_assignments = np.random.choice(N_LEVELS, NUM_REPLICATES, p = factor_frequencies);
            column_assignments = factor_frequencies[column_assignments];
            column_assignments = column_assignments.reshape((1,NUM_REPLICATES));
            rand_factors = (rand_factors < column_assignments).astype('int');
            rand_med_dissimilarity = np.median(1-np.dot(weights, rand_factors), axis=0);

            mu = np.mean(rand_med_dissimilarity);
            sigma = np.std(rand_med_dissimilarity);
            if(sigma == 0):
                p_value = 1;
            else:
                p_value = norm.cdf((med_dissimilarity - mu)/sigma); #Runtime error here possibly?  Divide by zero?

            factor_sig_proj_matrix[j,i] = med_dissimilarity;
            factor_sig_proj_matrix_p[j,i] = p_value;

        pp.update();

    #Concatenate the Factor sig-proj entires back in
    sig_proj_matrix = np.concatenate((sig_proj_matrix, factor_sig_proj_matrix), axis=0);
    sig_proj_matrix_p = np.concatenate((sig_proj_matrix_p, factor_sig_proj_matrix_p), axis=0);
    sp_row_labels = sp_row_labels + sp_row_labels_factors;

    sig_proj_matrix_p = p_to_q(sig_proj_matrix_p);
    sig_proj_matrix_p[sig_proj_matrix_p == 0] = 1e-300; #Correct for -inf
    sig_proj_matrix_p = np.log10(sig_proj_matrix_p);


    pp.complete();

    return (sp_row_labels, sp_col_labels, sig_proj_matrix, sig_proj_matrix_p);

def p_to_q(p_values):
    """
    Uses the Benjamini-Hochberg procedure to convert p_values to q_values

    :param p_values:  numpy.ndarray of p_values
    :return q_values: numpy.ndarray of q_values, same shape as p_values
    """
    original_shape = p_values.shape;
    p_vals_flat = p_values.flatten();
    rank = p_vals_flat.argsort().argsort()+1;
    num_tests = p_values.size;
    q_vals = p_vals_flat * num_tests / rank;
    q_vals.shape = original_shape;

    return q_vals;

def load_precomputed(filename, sample_labels):
    """
    Reads precomputed signature values form a tab-delimited text file
    First row of the file contains sample labels that the signatures correspond with
    Each subsequent row contains a signature name in the first column,
         followed by the signature type (either 'numerical' or 'factor')
         followed by the signature values, one for each sample label in the file

    :param filename: signature score file name
    :param sample_labels: labels for which we want the signature scores
    :return: Dictionary of Signature Name (String) -> Signature (SignatureScores)
    """

    with open(filename, 'r') as fin:
        #Determine column labels to apply
        line1 = fin.readline().rstrip().split('\t');
        line2 = fin.readline().rstrip().split('\t');

        if(len(line1) == len(line2)): #First two entries must be empty or column headers
            line1 = line1[2:];
        elif(len(line1) == len(line2)-2):
            line1 = line1;
        else: #Other arrangements signify some unusual formatting
            raise ValueError("Error in header line of precomputed signature file.\n"
                 + "First row should contain tab-separated list of samples");


        #match indices between signatures in file and sample_labels
        #want x such that file_cols[x] == sample labels
        target_l = [sl.lower() for sl in sample_labels];
        source_l = [sl.lower() for sl in line1];
        translation_indices = np.zeros(len(target_l), dtype=np.int32);
        for i in xrange(translation_indices.size):
            try:
                translation_indices[i] = source_l.index(target_l[i]);
            except ValueError:
                raise ValueError("Error: Missing value in precomputed signatures for sample " + target_l[i]);

        fin.seek(0);
        xx = fin.readline();

        #Gather signatures
        sig_scores = dict();
        for line in fin:
            line = line.strip();
            if(line == ""): continue;
            s_line = line.split("\t");
            sig_name = s_line[0];
            sig_type = s_line[1].strip().lower();
            sig_val_cells = s_line[2:];

            if(sig_type == 'numerical'):
                sig_isFactor = False;
                try:
                    sig_vals = np.array([float(x) for x in sig_val_cells]);
                    sig_vals = sig_vals[translation_indices];
                except ValueError as e:
                    print(e.message);
                    print('Error in precomputed signature:', sig_name);
                    for i,x in enumerate(sig_val_cells):
                        try:
                            y = float(x);
                        except ValueError:
                            print("Error in column", i);
                            print("Bad value:", x);
                    raise Exception('Failed to load precomputed signature. Correct file format and re-run.');
            elif(sig_type == 'factor'):
                sig_isFactor = True;
                sig_vals = [sig_val_cells[i] for i in translation_indices];
            else:
                raise ValueError('Column 2 of precomputed signature file should specify either "numerical" or "factor"');

            sig_scores[sig_name] = SignatureScores(sig_vals, sig_name, sample_labels, sig_isFactor, isPrecomputed=True);

        return sig_scores;

class Signature:
    
    def __init__(self, sig_dict, signed, filename, name):
        self.sig_dict = sig_dict;
        self.signed = signed;
        self.source = filename;
        self.name = name;

    def sig_indices(self, genes):
        """Helper method
        
        Returns an array with length = len(genes)
        Entries in the array are 0 if the gene is not in the signature
        Otherwise, value is determined by the signature type"""
        
        out = np.zeros((len(genes), 1), dtype=np.float64);
        
        for i, gene in enumerate(genes):
            if(self.sig_dict.has_key(gene)):
                val = self.sig_dict[gene];
                if(val == 1 or val == 0):
                    out[i] = 1;
                if(val == -1):
                    out[i] = -1;

        #If signed, weight the indices such that the sum of positive signatures
        #counts as much as the sum of negative signatures
        #Weights result in mean(data[pos_sig])/2 - mean(data[neg_sig])
        #                   =  mean(data*sig)
        #Results in mean(out) = 0
        #           mean(|out|) = 1 
        if(self.signed):
            num_pos = np.count_nonzero(out ==  1);
            num_neg = np.count_nonzero(out == -1);
            if(num_pos > 0 and num_neg > 0):
                num_total = num_pos + num_neg;
                
                pos_weight = num_total/num_pos/2;
                neg_weight = num_total/num_neg/2;
    
                out[out==1]  = pos_weight;
                out[out==-1] = neg_weight*-1;
         
        return out;

class SignatureScores:
    """
    Represents a Signature evaluated on a set of samples
    """

    @property
    def ranks(self):
        if(self.isFactor):
            raise Exception("Factor signature scores have no rank")

        if(self._ranks is None):
            self._ranks = rankdata(self.scores, method="average");

        return self._ranks;


    def __init__(self, scores, name, sample_labels, isFactor, isPrecomputed):
        self.scores = scores;
        self.name = name;
        self.sample_labels = sample_labels;
        self.isFactor = isFactor;
        self.isPrecomputed = isPrecomputed;
        self._ranks = None;

    def to_JSON(self):
        """
        Construct a dictionary of certain parameters.
        Parse that to JSON using HTMLViewer and return
        :return: String with JSON representation of the SignatureScores instance
        """

        out = dict({
            "name": self.name,
            "scores": self.scores,
            "isFactor": self.isFactor,
            "isPrecomputed": self.isPrecomputed,
        });

        if(not self.isFactor):
            out.update({"ranks": self.ranks});

        return HtmlViewer.toJS(out);



