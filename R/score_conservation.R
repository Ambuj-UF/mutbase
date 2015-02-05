
################################################################################################################
#                                                                                                              #
# Copyright (C) {2014}  {Ambuj Kumar, Kimball-Brain lab group, Biology Department, University of Florida}      #
#                       {Mrinal Mishra, University of Turku, Finland}                                          #
#                                                                                                              #
# This program is free software: you can redistribute it and/or modify                                         #
# it under the terms of the GNU General Public License as published by                                         #
# the Free Software Foundation, either version 3 of the License, or                                            #
# (at your option) any later version.                                                                          #
#                                                                                                              #
# This program is distributed in the hope that it will be useful,                                              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                #
# GNU General Public License for more details.                                                                 #
#                                                                                                              #
# This program comes with ABSOLUTELY NO WARRANTY;                                                              #
# This is free software, and you are welcome to redistribute it                                                #
# under certain conditions;                                                                                    #
#                                                                                                              #
################################################################################################################


# Based on:
# Capra JA and Singh M. Predicting functionally important residues from sequence conservation.
# Bioinformatics, 23(15):1875-82, 2007.
#
# Scoring protein sequence conservation using the Jensen-Shannon divergence.
#


library(hash)
PSEUDOCOUNT = .0000001
amino_acids = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-')
iupac_alphabet = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z", "X", "*", "-")

aa_to_index = hash(keys=amino_acids,values=1:length(amino_acids))




################################################################################
# Frequency count
################################################################################
weighted_freq_count_pseudocount=function(col, seq_weights, pc_amount){
    #""" Return the weighted frequency count for a column--with pseudocount."""
        
        # if the weights do not match, use equal weight
        if (length(seq_weights) != length(col)){
seq_weights = rep(as.double(1.0),length(col))}
    
    aa_num = 1
        freq_counts = rep(pc_amount,length(amino_acids)) # in order defined by amino_acids
            
            for (aa in amino_acids){
for (j in 1:length(col)){
    if (col[[j]] == aa){
        freq_counts[aa_num] = freq_counts[[aa_num]] + rep(1,seq_weights[[j]])}
        aa_num = aa_num + 1}}
            for (j in 1:length(freq_counts)){
freq_counts[j] = freq_counts[[j]] / (sum(seq_weights) + length(amino_acids) * pc_amount)}
    return (freq_counts)
    }


################################################################################
# Gap Penalty
################################################################################
weighted_gap_penalty=function(col, seq_weights){
  """ Calculate the simple gap penalty multiplier for the column. If the 
    sequences are weighted, the gaps, when penalized, are weighted 
    accordingly. """

  # if the weights do not match, use equal weight
  if (length(seq_weights) != length(col)){
     seq_weights = rep(as.double(1.0),length(col))}

  gap_sum =0
  for (i in 1:length(col)){
     if (col[i] == '-'){
        gap_sum += seq_weights[i]}}

  return (1 - (gap_sum / sum(seq_weights)))
}



gap_percentage <- function(col){

   num_gaps = 0
   for (aa in col){
     if (aa == '-'){ num_gaps = num_gaps + 1}
   }

   return (num_gaps / length(col))
}
  



################################################################################
# Shannon Entropy
################################################################################


shannon_entropy <- function(col, sim_matrix, bg_distr, seq_weights, gap_penalty=1) {
    """Calculates the Shannon entropy of the column col."""

    fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

    h = 0.
    for (i in 1:length(fc)) {
        if fc[[i]] != 0 {
            h = h + fc[[i]] * log(fc[[i]])
        }
    }
    
    h = h/log(min(length(fc), length(col)))
    inf_score = 1 - (-1 * h)

        
    if (gap_penalty == 1) {
        return inf_score * weighted_gap_penalty(col, seq_weights)
    }
    else {
        return inf_score
    }
}



################################################################################
# Property Entropy
################################################################################

property_entropy <- function(col, sim_matrix, bg_distr, seq_weights, gap_penalty=1) {
    """Calculate the entropy of a column col relative to a partition of the amino acids. Similar to Mirny '99. sim_matrix and bg_distr are ignored, but could be used to define the sets."""
    
    # Mirny and Shakn. '99
     property_partition = list(list('A','V','L','I','M','C'), list('F','W','Y','H'), list('S','T','N','Q'), list('K','R'), list('D', 'E'), list('G', 'P'), list('-'))
     
     # Williamson '95
    
     # property_partition = [['V','L', 'I','M'], ['F','W','Y'], ['S','T'], ['N','Q'], ['H','K','R'], ['D','E'], ['A','G'], ['P'], ['C'], ['-']]
     
     fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)
     
     # sum the aa frequencies to get the property frequencies
     fc = rep(0, length(property_partition))
     for (p in 1:length(property_partition)) {
         for (aa in property_partition[[p]]) {
             prop_fc[p] = prop_fc[[p]] + fc[[values(aa_to_index[aa])]]
         }
     }
     
     h = 0
     
     for (i in 1:length(prop_fc)) {
         if (prop_fc[[i]] != 0) {
             h = h + prop_fc[[i]] * log(prop_fc[[i]])
         }
     }
     
     h = h/log(min(length(property_partition), length(col)))
     
     if (gap_penalty == 1) {
         return(inf_score * weighted_gap_penalty(col, seq_weights))
     }
     else {
         return(inf_score)
     }
}



################################################################################
# Property Relative Entropy
################################################################################


property_relative_entropy <- function(col, sim_matrix, bg_distr, seq_weights, gap_penalty=1) {
    
    # Mirny and Shakn. '99
    #property_partition = [['A','V','L','I','M','C'], ['F','W','Y','H'], ['S','T','N','Q'], ['K','R'], ['D', 'E'], ['G', 'P'], ['-']]
    
    # Williamson '95
    property_partition = list(list('V','L', 'I','M'), list('F','W','Y'), list('S','T'), list('N','Q'), list('H','K','R'), list('D','E'), list('A','G'), list('P'), list('C'))
    
    prop_bg_freq = c()
    
    if len(bg_distr) == len(property_partition) {
        prop_bg_freq = bg_distr
    }
    else {
        prop_bg_freq = c(0.248, 0.092, 0.114, 0.075, 0.132, 0.111, 0.161, 0.043, 0.024, 0.000)
    }
    
    #fc = weighted_freq_count_ignore_gaps(col, seq_weights)
    
    fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)
    
    # sum the aa frequencies to get the property frequencies
    prop_fc = rep(0, length(property_partition))
    for (p in 1:length(property_partition)) {
        for (aa in property_partition[[p]]) {
            prop_fc[p] = prop_fc[[p]] + fc[[aa_to_index[[aa]]]]
        }
    }
    
    d = 0
    for (i in 1:length(prop_fc)) {
        if (prop_fc[[i]] != 0 & prop_bg_freq[[i]] != 0) {
            d = d + prop_fc[[i]] * log(prop_fc[[i]] / prop_bg_freq[[i]], 2)
        }
    }
    
    if (gap_penalty == 1) {
        return (d * weighted_gap_penalty(col, seq_weights))
    }
    
    else {
        return(d)
    }
}


################################################################################
# von Neumann Entropy
################################################################################


vn_entropy <- function(col, sim_matrix, bg_distr, seq_weights, gap_penalty=1) {
    
    aa_counts = rep(0, 20)
    for (aa in col) {
        if (aa != '-') {aa_counts[[values(aa_to_index[aa])]] = aa_counts[[values(aa_to_index[aa])]] + 1}
    }
    
    dm_size = 0
    dm_aas = c()
    
    for (i in 1:length(aa_counts)) {
        if (aa_counts[[i]] != 0) {
            dm_aas = c(dm_aas, i)
            dm_size = dm_size + 1
        }
    }
    
    if (dm_size == 0) { return(0.0)}
    
    row_i = 0
    col_i = 0
    
    dm = matrix(0, nrow=length(dm_aas), ncol=length(dm_aas))
    
    for (i in 1:dm_size) {
        row_i = dm_aas[[i]]
        for (j in 1:dm_size) {
            col_i = dm_aas[[j]]
            dm[i, j] = aa_counts[[row_i]] * sim_matrix[[row_i]][[col_i]]
        }
    }
    
    ev = eigen(dm)
    
    ev_vector = c()
    for (x in ev){for (j in x) {ev_vector <- c(ev_vector, j)}}
    
    temp = 0
    for (e in ev_vector) {temp = temp + e}
    
    if (temp != 0) {
        for (i in 1:length(ev_vector) {
            ev_vector[i] = ev_vector/temp
        }
    }
    
    vne = 0
    
    for (e in ev_vector) {
        if e > (10**-10) {
            vne = vne - e * log(e) / log(20)
        }
    }
    
    if (gap_penalty == 1) {
        #return (1-vne) * weighted_gap_penalty(col, seq_weights)
        return ((1-vne) * weighted_gap_penalty(col, rep(1, length(col))))
    }
    else {
        return(1-vne)
    }
}



################################################################################
# Relative Entropy
################################################################################

relative_entropy <- function(col, sim_matix, bg_distr, seq_weights, gap_penalty=1) {
    ###Calculate the relative entropy of the column distribution with a background distribution specified in bg_distr. This is similar to the approach proposed in Wang and Samudrala 06. sim_matrix is ignored."""
    
    distr = bg_distr
    fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)
    
    # remove gap count
    if (length(distr) == 20) {
        new_fc = fc[-length(fc)]
    }
    
    s = sum(new_fc)
    for (i in 1:length(new_fc)) {
        new_fc[i] = new_fc[[i]] / s
    }
    
    fc = new_fc
    if (length(fc) != length(distr) { return(-1) }
    
    d = 0
    for (i in 1:length(fc)) {
        if distr[[i]] != 0.0 {
            d = d + fc[[i]] * log(fc[[i]]/distr[[i]])
        }
    }
    
    d = d / log(length(fc))
    
    if (gap_penalty == 1) {
        return (d * weighted_gap_penalty(col, seq_weights))
    }
    else {
        return(d)
    }
}



################################################################################
# Jensen-Shannon Divergence
################################################################################

js_divergence <- function(col, sim_matrix, bg_distr, seq_weights, gap_penalty=1) {
    ###Return the Jensen-Shannon Divergence for the column with the background distribution bg_distr. sim_matrix is ignored. JSD is the default method."""
    
    distr = bg_distr
    
    fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)
    
    # if background distrubtion lacks a gap count, remove fc gap count
    if (length(distr) == 20) {
        new_fc = fc[-length(fc)]
    }
    
    s = sum(new_fc)
    for (i in 1:length(new_fc)) {
        new_fc[[i]] = new_fc[[i]] / s
    }
    
    fc = new_fc
    
    if (length(fc) != length(distr)) { return(-1) }
    
    # make r distriubtion
    r = c()
    for (i in 1:length(fc)) {
        r = c(r, 0.5 * fc[[i]] + 0.5 * distr[[i]])
    }
    
    d = 0
    for (i in 1:length(fc)) {
        if (r[[i]] != 0.0) {
            if (fc[[i]] == 0.0) {
                d = d + distr[[i]] * log(distr[[i]]/r[[i]], 2)
            }
            else if (distr[[i]] == 0.0) {
                d = d + fc[[i]] * log(fc[[i]]/r[[i]], 2)
            }
            else {
                d = d + fc[[i]] * log(fc[[i]]/r[[i]], 2) + distr[[i]] * math.log(distr[[i]]/r[[i]], 2)
            }
            
        }
    }
    
    # d /= 2 * math.log(len(fc))
    d = d / 2
    if (gap_penalty == 1) {
        return(d * weighted_gap_penalty(col, seq_weights))
    }
    else {
        return(d)
    }
}





################################################################################
# Mutation Weighted Pairwise Match
################################################################################

sum_of_pairs <- function(col, sim_matrix, bg_distr, seq_weights, gap_penalty=1) {
    """ Sum the similarity matrix values for all pairs in the column. This method is similar to those proposed in Valdar 02. bg_distr is ignored."""
    
    sum = 0
    max_sum = 0
    for (i in 1:length(col)) {
        for (j in 1:length(i)) {
            if (col[[i]] != '-' and col[[j]] != '-') {
                max_sum = max_sum + seq_weights[[i]] * seq_weights[[j]]
                sum = sum + seq_weights[[i]] * seq_weights[[j]] * sim_matrix[[values(aa_to_index[col[[i]]])]][[values(aa_to_index[col[[j]]])]]
            }
        }
    }
    
    if (max_sum != 0) {
        sum = sum / max_sum
    }
    else {
        sum = 0
    }
    
    if (gap_penalty == 1) {
        return(sum * weighted_gap_penalty(col, seq_weights))
    }
    else {
        return(sum)
    }
}



################################################################################
# Window Score
################################################################################

window_score <- function(scores, window_len, lam=.5) {
    #""" This function takes a list of scores and a length and transforms them so that each position is a weighted average of the surrounding positions. Positions with scores less than zero are not changed and are ignored in the calculation. Here window_len is interpreted to mean window_len residues on either side of the current residue. """
    
    w_scores = scores
    for (i in window_len:(length(scores) - window_len)) {
        if (scores[[i]] < 0) {
            next
        }
        
        sum = 0
        num_terms = 0
        for (j in (i - window_len): (i + window_len + 1)) {
            if (i != j & scores[[j]] >= 0) {
                num_terms = num_terms + 1
                sum = sum + scores[j]
            }
        }
        
        if (num_terms > 0) {
            w_scores[i] = (1 - lam) * (sum / num_terms) + lam * scores[[i]]
        }
    }
    
    return(w_scores)
}



calc_z_scores <- function(scores, score_cutoff) {
    #"""Calculates the z-scores for a set of scores. Scores below score_cutoff are not included."""
    
    average = 0
    std_dev = 0
    z_scores = c()
    num_scores = 0
    
    for (s in scores) {
        if (s > score_cutoff) {
            average = average + s
            num_scores = num_scores + 1
    }
    }
    
    if (num_scores != 0) {
        average = average / num_scores
    }
    
    for (s in scores) {
        if (s > score_cutoff) {
            std_dev = std_dev + ((s - average)**2) / num_scores
}
}
    
    std_dev = sqrt(std_dev)
    for (s in scores) {
        if (s > score_cutoff & std_dev != 0) {
            z_scores = c(z_scores, (s-average)/std_dev)
        }
        else {
            z_scores = c(z_scores, -1000.0)
    }
    }
    
    return(z_scores)
}


################################################################################
################################################################################
################################################################################
#  END CONSERVATION SCORES
################################################################################
################################################################################
################################################################################


read_scoring_matrix <- function(sm_file) {
    #""" Read in a scoring matrix from a file, e.g., blosum80.bla, and return it as an array. """
    
    aa_index = 0
    first_line = 1
    row = c()
    list_sm = c()

    tryCatch( {
        
        text <- readLines(sm_file,encoding="UTF-8")
        
        for (line in text) {
            if (strsplit(line, "")[[1]][1] != '#' & first_line) {
                first_line = 0
                if (length(amino_acids) == 0) {
                    for (c in strsplit(line, " ")) {
                        values(aa_to_index[tolower(c)]) = aa_index
                        amino_acids = c(amino_acids, tolower(c))
                        aa_index = aa_index + 1
                    }
                }
            }
            else if (strsplit(line, "")[[1]][1] != '#' & first_line == 0) {
                if (nchar(line) > 1) {
                    row = strsplit(line, " ")
                    list_sm = c(list_sm, row)
                }
            }
        } 
        error = function(err) {
            return(identity(20))
        } 
        finally = {
            sprintf("Could not load similarity matrix: %s. Using identity matrix...", sm_file)
        }
    } )
    
    # if matrix is stored in lower tri form, copy to upper
    if (length(list_sm[[1]]) < 20) {
        for (i in 1:19) {
            for (j in (i+1):20) {
                list_sm[i] = c(list_sm[i], list_sm[[j]][i])
            }
        }
    }
    
    return(list_sm)
}



calculate_sequence_weights <- function(msa) {
    """ Calculate the sequence weights using the Henikoff '94 method for the given msa. """
    
    seq_weights = rep(0, length(msa))
    
    for (i in 1:length(msa[[1]])) {
        freq_counts = rep(0, length(amino_acids))
        col = c()
        for (j in 1:length(msa)) {
            if (msa[[j]][i] != '-') {
                freq_counts[[values(aa_to_index[msa[[j]][i]])]] = freq_counts[[values(aa_to_index[msa[[j]][i]])]] + 1
            }
        }
        
        num_observed_types = 0
        for (j in 1:length(freq_counts)) {
            if (freq_counts[[j]] > 0) {
                num_observed_types = num_observed_types + 1
            }
        }
        
        for (j in 1:length(msa)) {
            d = freq_counts[[aa_to_index[[msa[[j]][i]]]]] * num_observed_types
            if d > 0 {
                seq_weights[j] = seq_weights[[j]] + 1/d
            }
        }
    }
    
    for w in range(len(seq_weights)) {
        seq_weights[w] = seq_weights[[w]] / length(msa[[1]])
    }
    
    return(seq_weights)
}


load_sequence_weights <- function(fname) {
    seq_weights = c()
    f = readLines(fname, encoding="UTF-8")
    for (line in f) {
        l = strsplit(line, " ")
        if (line[[1]][1] != '#' & nchar(l[[1]]) == 2) {
            seq_weights = c(seq_weights, l[[2]])
        }
    }
    
    return(seq_weights)
}


get_column <- function(col_num, alignment) {
    col = c()
    for (seq in alignment[[1]]) {
        if (col_num <= nchar(seq)) {
            seq = strsplit(seq, "")[[1]]
            col = c(col, seq[[col_num]])
        }
    }
    
    return(col)
}


map_float <- function(vec_Obj) {
    for (i in 1:length(vec_Obj)) {
        vec_Obj[i] = as.numeric(vec_Obj[[i]])
    }
    
    return(vec_Obj)
}


get_distribution_from_file <- function(fname) {
    """ Read an amino acid distribution from a file. The probabilities should be on a single line separated by whitespace in alphabetical order as in amino_acids above. # is the comment character."""
    
    distribution = c()
    f = readLines(fname, encoding="UTF-8")
    for (line in f) {
        if (line[[1]] == '#') {
            line = line[-nchar(line)]
            distribution = strsplit(line)
            distribution = map_float(distribution)
        }
    }
    
    if (.997 > sum(distribution) or sum(distribution) > 1.003) {
        sprintf("Distribution does not sum to 1. Using default (BLOSUM62) background.")
        sprintf(sum(distribution))
        return(c())
    }
    
    return(distribution)
}



#########################################################################
# Utility functions
#########################################################################


replace_element <- function(stringElement, Obj1, Obj2) {
    new_vec = strsplit(stringElement, "")
    for (i in 1:length(new_vec[[1]])) {
        if (new_vec[[1]][[i]] == Obj1) {
            new_vec[[1]][i] = Obj2
        }
    }
    
    retString = paste(new_vec[[1]], collapse = "")
    
    return(retString)
}


rem_last <- function(stObj) {
    x = strsplit(stObj, "")[[1]]
    x = x[-length(x)]
    strRet = paste(x, collapse="")
    return(strRet)
}


is.defined <- function(x) !is.null(x)



trim_seq <- function(stObj, initPos=NULL, endPos=NULL) {
    newSeqList = strsplit(stObj, "")[[1]]
    if (is.defined(initPos) == TRUE) {
        removePosInit = c(1:initPos)
    }
    else { removePosInit = c() }
    if (is.defined(endPos) == TRUE) {
        removePosEnd = c(endPos:length(newSeqList))
    }
    else { removePosEnd = c() }
    
    removers = c(removePosInit, removePosEnd)
    
    newSeqList = newSeqList[-removers]
    newSeq = paste(newSeqList, collapse="")
    
    return(newSeq)
    
}



search_element <- function(stObj, search_obj) {
    if ( search_obj %in% strsplit(stObj, "")[[1]] ) {return(TRUE)}
    else {return(FALSE)}
}




###########################################################################################
# Read fasta files
###########################################################################################

read_fasta_alignment <- function(filename){

  f = readLines(filename,encoding="UTF-8")

  names = c()
  alignment = c()
  cur_seq = ''

  for (line in f){
    line = rem_last(line)
    if (nchar(line) == 0){ next }
    if (line[[1]] == ';'){ next }
    if (search_element(line[[1]], '>')) {
        names = c(names, trim_seq(line[[1]], initPos=1)[[1]])

        if (cur_seq != '') {
            cur_seq = toupper(cur_seq)
            aminos = strsplit(cur_seq, "")[[1]]
            for (i in 1:length(aminos)) {
                if (!(aminos[[i]] %in% iupac_alphabet)) {
                    aminos[i] = "-"
                }
            }
            
            cur_seq = paste(aminos, collapse="")
  
            cur_seq = replace_element(cur_seq,'B','D')
            cur_seq = replace_element(cur_seq,'Z','Q')
            cur_seq = replace_element(cur_seq,'X','-')
            alignment = c(alignment, cur_seq)
            cur_seq = ''
        }
    }
    
    else if (strsplit(line, "")[[1]][1] %in% iupac_alphabet){
        cur_seq = paste(cur_seq, line[[1]], sep="")
    }
    
  }
  
  cur_seq = toupper(cur_seq)
  aminos = strsplit(cur_seq, "")[[1]]
  for (i in 1:length(aminos)) {
    if (!(aminos[[i]] %in% iupac_alphabet)) {
        aminos[i] = "-"
    }
  }
  
  cur_seq = paste(aminos, collapse="")
  cur_seq = replace_element(cur_seq,'B','D')
  cur_seq = replace_element(cur_seq,'Z','Q')
  cur_seq = replace_element(cur_seq,'X','-')
  
  alignment = c(alignment, cur_seq)
  return (list(names,alignment))
  
}





################################################################################
# Begin execution
################################################################################



execute_conserve <- function(infile_name,
                             outfile_name,
                             window_size=3,
                             win_lam=0.5,
                             seq_weights=NULL,
                             s_matrix_file="/blosum62.bla",
                             bg_distribution=blosum_background_distr,
                             scoring_function=js_divergence,
                             use_seq_weights=TRUE,
                             background_name='blosum62',
                             gap_cutoff=0.3,
                             use_gap_penalty=1,
                             seq_specific_output=NULL,
                             normalize_scores=FALSE
                             ) {
    
    blosum_background_distr = c(0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072)
    
    align_file = infile_name
    align_suffix_elements = strsplit(align_file, "\\.")
    align_suffix = align_suffix_elements[[1]][1]
    s_matrix = read_scoring_matrix(s_matrix_file)
    names = c()
    alignment = c()
    
    if (is.defined(seq_weights) == TRUE) {
        seq_weights = c()
    }
    
    ali_out  = read_fasta_alignment(align_file)
    names = ali_out[[1]]; alignment = ali_out[[2]]
    
    seq_len = nchar(alignment[[1]][1])
    for (i in length(alignment[[1]])) {
        if (nchar(alignment[[i]][1]) != seq_len) {
            sprintf("ERROR: Sequences of different lengths: %s (%d) != %s (%d).\n", names[[1]], seq_len, names[[i]], nchar(seq))
            quit()
        }
    }
    
    if (is.defined(seq_weights) == TRUE) {
        seq_weights = load_sequence_weights(replace_element(align_file, align_suffix, '.weights'))
        if (seq_weights == c()) {
            seq_weights = calculate_sequence_weights(alignment)
        }
        if (length(seq_weights) != length(alignment)) { seq_weights = rep(1, length(alignment)) }
    }
    else { seq_weights = rep(1, length(alignment)) }
    
    
    # handle print of output relative to specific sequence
    ref_seq_num = "None"
    if (is.defined(seq_specific_output) == TRUE) {
        if (!(seq_specific_output %in% names[[1]])) {
            sprintf("Sequence %s not found in alignment. Using default output format...\n", seq_specific_output)
            seq_specific_output = 0
        }
    }
    
    else if (seq_specific_output %in% names[[1]]) {
        ref_seq_num = which(seq_specific_output %in% names[[1]])
    }
    
    # calculate scores
    scores = c()
    
    
    for (i in 1:nchar(alignment[[1]][1])) {
        col = get_column(i, alignment)
        if (length(col) == length(alignment)) {
            if (gap_percentage(col) <= gap_cutoff) {
                scores = c(scores, scoring_function(col, s_matrix, bg_distribution, seq_weights, use_gap_penalty))
            }
            else {
                scores = c(scores, -1000)
            }
        }
    }
    
    if (window_size > 0) {
        scores = window_score(scores, window_size, win_lam)
    }
    
    if (normalize_scores) {
        scores = calc_z_scores(scores, -999)
    }
    
    return(scores)
}













