###############################################################################
# Mrinal Mishra, University of Turku Finland
###############################################################################
##Multiple seq aligner for Blast output 
## Sequences with their respective ID parsed from Blast output and aligned
library(muscle)
library(seqinr)
align_seq<- function(filename){
  parsed_data=blast_parser(filename)
  seq=c()
  seqname=c()
  for (i in 1:length(parsed_data)){
    seq <- c(seq, parsed_data[[i]]$seq)
    seqname <- c(seqname, names(parsed_data)[i])
     }
  seqinr::write.fasta(sequences = as.list(seq), names = as.list(seqname), nbchar = 80, file.out = "All_seq.fasta",open="w")
  muscle(seqs = "All_seq.fasta", out = "align.fas")
  }
