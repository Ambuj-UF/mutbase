###############################################################################
# Ambuj Kumar, University of Florida
# Mrinal Mishra, University of Turku, Finland
###############################################################################
library(seqinr)
library(annotate)
blastSeqKK <- function (x, database = "nr", hitListSize = "100",filter = "L",expect = "100",program = "blastn",attempts = 10) {
  baseUrl <- "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
  query <- paste("QUERY=", as.character(x), "&DATABASE=", database,"&HITLIST_SIZE=", hitListSize, "&FILTER=", filter, "&EXPECT=",expect,"&PROGRAM=", program, sep = "")
  url0 <- sprintf("%s?%s&CMD=Put", baseUrl, query)
  results <- tempfile()
  Sys.sleep(5)
  require(XML)
  post <- htmlTreeParse(url0, useInternalNodes = TRUE)
  x <- post[["string(//comment()[contains(., \"QBlastInfoBegin\")])"]]
  rid <- sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
  rtoe <- as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1", 
                         x))
  url1 <- sprintf("%s?RID=%s&FORMAT_TYPE=XML&CMD=Get", baseUrl, 
                  rid)
  Sys.sleep(rtoe)
  .tryParseResult <- function(url, attempts){
    for (i in 1:(attempts+1)) {
      result <- tryCatch({
        xmlTreeParse(url, useInternalNodes=TRUE,
                     error = xmlErrorCumulator(immediate=FALSE))
      }, error=function(err) NULL)
      if (!is.null(result)) return(result)
      Sys.sleep(10)
    }
    stop(paste("no results after ", attempts, 
               " attempts; please try again later", sep = ""))
  }
  result <- .tryParseResult(url1, attempts)
  qseq <- xpathApply(result, "//Hsp_qseq", xmlValue)
  hseq <- xpathApply(result, "//Hsp_hseq", xmlValue)
  require(Biostrings)
  res <- list()
  for (i in seq_len(length(qseq))) {
    res[i] <- DNAMultipleAlignment(c(hseq[[i]], qseq[[i]]), 
                                   rowmask = as(IRanges(), "NormalIRanges"), colmask = as(IRanges(), 
                                                                                          "NormalIRanges"))
  }
  res
}


mySeq <- read.fasta("check.fas")

res <- lapply(mySeq, function(x){
  # collapse seq into string
  seqCollapse <- paste(toupper(as.character(x)), collapse = "")
  # run blast
  blastRes <- blastSeqKK(x = seqCollapse, database = "nr", hitListSize = 100,
                         attempts = 20)
  return(blastRes)
})

res
