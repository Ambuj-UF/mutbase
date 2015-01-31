source("Blast.R")

out<-capture.output(res)
cat(out,file="out.txt",sep="\n",append=TRUE)

