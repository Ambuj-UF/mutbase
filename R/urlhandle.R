################################################################################################################
#                                                                                                              #
# Copyright (C) {2014}  {Ambuj Kumar, Kimball-Brain lab group, Biology Department, University of Florida}      #
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




require(RCurl)
require(XML)


urlencode <- function(object_list){
    if (is.null(object_list[["db"]])) {cat("Database name required")}
    if (is.null(object_list[["id"]])) {object_list[[-2]]}
    if (is.null(object_list[["seq_start"]])) {object_list[[-3]]}
    if (is.null(object_list[["seq_stop"]])) {object_list[[-4]]}
    if (is.null(object_list[["strand"]])) {object_list[[-5]]}
    s="?"
    for (key in names(object_list)) {
        s = paste(s, key, "=", object_list[[key]], "&", sep="")
    }
    
    return(s)
    
}


urlopen <- function(object_list){
    cgi <- 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    variables <- urlencode(object_list)
    url_submit <- paste(cgi, variables, "rettype=fasta&retmode=text", sep="")
    print(url_submit)
    webpage <- getURL(url_submit)
    webpage <- readLines(tc <- textConnection(webpage))
    close(tc)
    return(webpage)

}


object_list <- vector(mode="list", length=5)
names(object_list) <- c("db", "id", "strand", "seq_start", "seq_stop")


#object_list[["db"]] = "nuccore"
#object_list[["id"]] = 568815597
#object_list[["strand"]] = 1
#object_list[["seq_start"]] = 197084323
#object_list[["seq_stop"]] = 197084425

#result=urlopen(object_list)
#sink("mydata.fas")
#cat(result[1:3][1],"\n",result[1:3][2], "\n")
#sink()


