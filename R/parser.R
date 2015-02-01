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




trim <- function (x) gsub("^\\s+|\\s+$", "", x)
trail <- function (x) paste("</", gsub("<", "", x), sep="")
trim_flag <- function (x, flagger) gsub(trail(flagger), "", gsub(flagger, "", x))

blast_parser <- function(filename) {
    text <- readLines(filename,encoding="UTF-8")
    object_vector <- list()
    
    hitFlag = 0
    for (lines in text) {
        if (isTRUE(grepl("<Hit_id>", lines))) {
            hitFlag <- 1
            id = trim(trim_flag(lines, "<Hit_id>"))
            object_vector[[id]] <- list()
        }
        
        if (hitFlag == 1) {
            if (isTRUE(grepl("<Hsp_evalue>", lines))) {
                object_vector[[id]]["eval"] = trim(trim_flag(lines, "<Hsp_evalue>"))
            }
            if (isTRUE(grepl("<Hsp_hit-from>", lines))) {
                object_vector[[id]]["start"] = trim(trim_flag(lines, "<Hsp_hit-from>"))
            }
            if (isTRUE(grepl("<Hsp_hit-to>", lines))) {
                object_vector[[id]]["stop"] = trim(trim_flag(lines, "<Hsp_hit-to>"))
            }
            if (isTRUE(grepl("<Hsp_hit-frame>", lines))) {
                object_vector[[id]]["frame"] = trim(trim_flag(lines, "<Hsp_hit-frame>"))
            }
            if (isTRUE(grepl("<Hsp_hseq>", lines))) {
                object_vector[[id]]["seq"] = trim(trim_flag(lines, "<Hsp_hseq>"))
                hitFlag = 0
            }
        }
    }
    
    return(object_vector)
}
