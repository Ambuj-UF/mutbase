text <- readLines("ResultCheck.104",encoding="UTF-8")

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
