################################################################################################################
#                                                                                                              #
# Copyright (C) {2014}  {Ambuj Kumar, Kimball-Brain lab group, Biology Department, University of Florida}      #
#                                                                                                              #
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

# Calculates coevlving amino acid sites

library(Biostrings)

Blossum <- function() {
    data_Matrix = c(
               4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2,  0,
               0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2,  0,
              -2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3,  0,
              -1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2,  0,
              -2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3,  0,
               0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3,  0,
              -2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2,  0,
              -1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1,  0,
              -1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2,  0,
              -1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1,  0,
              -1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1,  0,
              -2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2,  0,
              -1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3,  0,
              -1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1,  0,
              -1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2,  0,
               1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2,  0,
               0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2,  0,
               0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1,  0,
              -3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2,  0,
              -2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7,  0,
               0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
               )

    Blossum_Matrix = matrix(data_Matrix, nrow=21, ncol=21)
    nameVector = c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-')
    
    colnames(Blossum_Matrix) = nameVector
    rownames(Blossum_Matrix) = nameVector
    
    return(Blossum_Matrix)
}


Poisson_dist <- function (seq1, seq2) {
    dist = 0
    gap = 0
    for (i in 1:length(seq1)){
        if (seq1[[i]] == seq2[[i]] && seq1[[i]] != '-' && seq2[[i]] != '-') {
            dist = dist + 1
        }
        if (seq1[[i]] == '-' || seq2[[i]] == '-') {
            gap = gap + 1
        }
    }
    return(Poisson(dist, length(seq1) - gap))
}

Poisson <- function (distance, long) {
    pdist = 0
    if ((1 - (distance/long)) != 0) {
        pdist = -log(1 - (distance/long))
    }
    else { pdist = 0}
    
    return(pdist)
}



Relative_distance <- function (distance, sorted_dist) {
    for (i in 1:length(names(distance))) {
        if (sorted_dist[0] != 0) {
            distance[names(distance[i])] = distance[names(distance[i])]/sorted_dist[0]
        }
    }
    
    return(distance)
}



optimize <- function (seqObj) {
    distance = list()
    seqNameVec = c()
    for (i in 1:length(names(seqObj))) {
        seqNameVec = c(seqNameVec, names(seqObj[i]))
        inSeqObj = seqObj[!names(seqObj[i]) %in% seqNameVec]
        seq1 = seqObj[i][names(seqObj[i])]
        for (j in 1:length(names(seqNameVec))) {
            seq2 = seqObj[j][names(seqObj[j])]
            distance[[paste(names(seqObj[i]), names(seqObj[j]), sep = '-')]] = Poisson_dist(seq1, seq2)
        }
    }
    
    dist_val = c()
    for (i in 1:length(names(distance))) {
        dist_val = c(dist_val, distance[names(distance[i])])
    }
    
    sorted_dist = sort(dist_val)
    optimize_dist = relative_distance(distance, sorted_dist)
    
    return(optimize_dist)
}



#thetaCfunc <- function (optimize_dist) {
#    valTheta = c()
#    for (i in 1:length(names(optimize_dist))) {
#        valTheta = c(valTheta, optimize_dist[names(optimize_dist[i])])
#    }
#
#    return(mean(valTheta))
#}


thetaEK <- function (aliObj, optimize_dist) {
    bmat = Blossum()
    posList = list()
    for (i in 1:length(aliObj[names(aliObj[1])])) {
        storeName = c()
        data = c()
        for (m in 1:length(names(aliObj))) {
            storeName = c(storeName, names(aliObj[m]))
            inAliObj = aliObj[!names(aliObj) %in% storeName]
            for (n in 1:length(names(inAliObj))) {
                data = c(data, bmat[sapply(aliObj[[names(aliObj)[m]]], as.character), sapply(inAliObj[[names(inAliObj)[n]]],
                    as.character)]/optimize_dist[[paste(names(aliObj[m], names(inAliObj[n]), sep='-'))]])
            }
        }
        posList[[paste(pos, i, sep = '')]] = data
    }
    
    return(posList)
}


thetaC <- function (thetaEKVals) {
    allTheta = c()
    for (i in 1:length(names(thetaEKVals))) {
        for (theta in thetaEKVals[names(thetaEKVals)[i]]) {
            allTheta = c(allTheta, theta)
        }
    }
    return(mean(allTheta))
}


variability <- function (thetaEKVals, meanThetaEK) {
    for (i in 1:length(names(thetaEKVals))) {
        for (i in thetaEKVals[names(thetaEKVals)[i]]) {
            thetaEKVals[[names(thetaEKVals)[i]]][j] = (thetaEKVals[[names(thetaEKVals)[i]]][j] - meanThetaEK)^2
        }
    }
    
    return(thetaEKVals)
}


meanVariability <- function(variabilityData) {
    allVar = c()
    for (i in 1:length(names(variabilityData))) {
        for (var in variabilityData[names(variabilityData)[i]]) {
            allVar = c(allVar, var)
        }
    }
    
    return(mean(allVar))
}


correlation(variability, mvar) {
    corData = list()
    storeName = c()
    for (i in length(names(variability))) {
        storeName = c(storeName, names(variability[i]))
        inVariability = variability[!names(variability) %in% storeName]
        for (j in length(names(inVariability))) {
            corData[[paste(names(variability)[i], names(inVariability)[j]), sep='-']] =
                cor(variability[names(variability)[i]], inVariability[names(inVariability)[j]])
        }
    }
    
    return(corData)
}



corAllSite <- function (corData) {
    allCor = c()
    for (i in 1:length(names(corData))) {
        allCor = c(allCor, corData[names(corData)[i]])
    }
    
    return(c(mean(allCor), var(allCor)))
}






