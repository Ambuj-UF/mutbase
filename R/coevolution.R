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


rangeF <- function(pos1, stLength, split_num) {
    retObj = c()
    for (i in 1:stLength) {
        if (i%%split_num == 1) {
            retObj = c(retObj, i)
        }
    }
    
    return(retObj)
}


spliter <- function(stObj, num) {
    retObj = c()
    for (start in rangeF(1, nchar(stObj), num)) {
        if (start == 1) {
            retObj = c(retObj, trim_seq(stObj, endPos=start+num))
        }
        else if (start == nchar(stObj) - 2) {
            retObj = c(retObj, trim_seq(stObj, initPos=start-1))
        }
        else {
            retObj = c(retObj, trim_seq(stObj, start-1, start+num))
        }
        
    }
    
    return(retObj)
}



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




hydrophob <- list(
    'F' = c(100, 92, 165.2),
    'I' = c(99, 100, 131.2),
    'W' = c(97, 84, 204.2),
    'L' = c(97, 100, 131.2),
    'V' = c(76, 79, 117.2),
    'M' = c(74, 74, 149.2),
    'Y' = c(63, 49, 181.2),
    'C' = c(49, 52, 121.2),
    'A' = c(41, 47, 89.1),
    'T' = c(13, 13, 119.1),
    'H' = c(8, -42, 155.2),
    'G' = c(0, 0, 75.1),
    'S' = c(-5, -7, 105.1),
    'Q' = c(-10, -18, 146.2),
    'R' = c(-14, -26, 174.2),
    'K' = c(-23, -37, 146.2),
    'N' = c(-28, -41, 132.1),
    'E' = c(-31, 8, 147.1),
    'P' = c(-46, -46, 115.1),
    'D' = c(-55, -18, 133.1),
    '-' = c(0, 0, 0)
)
    


z_table <- list(
    '0' =  c(0.5,       0.504,	0.508,	0.512,	0.516,	0.52,	0.5239,	0.5279,	0.5319,	0.5359),
    '1' =	c(0.5398,	0.5438,	0.5478,	0.5517,	0.5557,	0.5596,	0.5636,	0.5675,	0.5714,	0.5753),
    '2' =	c(0.5793,	0.5832,	0.5871,	0.591,	0.5948,	0.5987,	0.6026,	0.6064,	0.6103,	0.6141),
    '3' =	c(0.6179,	0.6217,	0.6255,	0.6293,	0.6331,	0.6368,	0.6406,	0.6443,	0.648,	0.6517),
    '4' =	c(0.6554,	0.6591,	0.6628,	0.6664,	0.67,	0.6736,	0.6772,	0.6808,	0.6844,	0.6879),
    '5' =	c(0.6915,	0.695,	0.6985,	0.7019,	0.7054,	0.7088,	0.7123,	0.7157,	0.719,	0.7224),
    '6' =	c(0.7257,	0.7291,	0.7324,	0.7357,	0.7389,	0.7422,	0.7454,	0.7486,	0.7517,	0.7549),
    '7' =	c(0.758,     0.7611,	0.7642,	0.7673,	0.7704,	0.7734,	0.7764,	0.7794,	0.7823,	0.7852),
    '8' =	c(0.7881,	0.791,	0.7939,	0.7967,	0.7995,	0.8023,	0.8051,	0.8078,	0.8106,	0.8133),
    '9' =	c(0.8159,	0.8186,	0.8212,	0.8238,	0.8264,	0.8289,	0.8315,	0.834,	0.8365,	0.8389),
    '10' =  c(0.8413,	0.8438,	0.8461,	0.8485,	0.8508,	0.8531,	0.8554,	0.8577,	0.8599,	0.8621),
    '11' =	c(0.8643,	0.8665,	0.8686,	0.8708,	0.8729,	0.8749,	0.877,	0.879,	0.881,	0.883),
    '12' =	c(0.8849,	0.8869,	0.8888,	0.8907,	0.8925,	0.8944,	0.8962,	0.898,	0.8997,	0.9015),
    '13' =	c(0.9032,	0.9049,	0.9066,	0.9082,	0.9099,	0.9115,	0.9131,	0.9147,	0.9162,	0.9177),
    '14' =	c(0.9192,	0.9207,	0.9222,	0.9236,	0.9251,	0.9265,	0.9279,	0.9292,	0.9306,	0.9319),
    '15' =	c(0.9332,	0.9345,	0.9357,	0.937,	0.9382,	0.9394,	0.9406,	0.9418,	0.9429,	0.9441),
    '16' =	c(0.9452,	0.9463,	0.9474,	0.9484,	0.9495,	0.9505,	0.9515,	0.9525,	0.9535,	0.9545),
    '17' =	c(0.9554,	0.9564,	0.9573,	0.9582,	0.9591,	0.9599,	0.9608,	0.9616,	0.9625,	0.9633),
    '18' =	c(0.9641,	0.9649,	0.9656,	0.9664,	0.9671,	0.9678,	0.9686,	0.9693,	0.9699,	0.9706),
    '19' =	c(0.9713,	0.9719,	0.9726,	0.9732,	0.9738,	0.9744,	0.975,	0.9756,	0.9761,	0.9767),
    '20' = c(0.9772,	0.9778,	0.9783,	0.9788,	0.9793,	0.9798,	0.9803,	0.9808,	0.9812,	0.9817),
    '21' =	c(0.9821,	0.9826,	0.983,	0.9834,	0.9838,	0.9842,	0.9846,	0.985,	0.9854,	0.9857),
    '22' =	c(0.9861,	0.9864,	0.9868,	0.9871,	0.9875,	0.9878,	0.9881,	0.9884,	0.9887,	0.989),
    '23' =	c(0.9893,	0.9896,	0.9898,	0.9901,	0.9904,	0.9906,	0.9909,	0.9911,	0.9913,	0.9916),
    '24' =	c(0.9918,	0.992,	0.9922,	0.9925,	0.9927,	0.9929,	0.9931,	0.9932,	0.9934,	0.9936),
    '25' =	c(0.9938,	0.994,	0.9941,	0.9943,	0.9945,	0.9946,	0.9948,	0.9949,	0.9951,	0.9952),
    '26' =	c(0.9953,	0.9955,	0.9956,	0.9957,	0.9959,	0.996,	0.9961,	0.9962,	0.9963,	0.9964),
    '27' =	c(0.9965,	0.9966,	0.9967,	0.9968,	0.9969,	0.997,	0.9971,	0.9972,	0.9973,	0.9974),
    '28' =	c(0.9974,	0.9975,	0.9976,	0.9977,	0.9977,	0.9978,	0.9979,	0.9979,	0.998,	0.9981),
    '29' =	c(0.9981,	0.9982,	0.9982,	0.9983,	0.9984,	0.9984,	0.9985,	0.9985,	0.9986,	0.9986),
    '30' = c(0.9987,	0.9987,	0.9987,	0.9988,	0.9988,	0.9989,	0.9989,	0.9989,	0.999,	0.999)
)










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
        seqNameVec = c(seqNameVec, names(seqObj[[i]]))
        inSeqObj = seqObj[!names(seqObj[i]) %in% seqNameVec]
        seq1 = seqObj[[i]][[names(seqObj[i])]]
        for (j in 1:length(names(seqNameVec))) {
            seq2 = seqObj[[j]][names(seqObj[j])]
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


meanTheta <- function (thetaEKVals) {
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


meanVar <- function(variabilityData) {
    allVar = c()
    for (i in 1:length(variabilityData)) {
        variabilityData[names(variabilityData)[[i]]] = mean(variabilityData[[names(variabilityData)[[i]]]])
    }
    
    return(variabilityData)
}


listAverage <- function(listObj) {
    data = c()
    for (i in 1: length(listObj)) {
        data = c(data, listObj[[names(listObj)[[i]]]])
    }
    
    return(mean(data))
}


thetaParam <- function(meanVarDict) {
    lengthAlign = length(meanVarDict)
    thetaParamVal = listAverage(meanVarDict)/length(meanVarDict)
    posCorExec = c()
    for (i in length(meanVarDict)) {
        if (meanVarDict[[names(meanVarDict)[[i]]]] > thetaParamVal) {
            posCorExec = c(posCorExec, names(meanVarDict)[[i]])
        }
    }
    
    return posCorExec
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







