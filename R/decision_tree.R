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


revOrder_list_key <- function(listObj) {
    """Reverse order list object by keys"""
    
    new_list = list()
    name_sorted_rev = sort(names(listObj), decreasing =TRUE)
    for (name in name_sorted_rev) {
        new_list[name] = listObj[[name]]
    }
    
    return(new_list)
}


majority_count <- function(classlist) {
    classcount = list()
    for (vote in classlist) {
        if !(vote %in% names(classcount)) {
            classcount[vote] = 0
        }
        
        classcount[vote] = classcount[[vote]] + 1
    }
    
    sortedClassCount = revOrder_list_key(classcount)
    
    return(sortedClassCount[[1]][[1]])
}


entropy <- function(dataset) {
    n = length(dataset)
    labels = list()
    for (record in dataset) {
        label = record[-length(record)]
        if !(label %in% names(labels)) {
            labels[label] = 0
        }
        labels[label] = labels[label] + 1
    }
    
    ent = 0
    
    for (key in names(labels)) {
        prob = labels[[key]] / n
        ent = -prob * log(prob, 2)
    }
    
    return(ent)
}















