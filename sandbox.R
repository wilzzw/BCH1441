#A sandbox. The extension helps my text editor (Visual Studio Code) to format the code.
txt <- ""
count_down <- 3
position <- 1
while (count_down >= 0) {
    if (count_down > 0) {
        txt[position] <- as.character(count_down)
        #print(as.character(count_down))
    } else {
        txt[position] <- "Lift Off!"
    }
    count_down <- count_down - 1
    position <- position + 1
}

print(txt)

countDown <- function(x=3) {
    txt <- ""
    start <- x
    position <- 1
    while (start >= 0) {
        if (start > 0) {
            txt[position] <- as.character(start)
        } else {
            txt[position] <- "Lift Off!"
        }
        start <- start - 1
        position <- position + 1
    }
    return(txt)
}

myLifeDays <- function(birthday, lday) {
    if (missing(birthday)) {
        print ("Enter your birthday as a string in \"YYYY-MM-DD\" format.")
        return()
    }
    bd <- strptime(birthday, "%Y-%m-%d") # convert string to time
    now <- format(Sys.time(), "%Y-%m-%d") # convert "now" to time
    diff <- round(as.numeric(difftime(now, bd, unit="days")))
    print(sprintf("This date was %d days ago.", diff))

    futureDay <- format(bd + as.difftime(lday, unit="days"), "%Y-%m-%d") #Consulted ?difftime
    print(sprintf("Celebrate your life day on: %s", futureDay))
    return()
}

#setwd("C:/Users/Wilson/Documents/BCH1441/ABC-units")
#source("FND-Genetic_code.R")

# Task: What do you need to change to print the table with U instead
#         of T? Try it.
for (i in seq_along(dim(cCube))) {
    dimnames(cCube)[[i]][4] <- "U"
}

# Task: Point mutations are more often transitions (purine -> purine;
#         pyrimidine -> pyrimidine) than transversions (purine -> pyrimidine;
#         pyrimidine -> purine), even though twice as many transversions
#         are possible in the code. This is most likely due a deamination /
#         tautomerization process that favours C -> T changes. If the code
#         indeed minimizes the effect of mutations, you would expect that
#         codons that differ by a transition code for more similar amino acids
#         than codons that differ by a transversion. Is that true? List the set
#         of all amino acid pairs that are encoded by codons with a C -> T
#         transition. Then list the set of amino acid pairs with a C -> A
#         transversion. Which set of pairs is more similar?

#I'm not sure whether I understand this question correctly but:
#All of them that contain "C"
codonWithC <- genCode[grep("C", names(genCode))]
codonVector <- genCode[1:length(genCode)] #For some reason...
#codonVector <- as.vector(genCode) #For some reason...
#names(codonVector) <- names(genCode)
#dictWith1C <- dictWithC[-grep("CC", names(dictWithC))]

mutate <- function(codon, from, to) {
    codonSplit <- strsplit(codon, split="")[[1]] #
    positionMutate <- which(codonSplit == from) #
    allMutated <- character(length(positionMutate))
    for (i in seq_along(positionMutate)) {
        mutated <- codonSplit
        mutated[positionMutate[i]] <- to
        mutated <- paste(mutated, sep="", collapse="")
        allMutated[i] <- mutated
    }
    return(allMutated)
}

translate <- function(codon) {
    return(codonVector[codon])
}

effectOfMutation <- function(codonBefore, from, to) {
    codonAfter <- mutate(codonBefore, from, to)
    aaBefore <- translate(codonBefore)[codonBefore]
    aaBefore <- as.character(aaBefore)
    aaAfter <- lapply(codonAfter, translate)#[codonAfter] #
    aaAfter <- as.character(aaAfter)
    mutationEffect <- data.frame(aaAfter, rownames=codonAfter)
    column <- paste(aaBefore, " (", as.character(codonBefore), ")", sep="", collapse="")
    #print(ncol(mutationEffect))
    colnames(mutationEffect) <- c(column, "Codon After")
    return(mutationEffect)
}

#CtoT <- function(codon) {
#    effect <- effectOfMutation(codon, "C", "T")
#    return(effect)
#}

#CtoA <- function(codon) {
#    effect <- effectOfMutation(codon, "C", "A")
#    return(effect)
#}

#CtoT_EffectTable <- CtoT(names(codonWithC)[1])
#for (i in seq_along(codonWithC[-1])) {
    #https://stackoverflow.com/questions/7739578/merge-data-frames-based-on-rownames-in-r
#    CtoT_EffectTable <- merge(CtoT_EffectTable, CtoT(names(codonWithC)[i]), all=TRUE) 
#}

outcomesOfMutation <- function(codonVector, from, to) {
    effectTable <- effectOfMutation(names(codonVector)[1], from, to)
    for (i in seq_along(codonVector[-1])) {
        #The use of merge() is inspired by https://stackoverflow.com/questions/7739578/merge-data-frames-based-on-rownames-in-r
        effectTable <- merge(effectTable, effectOfMutation(names(codonVector)[i], from, to), all=TRUE)
    }
    rownames(effectTable) <- effectTable[,"Codon After"]
    return(effectTable)
}

#The following nested for loops don't work to replace NA's as empty string. I'm not sure what's wrong yet.
for (aftermut in rownames(datf)) {
    for (beforemut in colnames(datf)) {
        if (is.na(datf[aftermut, beforemut])) {
            datf[aftermut, beforemut] <- ""
        }
    }
}

write.csv(outcomesOfMutation(codonWithC, "C", "T"), file="C_to_T_table.csv")
write.csv(outcomesOfMutation(codonWithC, "C", "A"), file="C_to_A_table.csv")
#And then I replaced all NA to empty in Excel. I don't know how to do this in R yet...

#dim-reduction Find subsets, and do statistics on the subset

'''
for (aftermut in rownames(outcomesOfMutation(codonWithC, "C", "T"))) {
    for (beforemut in colnames(outcomesOfMutation(codonWithC, "C", "T"))) {
        if (!is.na(outcomesOfMutation(codonWithC, "C", "T")[aftermut, beforemut])) {
            
        }
    }
}
'''