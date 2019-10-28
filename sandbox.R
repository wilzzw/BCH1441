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


#Line 123: database_name(List_name)$entity_name(Dataframe_name)$attribute_name(Column_name)[ID(Row_name)]
philDB$person$name[1]

# task: Write an expression that returns all "school" entries from the
#       person table.
philDB$person$school

#Line 136: will throw an error:
#Error in match.names(clabs, names(xi)) : 
#  names do not match previous names

#Make sure to turn of stringsAsFactor

#Line 264: If not cat() but rather print() is used, it will print out the vector object in its barebone representation..

#1.3 Task
#source the file from Line 1 to Line 269
autoincrement <- function(table, incr=1) {
  return(max(table$id) + incr)
}

newPersonID <- autoincrement(philDB$person)
newPhil <- data.frame(id = newPersonID,
                      name = "Immanuel Kant",
                      born = "1724",
                      died = "1804",
                      school = "Enlightenment",
                      stringsAsFactors = FALSE)
philDB$person <- rbind(philDB$person, newPhil)

newBooksID <- c(autoincrement(philDB$books), autoincrement(philDB$books, 2))
newBooks <- data.frame(id = newBooksID,
                       title = c("Critique of Pure Reason", "Critique of Judgement"),
                       published = c("1781", "1790"),
                       stringsAsFactors = FALSE)
philDB$books <- rbind(philDB$books, newBooks)

for (i in seq_along(newBooksID)) {
    workAssign <- data.frame(id = autoincrement(philDB$works), personID = newPersonID, bookID = newBooksID[i])
    philDB$works <- rbind(philDB$works, workAssign)
}


alphabetSort <- order(philDB$books$title)
df_sortedBooks <- philDB$books[alphabetSort,]

for (i in seq_along(df_sortedBooks$id)) {
    bookid <- df_sortedBooks$id[i]
    title <- df_sortedBooks$title[i]
    yearPublished <- df_sortedBooks$published[i]

    selectAssign <- philDB$works$bookID == bookid
    personid <- philDB$works$personID[selectAssign]
    selectPhil <- philDB$person$id == personid
    author <- philDB$person$name[selectPhil]

    cat(sprintf("\"%s\" - ", title))
    cat(sprintf("%s", author))
    cat(sprintf(" (%s)", yearPublished))
    cat("\n")
}

#Line 426: 
#source("./scripts/ABC-createRefDB.R")
# Error in fromJSON("./data/MBP1_SACCE.json") : 
#  could not find function "fromJSON"
#Solution: library(jsonlite)

#Line 449: all(myDB$protein$taxonomyID %in% myDB$taxonomy$ID)
#all() returns TRUE if all elements of the vector argument are TRUE

#Modifying data is done directly in the JSON files

lengthToLookNice <- 70

MYSEQ <- "        1 msstsvasrd qiysakysgv evyefihptg simkrkaddw vnathilkaa kfakakrtri
       61 lekevikdih ekvqggfgky qgtwvpldia rrlaekfdvl eelrplfdfs qrdgsasppq
      121 apkhhhasrs dstkkgtgks psgalnnasg svipkrrgrp prskkldrip asgdaalqrs
      181 rsdvtgfhkp sitistissh nlpsiqstlq rgvnideqqh yqdklqqqis qqkyeeldie
      241 dglssdietn layiaegpvg snrlntqlmt gkeepvssss slpsspsdfs apvpfdtqrv
      301 gsatspigam lprysmqsrp ptsdldqkvn dylaklvdyf insemqntna vpiellnpph
      361 stpyidawid sehhtafhwa camgnlpive allqagashr avnhagetpl mrasmfhnsy
      421 tkrtyprifq llqdtvfdid sqsqtvvhhi vkrrsntqsa lyyldvllsk ikdfspqyri
      481 etlintqddk gntplhiaai ngdkkffqtl lgsgalstlk nydgvtadvf innkfsrtln
      541 ysehsyhygn gtthspasts tgavitgpag aaaasasasf ihtgdmfpsq aatsvsraip
      601 evinlmkdma dsyqglyqdr sqelqsikkm lksmnntvas vdikiletld ikkyeqigqt
      661 meditqaide lqsrftvkqk clmnilekgq riqlqrline qeqeidkhqe esesksgpsi
      721 npnlitgike lailqlrrka kikqmlellc gnskvqkfrk misqgtdmel devdnfldvi
      781 lqqlnddnea kkinnpngvt"

cleanSeq <- dbSanitizeSequence(MYSEQ)

(seqLength <- nchar(cleanSeq))

toStore <- character()
i <- 0
while (i*lengthToLookNice+1 <= seqLength) {
    toStore[i+1] <- substr(cleanSeq, start=i*lengthToLookNice+1, stop=(i+1)*lengthToLookNice)
    i <- i + 1
}

toStore

#Line 566: Another functionality of cbind() I was not aware of: supplying extra arguments to create new columns. Great for calculated/derived properties from columns

#Line 589-590: Failed verification
#> sel <- myDB$taxonomy$species == MYSPE
#> myDB$taxonomy[sel, ]
#[1] ID      species
#<0 rows> (or 0-length row.names)

#Solution: changed "Eremothecium cymbalariae DBVPG#7215" to "Eremothecium cymbalariae" in MYSPEtaxonomy.json. I'm not sure what DBVPG#7215 indicates, but I'm removing it.
#Then, source("../makeProteinDB.R") again

#Task 3.4
# - On your submission page, note the E-value of your protein and link
#     to its NCBI protein database page.

#E-value rounded to 0.0
#NCBI page: https://www.ncbi.nlm.nih.gov/protein/XP_003645298.1?report=genbank&log$=protalign&blast_rank=1&RID=UT6SVB5U014

# - Copy and paste the contents of your two JSON files on your submission
#     page on the Student Wiki
MBP1_ERECY.json:
[
  { "name" : "MBP1_ERECY",
    "RefSeqID" : "XP_003645298",
    "UniProtID" : "G8JQ18",
    "taxonomyID" : 931890,
    "sequence" : [
       "MSSTSVASRDQIYSAKYSGVEVYEFIHPTGSIMKRKADDWVNATHILKAAKFAKAKRTRILEKEVIKDIH",
       "EKVQGGFGKYQGTWVPLDIARRLAEKFDVLEELRPLFDFSQRDGSASPPQAPKHHHASRSDSTKKGTGKS",
       "PSGALNNASGSVIPKRRGRPPRSKKLDRIPASGDAALQRSRSDVTGFHKPSITISTISSHNLPSIQSTLQ",
       "RGVNIDEQQHYQDKLQQQISQQKYEELDIEDGLSSDIETNLAYIAEGPVGSNRLNTQLMTGKEEPVSSSS",
       "SLPSSPSDFSAPVPFDTQRVGSATSPIGAMLPRYSMQSRPPTSDLDQKVNDYLAKLVDYFINSEMQNTNA",
       "VPIELLNPPHSTPYIDAWIDSEHHTAFHWACAMGNLPIVEALLQAGASHRAVNHAGETPLMRASMFHNSY",
       "TKRTYPRIFQLLQDTVFDIDSQSQTVVHHIVKRRSNTQSALYYLDVLLSKIKDFSPQYRIETLINTQDDK",
       "GNTPLHIAAINGDKKFFQTLLGSGALSTLKNYDGVTADVFINNKFSRTLNYSEHSYHYGNGTTHSPASTS",
       "TGAVITGPAGAAAASASASFIHTGDMFPSQAATSVSRAIPEVINLMKDMADSYQGLYQDRSQELQSIKKM",
       "LKSMNNTVASVDIKILETLDIKKYEQIGQTMEDITQAIDELQSRFTVKQKCLMNILEKGQRIQLQRLINE",
       "QEQEIDKHQEESESKSGPSINPNLITGIKELAILQLRRKAKIKQMLELLCGNSKVQKFRKMISQGTDMEL",
       "DEVDNFLDVILQQLNDDNEAKKINNPNGVT"]
  }
]

MYSPEtaxonomy.json:
[
  { "ID" : 931890,
    "species" : "Eremothecium cymbalariae"}
]
#Should I add to the original taxonomy json or just one entry?

# - Execute the two commands below and show the result on your submission page
#> biCode(myDB$taxonomy$species) %in% biCode(MYSPE)
# [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE

#> myDB$protein$taxonomyID %in% myDB$taxonomy$ID[(myDB$taxonomy$species == MYSPE)]
# [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#[14] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#[27] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#[40] FALSE FALSE FALSE FALSE FALSE  TRUE



#p.d.f.


# Task: Write a function calcGC() that calculates GC content in a sequence.
#       Hint: you could strsplit() the sequence into a vector, and count
#       G's and C's; or you could use gsub("[AT]", "", <sequence>) to remove
#       A's and T's, and use nchar() before and after to calculate the content
#       from the length difference.
#       Then write tests that:
#          confirm that calcGC("AATT") is 0;
#          confirm that calcGC("ATGC") is 0.5;
#          confirm that calcGC("AC")   is 0.5;
#          confirm that calcGC("CGCG") is 1;

calcGC <- function(sequence) {
    splitSeq <- unlist(strsplit(sequence, split=""))
    totalLength <- length(splitSeq)
    sel <- splitSeq == "G" | splitSeq == "C"
    GCcontent <- length(splitSeq[sel])
    GC <- GCcontent / totalLength
    return(GC)
}

library(testthat)

expect_equal(calcGC("AATT"), 0)
expect_true(all.equal(calcGC("ATGC"), 0.5))
expect_equal(calcGC("AC"), 0.5)
expect_equal(calcGC("CGCG"), 1)

#Alternative genetic code task
myDNA <- readLines("./data/S288C_YDL056W_MBP1_coding.fsa")[-1]
myDNA <- paste0(myDNA, collapse = "")
myDNA <- as.character(Biostrings::codons(Biostrings::DNAString(myDNA)))
myDNA <- myDNA[-length(myDNA)]

stdCode <- Biostrings::GENETIC_CODE

N <- 200
valSWGC <- numeric(N)

set.seed(112358)
for (i in 1:N) {
    swGC <- swappedGC(stdCode)
    x <- traRev(myAA, swGC)
    x <- randMut(x)
    x <- traFor(x, swGC)
    valSWGC[i] <- evalMut(myAA, x)
}
set.seed(NULL)

hist(valSTDC,
     breaks = 15,
     col = "palegoldenrod",
     xlim = c(0, 400),
     ylim = c(0, N/4),
     main = "Standard vs. Swapped Genetic Code",
     xlab = "Mutation penalty")

hist(valSWGC,
     col = "plum",
     breaks = 60,
     add = TRUE)