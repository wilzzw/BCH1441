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

#There are functions and variables (e.g. philDB) that I need from a state right before this question
#So I source BIN-Storing_data.R Line 1-269...
#EXCEPT: I commented out Line 135-136 because those are intentionally problematic for demonstration purposes

#source("1-269_BIN-Storing_data.R")

#Rewrite/overwrite autoincrement() with adjustable increment
#If default, it is the same as the autoincrement() before overwriting
autoincrement <- function(table, incr=1) {
  return(max(table$id) + incr)
}

#Adding new person and give it an ID by increment
newPersonID <- autoincrement(philDB$person) 
#Make a data frame object with the new philosopher's information
newPhil <- data.frame(id = newPersonID,
                      name = "Immanuel Kant",
                      born = "1724",
                      died = "1804",
                      school = "Enlightenment",
                      stringsAsFactors = FALSE)
#Append the new data frame to the person table in philDB
philDB$person <- rbind(philDB$person, newPhil)

#Adding two new books and give them IDs by increment
newBooksID <- c(autoincrement(philDB$books), autoincrement(philDB$books, 2))
#Make a data frame object with the new books' information
newBooks <- data.frame(id = newBooksID,
                       title = c("Critique of Pure Reason", "Critique of Judgement"),
                       published = c("1781", "1790"),
                       stringsAsFactors = FALSE)
#Append the new data frame to the book table in philDB
philDB$books <- rbind(philDB$books, newBooks)

#Create new junction (works) table entries
#And then append them to the works table in philDB
for (i in seq_along(newBooksID)) {
    workAssign <- data.frame(id = autoincrement(philDB$works), personID = newPersonID, bookID = newBooksID[i])
    philDB$works <- rbind(philDB$works, workAssign)
}

#Sort book titles alphabetically and rearrange books table accordingly in philDB
alphabetSort <- order(philDB$books$title)
df_sortedBooks <- philDB$books[alphabetSort,]

#Get ready for each book, gather their information by cross-referencing and print them out
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

#1.1  Task - fetchUniProtSeq() function
fetchUniProtSeq <- function(UniProtID) {
    URL <- sprintf("http://www.uniprot.org/uniprot/%s.fasta", UniProtID)
    response <- httr::GET(URL)
    if (httr::status_code(response) == 200) {
        sequence <- dbSanitizeSequence(strsplit(as.character(response), "\n"))
    } else {
        sequence <- character(0)
    }
    return(sequence)
}

#Task - fetchScanProsite() function
fetchScanProsite <- function(UniProtID) {
    URL <- "https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi"
    #idk whether all other options will be the same for all requests.. Fingers crossed.
    response <- httr::POST(URL,
                           body = list(meta = "opt1",
                                       meta1_protein = "opt1",
                                       seq = UniProtID,
                                       skip = "on",
                                       output = "tabular"))
    if (httr::status_code(response) == 200) {
        lines <- unlist(strsplit(httr::content(response, "text"), "\\n"))
        patt <- sprintf("\\|%s\\|", UniProtID)
        lines <- lines[grep(patt, lines)]
        features <- data.frame()
        for (line in lines) {
        tokens <- unlist(strsplit(line, "\\t|\\|"))
        features <- rbind(features,
                    data.frame(uID   =  tokens[2],
                               start =  as.numeric(tokens[4]),
                               end   =  as.numeric(tokens[5]),
                               psID  =  tokens[6],
                               psName = tokens[7],
                               stringsAsFactors = FALSE))
}
    } else {
        features <- character(0)
    }
    return(features)
}

#Task - fetchNCBItaxData() function
fetchNCBItaxData <- function(refSeqID) {
    eUtilsBase <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    #URL to get GI number
    URL <- paste(eUtilsBase,
             "esearch.fcgi?", 
             "db=protein",    
             "&term=", refSeqID,
             sep="")
    myXML <- xml2::read_xml(URL)

    #I have to error check a bit differently than the ones above, as the request structure changed
    #source("RPR-eUtils_XML.R") #Dependency on the node2text() function in the sourced code
    GID <- node2text(myXML, "Id")
    URL <- paste0(eUtilsBase,
            "esummary.fcgi?",
            "db=protein",
            "&id=",
            GID,
            "&version=2.0")
    myXML <- xml2::read_xml(URL)
    taxID <- node2text(myXML, "TaxId")
    organism <- node2text(myXML, "Organism")

    #Turns out to meet the requirement, I do not need to check for error.
    #If the retrieval is unsuccessful, node2text() anything will return character(0)
    #So that works for me here.
    outputList <- list(taxID, organism)
    names(outputList) <- c("TaxId", "Organism")
    
    return(outputList)
}

#Oral exam code
N <- 100000

GENCODE <- Biostrings::GENETIC_CODE #A named vector whose names are codons and elements are one-letter-code of corresponding amino acids
NUCLEOTIDES <- c("A", "G", "C", "T")
MUTATION_TYPES <- c("synonymous", "missense", "nonsense")

#mRNA sequences
load(file = "./data/ABC-INT-Mutation_impact.RData") #Get mRNA in terms of +strand DNA sequence data for 3 genes

#Translate codons into one-letter code of amino acids
translateCodon <- function(codon) {
    aa <- GENCODE[codon]
    return(aa)
}

#Skipped checking the case where before == after
typeOfMutation <- function(before, after) {
    aaBefore <- translateCodon(before)
    aaAfter <- translateCodon(after)
    if (aaBefore == aaAfter) {
        return("synonymous")
    } else if (aaAfter == "*") {
        return("nonsense")
    } else {
        return("missense")
    }
}

pointMutate <- function(nucleoSeq) {
    codonMutate <- sample(1:length(nucleoSeq), 1)                   #Find a codon position to mutate
    posMutate <- sample(1:3, 1)                                     #Find a position within the selected codon to mutate
    beforeCodon <- nucleoSeq[codonMutate]                           #Save the codon before mutation
    splitCodon <- unlist(strsplit(beforeCodon, split=""))           #Split nucleotide letters for substitution. unlist() turns it into a vector of characters
    beforeNucleotide <- splitCodon[posMutate]                       #Record the nucleotide to be mutated
    mutatePossible <- NUCLEOTIDES[NUCLEOTIDES != beforeNucleotide]  #Make sure mutation happens by excluding the same nucleotide as before mutation
    afterNucleotide <- sample(mutatePossible, 1)                    #Randomly mutate into ...
    splitCodon[posMutate] <- afterNucleotide                        #Mutate!
    afterCodon <- paste0(splitCodon, collapse="")                   #concat into codon after mutation
    #print(beforeCodon)
    #print(afterCodon)
    mutType <- typeOfMutation(beforeCodon, afterCodon)              #Categorize the type of mutation it just underwent
    return(mutType)

    #nucleoSeq[codonMutate] <- afterCodon
    #return(nucleoSeq)
}

#Start Mutating and collecting results!
resultKRAS <- character(0)
resultPTPN11 <- character(0)
resultOR1A1 <- character(0)

set.seed(1000)
for (i in 1:N) {
    resultKRAS[i] <- pointMutate(KRascodons)
    resultPTPN11[i] <- pointMutate(PTPN11codons)
    resultOR1A1[i] <- pointMutate(OR1A1codons)
}
set.seed(NULL)

#Calculate the frequency of each type of mutation
freqMut <- function(result, mutType) {
    numMutType <- length(result[result == mutType])
    proportionMutType <- numMutType / N
    return(proportionMutType)
}

#Collect frequencies
freqDistributions <- data.frame(KRas = numeric(length(MUTATION_TYPES)),
                                PTPN11 = numeric(length(MUTATION_TYPES)),
                                OR1A1 = numeric(length(MUTATION_TYPES)),
                                row.names = MUTATION_TYPES)

#Frequency data collected from IntOGen
#Exclude any splice mutations
databaseFreqDistributions <- data.frame(KRas = c(20/2310, 2287/2310, 3/2310),
                                        PTPN11 = c(33/248, 210/248, 5/248),
                                        OR1A1 = c(64/194, 121/194, 9/194),
                                        row.names = MUTATION_TYPES)

cat("KRas mutations:")
cat("\n")
for (m in MUTATION_TYPES) {
    freq <- freqMut(resultKRAS, m)
    freqDistributions[m, "KRas"] <- freq
    cat(sprintf("%s: %s", m, freq))
    cat("\n")
}

cat("PTPN11 mutations:")
cat("\n")
for (m in MUTATION_TYPES) {
    freq <- freqMut(resultPTPN11, m)
    freqDistributions[m, "PTPN11"] <- freq
    cat(sprintf("%s: %s", m, freqMut(resultPTPN11, m)))
    cat("\n")
}

cat("OR1A1 mutations:")
cat("\n")
for (m in MUTATION_TYPES) {
    freq <- freqMut(resultOR1A1, m)
    freqDistributions[m, "OR1A1"] <- freq
    cat(sprintf("%s: %s", m, freqMut(resultOR1A1, m)))
    cat("\n")
}

#Plot
#Method adopted from Sven Hohenstein's answer at https://stackoverflow.com/questions/20349929/stacked-bar-plot-in-r
barplot(as.matrix(freqDistributions))
barplot(as.matrix(databaseFreqDistributions))


#Compute KL-divergence. Without considering the lengths of p and q
KL_Div <- function(p, q) {
    val <- sum(p * log(p / q))
    return(val)
}

obsDiv <- numeric(0)

#For KRas...
(obsDiv["KRas"] <- KL_Div(freqDistributions[,"KRas"], databaseFreqDistributions[,"KRas"])) #0.6303814

#For PTPN11...
(obsDiv["PTPN11"] <- KL_Div(freqDistributions[,"PTPN11"], databaseFreqDistributions[,"PTPN11"])) #0.0399703

#For OR1A1...
(obsDiv["OR1A1"] <- KL_Div(freqDistributions[,"OR1A1"], databaseFreqDistributions[,"OR1A1"])) #0.02687701

#Simulate P-values
#First, I will simulate the distribution of the KL-Divergence
gene <- "KRas"
N <- 1000
nSample <- 100
truePMF <- freqDistributions[,gene]

divs <- numeric(N)
for (i in 1:N) {
  x <- table(sample(MUTATION_TYPES, nSample, replace = TRUE, prob=truePMF))
  divs[i] <- KLdiv(truePMF, pmfPC(x, MUTATION_TYPES))
}

#Now, how many are as extreme as the value, or more extreme?
(moreExtreme <- sum(divs >= obsDiv[gene]))
#P-value
(Pval <- moreExtreme / length(divs))

###Get SGD features into a data.frame 19-11-30

# This line takes a very long time to done correctly
# At first I assumed "SGD_features.tab" is formatted poorly, so I use fill = TRUE
# I could not get the nrow(SGD_features) correctly, but opening as csv works properly.
# I pinned down the first instance of problematic line to be Line 75 in "SGD_features.tab"
# There is a single quote ' in the description of Line 75, where "3'" is mentioned
# "SGD_features.README.txt" had some quote issues as well. In ?read.table, looks like the argument quote solves the problem
# fill = TRUE is not necessary. Turns out the format is done nicely (I guess...) for parsing 
SGD_features <- read.table("./data/SGD_features.tab", sep = "\t", quote = "", stringsAsFactors=FALSE)

# Name its columns according to "SGD_features.README.txt"
colnames(SGD_features) <- c("Primary SGDID", 
                            "Feature type", 
                            "Feature qualifier", 
                            "Feature name", 
                            "Standard gene name", 
                            "Alias", 
                            "Parent feature name",
                            "Secondary SGDID",
                            "Chromosome",
                            "Start_coordinate",
                            "Stop_coordinate",
                            "Strand",
                            "Genetic position",
                            "Coordinate version",
                            "Sequence version",
                            "Description")

#Select essential information and rename to be consistent with provided source code
subset <- c(1:5, 16)
SGD_features <- SGD_features[,subset]
colnames(SGD_features) <- c("SGDID",
                            "type",
                            "qual",
                            "sysName",
                            "name",
                            "description")

# Remove all rows that don't have a systematic name
hasSysName <- as.logical(nchar(SGD_features[, "sysName"]))
SGD_features <- SGD_features[hasSysName,]

# All systematic names are unique?
(! any(duplicated(SGD_features[,"sysName"])))

# Name rows after systematic names
rownames(SGD_features) <- SGD_features[,"sysName"]

# How many "genes" in the expression dataset is not in SGD features (yeast genome)?
notInFeatures <- setdiff(rownames(Biobase::exprs(GSE3635)), rownames(SGD_features))
(length(notInFeatures))
# Above is not actually the number of genes not in yeast genome. Some of them were actually names of controls
# The number of yeast reading frame that is not in SGD features database is:
(length(notInFeatures[-1:-11])) #63

# How many / which genes in the feature table do not have expression data?
noData <- setdiff(rownames(SGD_features), rownames(Biobase::exprs(GSE3635))) 
#1907 entries. Though not all of them are genes though e.g. centromere. 
# grep() according to SGD_features[,"type"] can be helpful for that purpose

#TASK>  Plot expression profiles for these genes and study them. What do you
#TASK>  expect the profiles to look like, given the role of these genes? What
#TASK>  do you find? (Hint: you will make your life much easier if you define
#TASK>  a function that plots and prints descriptions with a gene name as input.
#TASK>  Also: are the gene names in the feature table upper case, or lower case?
#TASK>  Also: note that the absolute values of change are quite different.
#TASK>  Also: note that some values may be outliers i.e. failed experiments.)

plotGene <- function(gName, color) {
    #Make sure it is uppercase
    gName <- toupper(gName)
    # Get index in SGD_features where the input gene occurs
    iFeature <- which(SGD_features$name == gName)
    # Get index in expression data GSE3635 where the input gene occurs
    iExprs <- which(Biobase::featureNames(GSE3635) == SGD_features$sysName[iFeature])

    # Plot it
    plot(seq(0, 120, by = 10),
         Biobase::exprs(GSE3635)[iExprs, ],
         main = "Expression profile",
         xlab = "time (min)",
         ylab = "expression",
         ylim = range(c(-2,2)),
         type = "b",
         col = color)

    # Description
    cat(sprintf("Gene Description: %s", SGD_features$description[iFeature]))
}

getColors <- function(numColors) {
    # Divide into intervals, whose values will be assigned to colors
    unitVals <- 1 / numColors
    comp1 <- seq(1, 0, length.out=numColors)
    comp2 <- seq(0, 1, length.out=numColors)

    colorOut <- c()
    for (i in 1:numColors) {
        colorOut[i] <- rgb(comp1[i], comp2[i], 0) # This is just the red mode. I am not familiar with how to select blue modes elegently in R
    }
    return(colorOut)
}

plotGeneSeries <- function(gNames) {
    allColors <- getColors(length(gNames))
    for (i in seq_along(gNames)) {
        color <- allColors[i]
        plotGene(gNames[i], color)
        # https://stackoverflow.com/questions/2564258/plot-two-graphs-in-same-plot-in-r
        par(new=TRUE)
    }
    legend("topleft", gNames, fill = allColors, bty = "n")
    abline(h =  0, col = "#00000055")
    abline(v = 60, col = "#00000055") 
}

onGenes <- c('Cdc14', 'Mbp1', 'Swi6', 'Swi4', 'Whi5', 'Cdc28', 'Cln1', 'Cln2', 'Cln3')
plotGeneSeries(onGenes)

offGenes <- c('Rad53', 'Cdc28', 'Clb1', 'Clb2', 'Clb6', 'Nrm1')
plotGeneSeries(offGenes)

hkGenes <- c('Act1', 'Alg9')
plotGeneSeries(hkGenes)


# ==   5.1  Final task: Gene descriptions  =====================================

#    Print the descriptions of the top ten differentially expressed genes
#    and comment on what they have in common (or not).

#Top 10 differentially expressed identified by R-package limma tools
topTenSysName <- myTable[1:10 , "ID"]
SGD_descriptions <- SGD_features$description
names(SGD_descriptions) <- SGD_features$sysName

#Retrieve descriptions by names
topDescriptions <- SGD_descriptions[topTenSysName]

for (i in seq_along(topDescriptions)) {
    print(topDescriptions[i])
}




#CPdefs <- list()
for (ID in ENSPsel) {
	cat(ID)
    cat("\n")
    infoID <- biomaRt::getBM(filters = "ensembl_peptide_id",
                             attributes = c("hgnc_symbol",
                                                  "wikigene_description",
                                                  "phenotype_description"),
                             values = ID,
                             mart = myMart)
	cat(sprintf("hgnc_symbol: %s", infoID[1,"hgnc_symbol"]))
    cat("\n")
    cat(sprintf("wikigene_description: %s", infoID[1,"wikigene_description"]))
    cat("\n")
    cat(sprintf("phenotype_description: %s", infoID[1,"phenotype_description"]))
    cat("\n")
    cat("\n")
}


### Phylogeny

# I just realized in "makeProteinDB.R", by including Line 9:
# myDB <- dbAddProtein(myDB, fromJSON("../ERECY_APSES_PSI-BLAST.json"))
# Sequence duplicates are added to myDB that are in both "ERECY_APSES_PSI-BLAST.json" and "./data/refAPSES_PSI-BLAST.json"
# I don't know whether having them will affect my multiple sequence alignment. Or automatically removes 100% duplicates
# So I tried the following modified prep script as well. The difference is removing the duplicates from mySeq

# Write and save fasta MSA format first
withDups <- writeMFA(APSESmsa)

library(msa)
 
# Align all sequences in the database + KILA_ESSCO
mySeq <- myDB$protein$sequence
names(mySeq) <- myDB$protein$name
mySeq <- c(mySeq,
           "IDGEIIHLRAKDGYINATSMCRTAGKLLSDYTRLKTTQEFFDELSRDMGIPISELIQSFKGGRPENQGTWVHPDIAINLAQ")
names(mySeq)[length(mySeq)] <- "KILA_ESCCO"

# Remove duplicates!!!!
mySeq <- mySeq[! duplicated(mySeq)]
# Check
any(duplicated(mySeq))
any(duplicated(names(mySeq)))


mySeqMSA <- msaClustalOmega(AAStringSet(mySeq)) # too many sequences for MUSCLE
 
# get the sequence of the SACCE APSES domain
sel <- myDB$protein$name == "MBP1_SACCE"
proID <- myDB$protein$ID[sel]
 
sel <- myDB$feature$ID[myDB$feature$name == "APSES fold"]
fanID <- myDB$annotation$ID[myDB$annotation$proteinID == proID &
                            myDB$annotation$featureID == sel]
start <- myDB$annotation$start[fanID]
end   <- myDB$annotation$end[fanID]
 
SACCEapses <- substring(myDB$protein$sequence[proID], start, end)
 
# extract the APSES domains from the MSA
APSESmsa <- fetchMSAmotif(mySeqMSA, SACCEapses)

# Write fasta MSA again
withoutDups <- writeMFA(APSESmsa)

# The alignments look different!!!
# From here I have to use APSESmsa reflected in withoutDups, which has duplicates removed

# Filter out E. coli KilA-N and all Mbp1 orthologs 
not.ECOLI.or.MBP1 <- ! (grepl("KILA_ESCCO", APSESmsa@ranges@NAMES) | grepl("^MBP1_", APSESmsa@ranges@NAMES))

#Randomly sample 10 sequences that are not E.coli KilA or Mbp1 orthologs
set.seed(112358)
APSES.seqSample <- sample(APSESmsa[not.ECOLI.or.MBP1], 10)
set.seed(NULL)

# I can combine sampled sequence s with E.coli KilA and Mbp1 orthologs like this:
APSES.fullSeqSample <- c(APSESmsa[! not.ECOLI.or.MBP1], APSES.seqSample)


# Masking >80% gap columns.
# The following code is largely inspired by 4.1 section in "BIN-PHYLO-Data_preparation.R"
lengthAli <- nchar(APSES.fullSeqSample[1])
msaMatrix <- matrix(character(length(APSES.fullSeqSample) * lengthAli), ncol = lengthAli)

# Sanity check. Should print TRUE
(nrow(msaMatrix) == length(APSES.fullSeqSample))

# Put sequences into the matrix with individual characters
for (i in 1:nrow(msaMatrix)) {
    msaMatrix[i,] <- seqinr::s2c(as.character(APSES.fullSeqSample[[i]]))
}

# Create a masking vector
colMask <- logical(ncol(msaMatrix))
limit <- round(nrow(msaMatrix) * (80/100))

for (i in 1:ncol(msaMatrix)) {
    count <- sum(msaMatrix[ , i] == "-") # Number of characters in the column that would be "-"
    colMask[i] <- (count <= limit) # TRUE if fewer than 80% are gaps
}

# Masking! Also preparing input
maskedMatrix <- msaMatrix[ , colMask]

phyloInput <- character()
for (i in 1:nrow(maskedMatrix)) {
    phyloInput[i] <- paste(maskedMatrix[i, ], collapse="")
}
names(phyloInput) <- APSES.fullSeqSample@ranges@NAMES

# Perform MSA
sampledSeqMSA <- msaClustalOmega(Biostrings::AAStringSet(phyloInput))