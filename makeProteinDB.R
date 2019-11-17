# This command needs to be executed whenever you recreate the database. 
# In particular, whenver you have added or modified data in any of the JSON files.

source("./scripts/ABC-createRefDB.R")

myDB <- dbAddProtein(myDB, fromJSON("../MBP1_ERECY.json"))
myDB <- dbAddTaxonomy(myDB, fromJSON("../MYSPEtaxonomy.json"))
myDB <- dbAddAnnotation(myDB, fromJSON("../ERECY-Annotations.json")) # Added for the sake of activity 4 in "BIN-ALI-Optimal_sequence_alignment.R"

#Alternative code
#source("C:/Users/Wilson/Documents/BCH1441/ABC-units/scripts/ABC-createRefDB.R")
#setwd("C:/Users/Wilson/Documents/BCH1441")

#myDB <- dbAddProtein(myDB, fromJSON("MBP1_ERECY.json"))
#myDB <- dbAddTaxonomy(myDB, fromJSON("MYSPEtaxonomy.json"))