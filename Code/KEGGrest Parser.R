Reactions <- read.delim("~/PotatoRecon/Data/Kegg_reaction.lst",header = TRUE)
Reactions <- as.vector(Reactions$keggList..reaction..)
reaction <- Reactions[64]

extract <- function(reaction){
  if (grepl(";",reaction)){
    parts <- unlist(strsplit(reaction,";[[:blank:]]+"))
    return(parts[length(parts)])
  }else{
    return(unlist(strsplit(reaction,"rn:R[[:digit:]]+[[:blank:]]"))[2])
  }
}
ids <- regmatches(Reactions,regexpr("R[[:digit:]]+",Reactions))
rx<-as.vector(sapply(Reactions, extract))
out <- cbind(ids,rx)
write.table(out,quote = FALSE,col.names = FALSE,row.names = FALSE,file = "KEGG.out")

