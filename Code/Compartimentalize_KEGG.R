kegg_sot <- read.csv("~/PotatoRecon/Results/kegg_sot")
kegg_rx <- as.vector(kegg_sot$reaction)
kegg_rx <- as.vector(sapply(kegg_rx, function(reaction){
  reaction <- unlist(strsplit(reaction,"[[:blank:]]+<=>[[:blank:]]+"))
  reactants <- unlist(strsplit(reaction[1],"[[:blank:]]\\+[[:blank:]]"))
  products <- unlist(strsplit(reaction[2],"[[:blank:]]\\+[[:blank:]]"))
  reactants <- paste0(as.vector(sapply(reactants, function(metabolite){paste0(metabolite,"[c]")},simplify = TRUE)),collapse = " + ")
  products <- paste0(as.vector(sapply(products, function(metabolite){paste0(metabolite,"[c]")},simplify = TRUE)),collapse = " + ")
  return(paste(reactants,products,sep = " <=> "))
}))
kegg_sot <- cbind(as.vector(kegg_sot$id),kegg_rx)
setwd("~/PotatoRecon/Results/")
write.csv2(kegg_sot,file = "compartimentalizado")


