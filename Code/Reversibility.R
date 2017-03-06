###Curar reversibilidad

# MetaCyc
dGMetaCyc <- read.csv("~/PotatoRecon/Data/All_reactions_of_MetaCyc.csv")
dGMetaCyc <- dGMetaCyc[!dGMetaCyc$KEGG=="",]
dGMetaCyc <- dGMetaCyc[!is.na(dGMetaCyc$Gibbs.0),]
id <- regmatches(dGMetaCyc$KEGG,regexpr("R[[:digit:]]+",dGMetaCyc$KEGG))
dG_MetaCyc <- dGMetaCyc$Gibbs.0
dGMetaCyc <- as.data.frame(cbind(id,dG_MetaCyc))

# Equilibrator
EquilibratorpH7 <- read.csv("~/PotatoRecon/Data/kegg_reactions_CC_ph7.0.csv")
id <- as.vector(EquilibratorpH7$X.MiriamID..urn.miriam.kegg.reaction)
dG_Equilibrator <- EquilibratorpH7$X.dG0_prime..kJ.mol.
dGEquilibrator <- as.data.frame(cbind(id,dG_Equilibrator))
dGEquilibrator <- dGEquilibrator[!is.na(dGEquilibrator$dG_Equilibrator),]

# Merge
dG <- as.data.frame.array(merge(dGEquilibrator,dGMetaCyc,by.x = "id",by.y = "id"))
dG[,4] <- NA
dG[as.numeric(dG[,2])< 0 & as.numeric(dG[,3])<0, 4] <-"=>"
dG[as.numeric(dG[,2])>=(-0.9) & as.numeric(dG[,3])>=(-0.9),4] <-"<=>"
dG[is.na(dG[,4]),4] <- "<=>"
colnames(dG)[4] <- "Direction"


# Cargando modelo
library(gdata)
Potato <- read.xls("~/PotatoRecon/Data/RXN-reversibilidad.xlsx",sheet = 3)
#Potato <- Potato[!Potato$X.=="#",]
sot <- as.data.frame(cbind(as.vector(Potato$ID),as.vector(Potato$EQUATION)))
colnames(sot) <- c("id","reaction")


change.reversibility <- function(id, reaction){
  original <- regmatches(reaction,regexpr("<?=>?",reaction))
  if (id%in%dG$id){
  new <- dG[dG$id%in%id,4]
  reaction <-gsub(original,new,reaction)
    }
  return(reaction)
}

newreact <- NULL
for(i in 1:dim(sot)[1]){
  newreact<-c(newreact,change.reversibility(as.vector(sot$id)[i],as.vector(sot$reaction)[i]))
}

newreact


Potato$EQUATION <- newreact

# Cambio del Flujo dependiendo de la flecha

Potato$LOWER.BOUND[grepl(" => ",Potato$EQUATION)] <- 0
Potato$UPPER.BOUND[grepl(" => ",Potato$EQUATION)] <- 1000
Potato$LOWER.BOUND[grepl(" <=> ",Potato$EQUATION)] <- -1000
Potato$UPPER.BOUND[grepl(" <=> ",Potato$EQUATION)] <- 1000

setwd("~/PotatoRecon/Results/")
write.csv2(Potato,file = "potato_reversibility")


