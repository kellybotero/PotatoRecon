library(gdata)
library(minval)
Potato3 <- read.csv("~/PotatoRecon/Results/Potato3.6-3.csv")
Potato3 <- Potato3[,c(1,2,3,4,5,6,7,8,9)]
names(Potato3) <- c("ID","NAME","REACTION","EC.NUMBER","GPR","LOWER.BOUND","UPPER.BOUND","OBJECTIVE","COMPARTMENT")
Potato3 <- as.data.frame.array(Potato3)

##REACTIONS
#asiganar una variable como un vector con todas las ecuciones de potato 
potato_equation<-as.vector(Potato3$REACTION)
potato_equation

##METABOLITOS
#extraer los metabolitos
metabolites(potato_equation)
#Quitar los compartimentos a metabolitos
metabolites(potato_equation,woCompartment = TRUE)

##EC.NUMBERS
ecs <- unique(unlist(strsplit(unique(as.vector(Potato3$EC.NUMBER)),"[[:space:]]",perl = TRUE)))
length(ecs)
#Conteo de ECs por actividad catalica
unlist(strsplit(Potato3$EC.NUMBER," "))
#conteo pora las que comienzan con 1
length(grep("^1",unlist(strsplit(Potato3$EC.NUMBER," "))))

##GENES
#Primero se debe cargar el modelo (UploadinPotatoR)
length(unique(model1@allGenes))

##GPR
gpr <- gsub("[()]","",Potato3$GPR)
gpr <- gsub("[[:space:]]","",Potato3$GPR)
length(grep(pattern = "[()]", gpr ))



###PATHWAYS
#Manual


