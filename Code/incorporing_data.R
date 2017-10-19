##Anotación de ESTs
library(gdata)
library(seqinr)
library(DESeq)
library(UniProt.ws)
#Leer datos
Exp_set <- read.xls("~/PotatoRecon/Data/exp_set.xlsx", sheet = 1)
dim(Exp_set)[1]
# Creando multifasta para BLAST
write.fasta(sequences = strsplit(as.vector(Exp_set[1:dim(Exp_set)[1],1]),""),
            names = 1:dim(Exp_set)[1],
            file.out = "~/PotatoRecon/Data/exp_set.fasta")
#BLAST
#Creación de bases de datos
##Fuente: Index of ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000226075.1_SolTub_3.0/
#makeblastdb -in GCF_000226075.1_SolTub_3.0_protein.faa -parse_seqids -dbtype prot # Con esta DB no se obtinen anotación probablemente por la corta longitud del query
#makeblastdb -in GCF_000226075.1_SolTub_3.0_rna.fna -parse_seqids -dbtype nucl # Final
#makeblastdb -in GCF_000226075.1_SolTub_3.0_genomic.fna -parse_seqids -dbtype nucl  # Esta DB no sirve para anotar genes porque los ID son de Scafolds completos

# BLASTN con cobertura 100% y E-Value menor a 0.1 (final)
#Se realiza en el servidor de la universidad
# nohup blastn -db ../DBs/rna_potato/GCF_000226075.1_SolTub_3.0_rna.fna -query exp_set.fasta -task 'blastn-short' -perc_identity 100 -num_alignments 1 -outfmt 6 -evalue 0.1 -out anotation_ESTs.out  -num_threads 16

#BLASTX (test)
#nohup blastx -db ../DBs/protein_potato/GCF_000226075.1_SolTub_3.0_protein.faa -query exp_set.fasta -outfmt 6 -evalue 0.0000001 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -max_target_seqs 1 -out annotation_ESTs_prot.out -num_threads 12
#No hay hit para blastx, ni siquiera dejando los parametros menos restrictivos

# BLASTN con cobertura 100% y E-Value menor a 0.1 (test)
#blastn -db ../DBs/genomic_potato/GCF_000226075.1_SolTub_3.0_genomic.fna -query exp_set.fasta -task 'blastn-short' -perc_identity 100 -num_alignments 1 -outfmt 6 -evalue 0.1 -out anotation_ESTs_genomic.out  -num_threads 12 

#Cargar anotaciones
Ano_exp_set<-read.csv("~/PotatoRecon/Data/annotation_ESTs.out",sep = "\t",header = FALSE)

# Buscando genes únicos
names<-as.vector(Ano_exp_set[,2])
names
length(names)
unique_names <-names(table(Ano_exp_set[,2]))
length(unique_names)

# Declarando función para obtener el número máximo de conteos
colMax <- function(data){sapply(data, max, na.rm = TRUE)}

# Función para extraer el número máximo de conteos desde nuestros datos
maxcounts<- function(id){
  colMax(Exp_set[Ano_exp_set[Ano_exp_set[,2]==id,1],2:10])
}

# Extrayendo datos finales
Exp_set<-t(sapply(unique_names[grep("_",unique_names)], maxcounts))
head(Exp_set)

# Escribiendo archivo de resultados
setwd("~/PotatoRecon/Results/")
write.table(x = Exp_set, quote = FALSE, sep = "\t", file = "Expr_set_anotado")

##Obtener GENEID para los genes anotados
#Manualmente se modifica "Expr_set_anotado"
# 1. Se pone titulo a la columna de los genes anotados
# 2. Se obtiene el promedio de los conteos por cada tratamiento
# 3. Se guarda como xlsx

Exp_set <- read.xls("~/PotatoRecon/Results/Expr_set_anotado.xlsx", sheet = 1)

sotUp <- UniProt.ws(taxId=4113)
keytypes(sotUp)
columns(sotUp)
egs = keys(sotUp, "ENTREZ_GENE")
ref = keys(sotUp, "REFSEQ_NUCLEOTIDE")

gen_id<-select(sotUp, keys = as.vector(Exp_set$Genes), columns = "ENTREZ_GENE", keytype = "REFSEQ_NUCLEOTIDE")
write.table(x = gen_id, file = "id_entrez-gene",quote = FALSE, sep = "\t")
#En este punto en excel se eliminan los IDs duplicados y se combinan las tablas (Expr_set_anotado y id_entrez_gene) para obtener los IDs Entrez Gene a los conteos

###Incorporacion de datos al modelo
source("https://bioconductor.org/biocLite.R")
biocLite("Biobase")
library(devtools)
install_github("dosorio/minval")
library(minval)
library(Biobase)
library(gdata)
library(sybil)
library(sybilSBML)
source("https://bioconductor.org/biocLite.R")
biocLite("gage")
library(gage)
install_github("gibbslab/exp2flux")
library(exp2flux)
#Leer datos 
data_Exp <- read.xls("~/Dropbox/model_analysis/Results/Expr_set_anotado_inc.xlsx")
ID_Gen <- read.csv("~/Dropbox/model_analysis/Data/ID_ENTREZ_GENE.csv",sep = "\t",header = TRUE)
#ID_ENTREZ_GENE.csv es un archivo con solo los ENTREZ_GENE de Expr_set_anotado_inc.xlsx

##Sumarización de ENTREZ_GENE
# Buscando genes únicos
namesEG<-as.vector(ID_Gen[,2])
namesEG
length(namesEG)
unique_namesEG <-names(table(ID_Gen[,2]))
length(unique_namesEG)

# Declarando función para obtener el número máximo de conteos
colMax <- function(data){sapply(data, max, na.rm = TRUE)}

# Función para extraer el número máximo de conteos desde nuestros datos
maxcounts<- function(id){
  colMax(data_Exp[ID_Gen[ID_Gen[,2]==id,1],2:13])
}

# Extrayendo datos finales
Data_Exp<-t(sapply(unique_namesEG, maxcounts))
head(Data_Exp)
dim(Data_Exp)
setwd("~/PotatoRecon/Results/")
write.table(x = Data_Exp, quote = FALSE, sep = "\t", file = "Exp_set_inc")
#En este punto se organizan los titulos del dataset manualmente

##Creando los expresiónSet
data <- read.xls("~/PotatoRecon/Results/Exp_set_inc.xlsx")
data <- unique(data[!is.na(data[,"ENTREZ_GENE"]),])
Expr_set0<- ExpressionSet(matrix(data[,5],dimnames = list(data[,"ENTREZ_GENE"],c())))
Expr_set24<- ExpressionSet(matrix(data[,9],dimnames = list(data[,"ENTREZ_GENE"],c())))
Expr_set72<- ExpressionSet(matrix(data[,13],dimnames = list(data[,"ENTREZ_GENE"],c())))

##Generar y cargar modelo SBML
setwd("~/PotatoRecon/Results/")
potato <- read.csv("Potato3.6-3.csv")
potato <- as.data.frame.array(potato)
colnames(potato) <- c("ID","NAME","REACTION","EC","GPR","LOWER.BOUND","UPPER.BOUND","OBJECTIVE","COMPARTMENT","SUBSYSTEM")
writeSBMLmod(modelData = potato,modelID = "Potato 3.6.3",outputFile = "Potato.xml")
po <- readSBMLmod("Potato.xml")
po <- addExchReact(po,"Biomass[c]",lb = 0,ub = 1000)
optimizeProb(po)
writeSBMLmod(modelData = po,outputFile = "Potato.xml")
optimizeProb(readSBMLmod("Potato.xml"))
#Generar modelo con Cobrapy (Code: Potato_v.3.6.ipynb)
#Cargar modelo generado con minval 
#model1 <- readSBMLmod("Potato.xml")
#Cargar modelo generado con cobrapy
model1 <- readSBMLmod("Potato3.6xml")
model1 <- addExchReact(model1,"Biomass[c]",ub = 1000)
findExchReact(model1)
optimizeProb(model1)
model1@lowbnd[75]
FBA_modelin <- optimizeProb(model1, algorithm = "fba", retOptSol = FALSE)

#numero de reacciones con flujo.
length(FBA_modelin$fldind[FBA_modelin$fluxes!=0])
rxn_flux_in<- FBA_modelin$fldind[FBA_modelin$fluxes!=0]
write.table(rxn_flux_in,"rxn_flux_in", sep = "\t")

FBA_modelinit <- optimizeProb(model1, algorithm = "fba", retOptSol = TRUE)

####MEAN###
### Estos son los modelos determinados para la interacción P. infestans y s. tuberosum
#Incorporar Datos de expresión al modelo tomando  el promedio de los datos de expresión como valor de flujo para las reacciones que no tienen datos de expresión
#par(mfcol=c(1,3))
par(lab=c(5,20,10))
model0h<- exp2flux(model = model1, expression = Expr_set0, organism = "sot", missing = "mean")
FBA_model0h <- optimizeProb(model0h)
plot(getFluxDist(FBA_model0h),cex=0.1,ylab="Flux 0 h",xlab="Reaction",main=round(FBA_model0h@lp_obj,digits = 5))
table(getFluxDist(FBA_model0h))
FBA_model0h <- optimizeProb(model0h, algorithm = "fba", retOptSol = FALSE)
length(FBA_model0h$fldind[FBA_model0h$fluxes!=0])
rxn_flux0<- FBA_model0h$fldind[FBA_model0h$fluxes!=0]
write.table(rxn_flux0, file = "rxn_flux0", sep = "\t")

model24h<- exp2flux(model = model1, expression = Expr_set24, organism = "sot", missing = "mean")
FBA_model24h <- optimizeProb(model24h)
plot(getFluxDist(FBA_model24h),cex=0.1,ylab="Flux 24 h",xlab="Reaction",main=round(FBA_model24h@lp_obj,digits = 5))
FBA_model24h <- optimizeProb(model24h, algorithm = "fba", retOptSol = FALSE)
length(FBA_model24h$fldind[FBA_model24h$fluxes!=0])
rxn_flux24<-FBA_model24h$fldind[FBA_model24h$fluxes!=0]
write.table(rxn_flux24, file = "rxn_flux24", sep = "\t")

model72h<- exp2flux(model = model1, expression = Expr_set72, organism = "sot", missing = "mean")
  #par(lab=c(2,10,1))
FBA_model72h <- optimizeProb(model72h)
plot(getFluxDist(FBA_model72h),cex=0.1,ylab="Flux 72 h",xlab="Reaction",main=round(FBA_model72h@lp_obj,digits = 5))
FBA_model72h <- optimizeProb(model72h, algorithm = "fba", retOptSol = FALSE)
length(FBA_model72h$fldind[FBA_model72h$fluxes!=0])
rxn_flux72<-FBA_model72h$fldind[FBA_model72h$fluxes!=0]
write.table(rxn_flux72, file = "rxn_flux72", sep = "\t")

#Numero de reacciones con flujo en todos los modelos
rxn_flux <- read.xls(xls = "~/Dropbox/model_analysis/Results/rxn_flux_general_model.xlsx")
length(unique(rxn_flux$X))

###Results analysis
##Obtener flujo de reacciones claves modelo metabolico inicial sin resticciones valores de expresión
FBA_modelin$fluxes[2046]#RBS01
FBA_modelin$fluxes[74]#RK0004
FBA_modelin$fluxes[75]#RK0005
FBA_modelin$fluxes[71]#R00024
FBA_modelin$fluxes[285]#R01512
FBA_modelin$fluxes[282]#R01061
FBA_modelin$fluxes[283]#R01063
FBA_modelin$fluxes[281]#R01015
FBA_modelin$fluxes[345]#R01068
FBA_modelin$fluxes[344]#R00762
FBA_modelin$fluxes[770]#R01067
FBA_modelin$fluxes[72]#R01829
FBA_modelin$fluxes[73]#R01845
FBA_modelin$fluxes[795]#R01641
FBA_modelin$fluxes[399]#R01056
FBA_modelin$fluxes[793]#R01523
FBA_modelin$fluxes[1622]#R03140
FBA_modelin$fluxes[1623]#R01334
FBA_modelin$fluxes[1624]#R00475
FBA_modelin$fluxes[1625]#R00009
FBA_modelin$fluxes[1626]#R00372
FBA_modelin$fluxes[312]#R00945
FBA_modelin$fluxes[1627]#R01221
FBA_modelin$fluxes[313]#R00588
FBA_modelin$fluxes[315]#R01388
FBA_modelin$fluxes[314]#R01514
FBA_modelin$fluxes[1886]#R02110
FBA_modelin$fluxes[1866]#R05196
FBA_modelin$fluxes[1865]#R00028
FBA_modelin$fluxes[400]#R01529

###Obtener flujo de reacciones claves modelo metabolico 0hdp
FBA_model0h$fluxes[2046]#RBS01
FBA_model0h$fluxes[74]#RK0004
FBA_model0h$fluxes[75]#RK0005
FBA_model0h$fluxes[71]#R00024
FBA_model0h$fluxes[285]#R01512
FBA_model0h$fluxes[282]#R01061
FBA_model0h$fluxes[283]#R01063
FBA_model0h$fluxes[281]#R01015
FBA_model0h$fluxes[345]#R01068
FBA_model0h$fluxes[344]#R00762
FBA_model0h$fluxes[770]#R01067
FBA_model0h$fluxes[72]#R01829
FBA_model0h$fluxes[73]#R01845
FBA_model0h$fluxes[795]#R01641
FBA_model0h$fluxes[399]#R01056
FBA_model0h$fluxes[793]#R01523
FBA_model0h$fluxes[1622]#R03140
FBA_model0h$fluxes[1623]#R01334
FBA_model0h$fluxes[1624]#R00475
FBA_model0h$fluxes[1625]#R00009
FBA_model0h$fluxes[1626]#R00372
FBA_model0h$fluxes[312]#R00945
FBA_model0h$fluxes[1627]#R01221
FBA_model0h$fluxes[313]#R00588
FBA_model0h$fluxes[315]#R01388
FBA_model0h$fluxes[314]#R01514
FBA_model0h$fluxes[1878]#R02110
FBA_model0h$fluxes[1866]#R05196
FBA_model0h$fluxes[1865]#R00028
FBA_model0h$fluxes[400]#R01529
FBA_model0h$fluxes[398]#R01051
FBA_model0h$fluxes[1973]#RK0003
FBA_model0h$fluxes[566]#R00948

###Obtener flujo de reacciones claves modelo metabolico 2hdp
FBA_model24h$fluxes[2046]#RBS01
FBA_model24h$fluxes[74]#RK0004
FBA_model24h$fluxes[75]#RK0005
FBA_model24h$fluxes[71]#R00024
FBA_model24h$fluxes[285]#R01512
FBA_model24h$fluxes[282]#R01061
FBA_model24h$fluxes[283]#R01063
FBA_model24h$fluxes[281]#R01015
FBA_model24h$fluxes[345]#R01068
FBA_model24h$fluxes[344]#R00762
FBA_model24h$fluxes[770]#R01067
FBA_model24h$fluxes[72]#R01829
FBA_model24h$fluxes[73]#R01845
FBA_model24h$fluxes[795]#R01641
FBA_model24h$fluxes[399]#R01056
FBA_model24h$fluxes[793]#R01523
FBA_model24h$fluxes[1622]#R03140
FBA_model24h$fluxes[1623]#R01334
FBA_model24h$fluxes[1624]#R00475
FBA_model24h$fluxes[1625]#R00009
FBA_model24h$fluxes[1626]#R00372
FBA_model24h$fluxes[312]#R00945
FBA_model24h$fluxes[1627]#R01221
FBA_model24h$fluxes[313]#R00588
FBA_model24h$fluxes[315]#R01388
FBA_model24h$fluxes[314]#R01514
FBA_model24h$fluxes[1878]#R02110
FBA_model24h$fluxes[1866]#R05196
FBA_model24h$fluxes[1865]#R00028
FBA_model24h$fluxes[400]#R01529
FBA_model24h$fluxes[398]#R01051
FBA_model24h$fluxes[1973]#R01051
FBA_model24h$fluxes[566]#R00948

###Obtener flujo de reacciones claves modelo metabolico 3hdp
FBA_model72h$fluxes[2046]#RBS01
FBA_model72h$fluxes[74]#RK0004
FBA_model72h$fluxes[75]#RK0005
FBA_model72h$fluxes[71]#R00024
FBA_model72h$fluxes[285]#R01512
FBA_model72h$fluxes[282]#R01061
FBA_model72h$fluxes[283]#R01063
FBA_model72h$fluxes[281]#R01015
FBA_model72h$fluxes[345]#R01068
FBA_model72h$fluxes[344]#R00762
FBA_model72h$fluxes[770]#R01067
FBA_model72h$fluxes[72]#R01829
FBA_model72h$fluxes[73]#R01845
FBA_model72h$fluxes[795]#R01641
FBA_model72h$fluxes[399]#R01056
FBA_model72h$fluxes[793]#R01523
FBA_model72h$fluxes[1622]#R03140
FBA_model72h$fluxes[1623]#R01334
FBA_model72h$fluxes[1624]#R00475
FBA_model72h$fluxes[1625]#R00009
FBA_model72h$fluxes[1626]#R00372
FBA_model72h$fluxes[312]#R00945
FBA_model72h$fluxes[1627]#R01221
FBA_model72h$fluxes[313]#R00588
FBA_model72h$fluxes[315]#R01388
FBA_model72h$fluxes[314]#R01514
FBA_model72h$fluxes[1878]#R02110
FBA_model72h$fluxes[1866]#R05196
FBA_model72h$fluxes[1865]#R00028
FBA_model72h$fluxes[400]#R01529
FBA_model72h$fluxes[398]#R01051
FBA_model72h$fluxes[398]#R01051
FBA_model72h$fluxes[1973]#RK0003
FBA_model72h$fluxes[566]#R00948
















