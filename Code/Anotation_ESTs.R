##Anotación de ESTs
library(gdata)
library(seqinr)
library(DESeq)
library(UniProt.ws)
library(gdata)
Exp_set <- read.xls("~/Dropbox/model_analysis/Data/exp_set.xlsx", sheet = 1)
dim(Exp_set)[1]
# Creando multifasta para BLAST
write.fasta(sequences = strsplit(as.vector(Exp_set[1:dim(Exp_set)[1],1]),""),
            names = 1:dim(Exp_set)[1],
            file.out = "Dropbox/model_analysis/Data/exp_set.fasta")
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
Ano_exp_set<-read.csv("~/Dropbox/model_analysis/Data/annotation_ESTs.out",sep = "\t",header = FALSE)

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
setwd("~/Dropbox/model_analysis/Results/")
write.table(x = Exp_set, quote = FALSE, sep = "\t", file = "Expr_set_anotado")

##Obtener GENEID para los genes anotados
#Manualmente se modifica "Expr_set_anotado"
# 1. Se pone titulo a la columna de los genes anotados
# 2. Se obtiene el promedio de los conteos por cada tratamiento
# 3. Se guarda como xlsx

Exp_set <- read.xls("~/Dropbox/model_analysis/Results/Expr_set_anotado.xlsx", sheet = 1)

sotUp <- UniProt.ws(taxId=4113)
keytypes(sotUp)
columns(sotUp)
egs = keys(sotUp, "ENTREZ_GENE")
ref = keys(sotUp, "REFSEQ_NUCLEOTIDE")

gen_id<-select(sotUp, keys = as.vector(Exp_set$Genes), columns = "ENTREZ_GENE", keytype = "REFSEQ_NUCLEOTIDE")
write.table(x = gen_id, file = "id_entrez-gene",quote = FALSE, sep = "\t")
#En este punto en excel se eliminan los IDs duplicados y se combinan las tablas (Expr_set_anotado y id_entrez_gene) para obtener los IDs Entrez Gene a los conteos


#### Pruebas

namesRef <- gen_id$REFSEQ_NUCLEOTID
namesRef
length(namesRef)
unique_namesRef<- unique(gen_id$REFSEQ_NUCLEOTIDE)

length(unique_namesRef)
Exp_set$Genes[Exp_set$Genes%in%unique_namesRef] <- gen_id$ENTREZ_GENE[Exp_set$Genes%in%unique_namesRef]
for(i in Exp_set$Genes){
  print(i)
  Exp_set$Genes <-factor(gen_id$ENTREZ_GENE[Exp_set$Genes%in%unique_namesRef])
  Exp_set[Exp_set$Genes]<-as.factor(gen_id$ENTREZ_GENE[unique_namesRef%in%Exp_set$Genes])
  print(Exp_set$Genes)
}
 
Exp_set$Genes[Exp_set$Genes%in%unique_namesRef] <- factor(gen_id$ENTREZ_GENE[unique_namesRef%in%Exp_set$Genes])
Potato1$LOWER.BOUND[Potato1$ID%in%id_rev]=-1000

replace(x = as.vector(Exp_set$Genes), list = as.vector(gen_id$ENTREZ_GENE), values = gen_id$ENTREZ_GENE)
  

replace(x = as.vector(gen_id$ENTREZ_GENE[unique_namesRef%in%Exp_set$Genes]), list = as.vector(Exp_set$Genes[Exp_set$Genes%in%unique_namesRef]), values = gen_id$ENTREZ_GENE)

for(i in Exp_set$Genes){
  grep(pattern = i, x = gen_id$REFSEQ_NUCLEOTIDE)
}


ids_ref<- as.vector(unique(gen_id$REFSEQ_NUCLEOTIDE))

ids_EG<-(unique(gen_id[gen_id[,1]%in%ids_ref,2]))

ids_EG<-(unique(gen_id[gen_id[,1]%in%Exp_set$Genes,2]))

ids_EG

ids_EG<-paste0(gen_id[gen_id[,1]%in%ids_ref,2],collapse = " ")

gen_id_unique<- as.data.frame(cbind(ids_ref,ids_EG))

summary.ec <- function(id) {
  data <- reaction_all[reaction_all[, "id"] %in% id, ]
  data[, "ec"] <- paste0(unique(data[, "ec"]), collapse = sep)
  ec <- as.vector(sapply(id,function(id){paste0(sot[sot[,1]%in%id,3],collapse = "   ")}))
  
  sel <- gen_id[-c(duplicated(gen_id$REFSEQ_NUCLEOTIDE)),]
  sel <- gen_id[ gen_id$REFSEQ_NUCLEOTIDE ,]

###Sumarizar y normalizar datos de expresión
#Obtener caracteristicas de los datos

#dim(Exp_set_anotado)
#head(Exp_set_anotado)

#Crear el vector de condiciones
#cond <- factor(c("0h","0h","0h","24h","24h", "24h", "72h", "72h","72h"))
#cond

#7. Creando countsTable
#rownames(Exp_set_anotado) <- Exp_set_anotado$Genes
#Exp_set_anotado[,-1]
#countsTable <- Exp_set_anotado[,-1]
#head(countsTable)
#countsTable <- as.data.frame(Exp_set[,-1])
#head(countsTable)

#8.Creando objeto para usar DeSeq

#cds <- newCountDataSet(countsTable, cond)
#cds 

#9. Ajustar factores de tamaño y normalización 
#cds <- estimateSizeFactors(cds)
#sizeFactors(cds)
#norm_count<-(counts(cds,normalized=TRUE))

# obtener la media de los conteos por tratamiento
#col1 <-conditions(cds)=="DMSO"
#head(counts(cds)[,col1])
#col2 <-conditions(cds)== "CA"
#head(counts(cds)[,col2])
#means_count_t1<-((rowMeans(t(t(counts(cds)[,col1])/sizeFactors(cds)[col1]))))
#head(means_count_t1)
#means_count_t2<-((rowMeans(t(t(counts(cds)[,col2])/sizeFactors(cds)[col2]))))



