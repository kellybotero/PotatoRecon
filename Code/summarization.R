##Anotación de ESTs
library(gdata)
library(seqinr)
library(DESeq)
library(UniProt.ws)
#Leer datos
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
