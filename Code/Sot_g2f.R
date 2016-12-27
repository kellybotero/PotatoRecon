# Librerias
library(gdata)
library(devtools)
library(minval)
library(KEGGREST)
library(sybil)
library(g2f)


###VERIFICACIÓN DE RUTAS METABOLICAS DE KEGG EN POTATO_RECON
path_sot<-unique(keggList("pathway", "sot")) 
path_sot

#ingresar pathways potatorecon
path_potato <- read.csv("~/PotatoRecon/Data/path_potato.csv",header=FALSE)# se hace con la potatorecon_3.7
path_potato <- unique(paste0(gsub(",$","",paste(path_potato$V1,path_potato$V2,sep = ","))," - Solanum tuberosum (potato)"))
path_potato

# Paths from SOT in POTATOrecon
setwd("~/PotatoRecon/Results/")
path_sot_in_potato<-as.vector(path_sot[path_sot%in%path_potato])
path_sot_in_potato
write.csv(path_sot_in_potato, 'path_sot_in_potato')

# Lo que esta en POTATOrecon que no esta en SOT
path_potato_no_sot <-as.vector(path_potato[!path_potato%in%path_sot])
path_potato_no_sot
path_potato_no_sot <- unique(path_potato_no_sot)
write.csv(path_potato_no_sot, 'Path_potato_no_sot')

### EN ESTE PUNTO SE REALIZA VALIDACIÓN MANUAL CON LA BASE DE DATOS DE POTATOCYC, PARA VERIFICAR LAS ENZIMAS REPORTADAS PARA LAS RUTAS
###SE ELIMINAN LAS RUTAS QUE NO TIENEN ENZIMAS REPORTADAS EN POTATOCYC, O RUTAS QUE, AUNQUE TIENEN LA ENZIMA REPORTADA, LA REACCIÓN NO ES COHERENTE

### GAP FILL Solanum

# Ingresar Datos
reactions <- read.xls("~/PotatoRecon/Data/potatorecon_vr3.7_sinEX.xlsx",sheet = 1)
unique(as.vector(reactions$ID))

# Referencia para Solanum tuberosum desde KEGG
sot <- getReference("sot")

# Elimino las reacciones ya incluidas en la reconstrucción (sot)
sot_new <- sot[!sot$id%in%as.vector(reactions$ID),]
## se identifican todas las reacciones que estan en KEGG reportadas para Solanum y no estan en potatorecon
#Esto para incluir todas las reacciones de sot kegg en potatorecon
setwd("~/PotatoRecon/Results/")
sot_new
write.table(sot_new, 'kegg_sot', sep = ',')
## se identifican todas las reacciones que estan en KEGG reportadas para Solanum y no estan en potato-recon
## En este punto se incluyen manualmente las reacciones identificadas a potato-recon
##Antes de incluir las reacciones a potato recon, se utiliza el script "compartimentalize-KEGG.R 

# Identifico cuales reacciones no son KEGG reportadas para Solanum tuberosum
noKEGG<-reactions[!as.vector(reactions$ID)%in%sot$id,]
length(noKEGG)


###CURACIÓN DE REVERSIBILIDAD
##Primero se filtran todas las reacciones curadas manualmente (se filtran todas las reacciones irrevesribles)
## la curacion se realiza con todas las reversibles, aunque esten curadas manualmente, para tener puntos de validación
## En este punto se utiliza el script reversibily.R 

###GAP FILL CON TODO KEGG
##Despues de incluir todas las rxns de kegg reportadas para S. tuberosum, se realiza gap fill con todo kegg 

# Ingresar Datos
reactions <- read.xls("~/PotatoRecon/Data/potatorecon_vr3.8.xlsx",sheet = 1)# Se usan todas las rxn de potatorecon, incluyendo el gapfill de potato y la curación de reversibilidad

#Referencia todo KEGG
kegg <- as.data.frame(getReference('all'))
dim(kegg)
# Elimino las reacciones ya incluidas en la reconstrucción (kegg)
kegg_new <- kegg[!kegg$id%in%as.vector(reactions$ID),]
dim(kegg_new)

# GAP FILL KEGG all
new_react_KEGG_all <- gapFill(reactions$EQUATION, kegg_new$reaction, limit = 0,25, woCompartment = TRUE)
new_react_KEGG_all
#write.csv(new_react_KEGG_all, 'new_react_Kegg_all')

#Obtener los IDS y los EC de las reacciones fill
adata<-associated.data(new_react_KEGG_all,kegg_new,2)
ecs <- unique(unlist(strsplit(unique(as.vector(reactions$EC.NUMBER)),"[[:space:]]",perl = TRUE))) #obtenr los ecs de la reconstrucción
newtoinclude <- adata[adata$ec%in%ecs,]

###Antes de ejecutar los siguientes comando  obtener sot KEGG 
newnew <- adata[adata$ec%in%sot_summary$ec,]
newnew <- newnew[!is.na(newnew$ec),]
toinclude <- unique(rbind(newnew,newtoinclude))
write.table(toinclude, file = "gapfill_kegg_all", sep = ",")

#En este punto se realiza validación manual de las reacciones a incluir
#las reacciones del gap fil se compartimentalizan y se incluyen a la Potato_recon
#Este proceso permite identificar errores en la reconstrucción que son curados manualmente

### INCORPORAR GPRs A LA RECONSTRUCCIÓN 
map.gpr <- function(id, organism){
  id <- as.vector(id)
  gprs <- get.reference(organism)
  map <- function(id){
    return (as.vector(gprs[gprs[,"id"]%in%id,"gpr"]))
  }
  return(sapply(id, map))
}

###Aca ingreso potato2 que fue obtendio despues de cargar el modelo y verificar que tiene flujo para la función de biomasa
gpr<-unlist(map.gpr(Potato2$ID, "sot"))
class(gpr)
Potato2 <- as.data.frame.array(Potato2)
Potato2$GENE.ASSOCIATION[Potato2$ID%in%names(gpr)] <- as.vector(gpr)
setwd("~/Dropbox/model_analysis/Results/")
write.csv(Potato2, "Potato2.csv")

 
