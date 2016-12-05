library(devtools)
install_github("dosorio/convert2sybil")
library(convert2sybil)
library(gdata)
library(glpkAPI)
library(sybil)
library(minval)

###CARGAR EL MODELO INICIAL
setwd("~/Dropbox/model_analysis/Results/")
Potato1 <- read.xls("~/Dropbox/model_analysis/Data/model_sybil/Potato-F1.2b.xlsx", sheet = 1)
Potato1 <- as.data.frame.array(Potato1)
#Potato1 <- Potato1[!Potato1$X.=='#',]
convert2sybil(reactions = Potato1$EQUATION,abbreviation = as.vector(Potato1$ID),file = "Potato1")
#convert2sybil(abbreviation = as.vector(Potato1$ID),reactions = Potato1$EQUATION,gene = as.vector(Potato1$GENE.ASSOCIATION),objective_function = 2054,file = "Potato1")
model <- readTSVmod(prefix = "Potato1",fpath = "~/Dropbox/model_analysis/Results",quoteChar = "\"")
warnings()
findExchReact(model)

#NUEVA VERSIÓN PARA GENERAR EL MODELO
Potato2 <- read.csv("~/Dropbox/model_analysis/Results/Potato2.csv")
Potato2 <- Potato2[,c(4,5,6,7,8,9,10,11,12)]
names(Potato2) <- c("ID","NAME","REACTION","EC.NUMBER","GPR","LOWER.BOUND","UPPER.BOUND","OBJECTIVE","COMPARTMENT")
Potato2 <- as.data.frame.array(Potato2)
to.sbml(Potato2,"Potato2.xml")
model1 <- readSBMLmod("Potato2.xml")
model1 <- addExchReact(model1,"Biomass[c]",ub = 1000)
findExchReact(model1)
optimizeProb(model1)


####EXPERIMENTO GENERAL1
### VALIDAR CAPACIDADES METABOLICAS

model <- addReact(model = model, id = "FM1",met = c("H2O[c]","CO2[c]","Photon[c]","Oxygen[c]"),Scoef = c(-1,-1,-1,1),obj = 0)
plot(getFluxDist(optimizeProb(model)),cex=0.1)
model <- addReact(model = model, id = "FM1",met = c("H2O[c]","CO2[c]","Photon[c]","Glycerone-phosphate[c]"),Scoef = c(-1,-1,-1,1),obj = 1)
plot(getFluxDist(optimizeProb(model)),cex=0.1)
model <- addReact(model = model, id = "FM1",met = c("H2O[c]","CO2[c]","Photon[c]","Starch[c]"),Scoef = c(-1,-1,-1,1),obj = 1)
plot(getFluxDist(optimizeProb(model)),cex=0.1)
model <- addReact(model = model, id = "FM1",met = c("H2O[c]","CO2[c]","Photon[c]","Sucrose[c]"),Scoef = c(-1,-1,-1,1),obj = 1)
plot(getFluxDist(optimizeProb(model)),cex=0.1)
model <- addReact(model = model, id = "FM1",met = c("H2O[c]","CO2[c]","Photon[c]","Cellulose[c]"),Scoef = c(-1,-1,-1,1),obj = 1)
plot(getFluxDist(optimizeProb(model)),cex=0.1)
model <- addReact(model = model, id = "FM1",met = c("H2O[c]","CO2[c]","Photon[c]","D-Xylose[c]"),Scoef = c(-1,-1,-1,1),obj = 1)
plot(getFluxDist(optimizeProb(model)),cex=0.1)
model <- addReact(model = model, id = "FM1",met = c("Glycerone-phosphate[c]","Starch[c]"),Scoef = c(-1,1),obj = 1)
plot(getFluxDist(optimizeProb(model)),cex=0.1)
model <- addReact(model = model, id = "FM1",met = c("H2O[c]","CO2[c]","Photon[c]", "ADP[c]", "H*[c]", "Orthophosphate[c]"),Scoef = c(-1,-1,-1,1,1,1),obj = 1)
plot(getFluxDist(optimizeProb(model)),cex=0.1)


###COMANDOS GENERALES

# Set Objective Function
model@obj_coef[model@react_id=="RBS01"]<-1 #activar
optimizeProb(model)

# Testing Biomass compounds
New <- NULL
model <- addReact(model,"T1",c("D-Xylose[c]"),Scoef = c(1),obj = 1) # adiciona reacciones y la pone como funcipon objetivo
optimizeProb(model)
New <- unique(c(New,model@react_id[getFluxDist(optimizeProb(model))!=0]))
getFluxDist(optimizeProb(model))[model@react_id=="T1"]
plot(getFluxDist(optimizeProb(model)))
model<-rmReact(model = model,react = "T") #elimina reacciones

####EXPERIMENTO GENERAL2
###CAMBIAR LAS REVERSIBILIDADES DE POTATO DE ACUERDO A LOS METABOLITOS QUE NO TIENE FLUJO EN LA FUNCIÓN OBJETIVO

##Aunque este paso se hizo, no fueron los experimentos exitosos para continuar con el proceso
##Sacar las reaaciones que tienen flujo en la recon normal de los metabolitos que no tienen flujo en la función de biomasa
#1 EXPERIMENTO POR METABOLITO
#Verificar flujo por metabolito
model <- addReact(model,"T1",c("4-Coumaryl-alcohol[c]"),Scoef = c(1),obj = 1) # adiciona reacciones y la pone como funcipon objetivo
optimizeProb(model)
# Aca se realizan los mismos pasos del experimento 2 despues de hacer el FBA

# 2 EXPERIMENTO POR FUNCIÓN METABOLICA (del metabolito sin flujo)
#verificar flujo por función metabolica
model_bloq <- addReact(model = model, id = "FM1",met = c("H2O[c]","CO2[c]","Photon[c]","4-Coumaryl-alcohol[c]"),Scoef = c(-1,-1,-1,1),obj = 1)
optimizeProb(model_bloq)## FBA
#obtener las reacciones que tienen flujo para la función objetibo
react_flux<-model_bloq@react_id[getFluxDist(optimizeProb(model_bloq))!=0]
comoestan <- model@react_id[getFluxDist(optimizeProb(model))>0]

##Cambiar reversibilidad de todas las reacciones
model_rev <-model
lowbnd(model_rev)
lowbnd(model_rev)<- rep(-1000,model_rev@react_num)
uppbnd(model_rev)<- rep(1000,model_rev@react_num)

##Sacar las reaaciones que tienen flujo en la recon_rev y la normal de los metabolitos que no tienen flujo en la función de biomasa
#para experimeto1
model_rev <- addReact(model_rev,"RSB01",c("4-Coumaryl-alcohol[c]"),Scoef = c(1),obj = 1)
#para experimento 2
#model_rev <- addReact(model = model_rev, id = "FM1",met = c("H2O[c]","CO2[c]","Photon[c]","4-Coumaryl-alcohol[c]"),Scoef = c(-1,-1,-1,1),obj = 1)
optimizeProb(model_rev)
#model_Rrev<-model_rev@react_id[getFluxDist(optimizeProb(model_rev))!=0]
comodeberia <- model_rev@react_id[getFluxDist(optimizeProb(model_rev))<0]

###Identificar reacciones para cambiar reversibilidad
id_rev<-as.vector(comodeberia[!comodeberia%in%react_flux])
#ver_rev<-comoestan[comoestan%in%comodeberia]
#write.table(comoestan[comoestan%in%comodeberia],file = "EstablecerReversibles.txt")

##Cambiar reversibilidad recon
Potato1$LOWER.BOUND[Potato1$ID%in%id_rev]#Verificación 
Potato1$LOWER.BOUND[Potato1$ID%in%id_rev]=-1000 #Cambio
Potato1$LOWER.BOUND[Potato1$ID%in%id_rev]#Verificación
setwd("~/Dropbox/model_analysis/Results/")
write.table(as.vector(Potato1$EQUATION[Potato1$ID%in%id_rev]), file = "equaciones_rev_Coumaryl-alcohol")

as.vector(Potato1$EQUATION[Potato1$ID%in%id_rev])

Potato1$EQUATION[Potato1$ID%in%id_rev] <- gsub(" => ", " <=> ", as.vector(Potato1$EQUATION[Potato1$ID%in%id_rev]))
Potato1$EQUATION[Potato1$ID%in%id_rev]

####EXPERIMENTO GENERAL3
##GENERAR MATRIZ ESTEQUIOMETRICA
#Con esta matriz 
matrix<-stoichiometric.matrix(Potato1$EQUATION)
bal_met<- as.vector(rowSums(matrix))
bal_met[bal_met!=0]

####EXPERIMENTO GENERAL4: Exp usado para generar el modelo final
### CAMBIAR REVERSIBILIDADES DE POTATO PARA OBTENER FLUJO EN LA REACCIÓN DE BIOMASA
# EXPERIMENTO POR FUNCIÓN OBJETIVO
#Cambiar las reversibilidades del modelo y dejarlo completamente reversible
lowbnd(model)<-rep(-1000,model@react_num)
uppbnd(model)<-rep(1000,model@react_num)
#Asignar nueva variable para potato
Potato2 <- Potato1
#Cambiar las reversibilidades en potato2 para las reacciones que deben ser reversibles cuando se corre el modelo completamente reversible
#lr <- model@react_id[getFluxDist(optimizeProb(model))>0]# No es necesario, Verificar(irrevesribles forward)
rl <- model@react_id[getFluxDist(optimizeProb(model))<0] #(reversibles)
#Potato2$EQUATION[Potato2$ID%in%lr] <- gsub("<?=>","=>",Potato2$EQUATION[Potato2$ID%in%lr])#No es necesario
Potato2$EQUATION[Potato2$ID%in%rl] <- gsub("<?=>","<=>",Potato2$EQUATION[Potato2$ID%in%rl])
Potato2$LOWER.BOUND[Potato2$ID%in%rl]<--1000
write.csv(Potato2, 'Potato2p.csv')

#Cargar el modelo para verificación
Potato2 <- read.csv("~/Dropbox/model_analysis/Results/Potato2p.csv")
Potato2 <- as.data.frame.array(Potato2)
#convert2sybil(reactions = Potato2$EQUATION,abbreviation = as.vector(Potato2$ID),file = "Potato2")
convert2sybil(abbreviation = as.vector(Potato2$ID),reactions = Potato2$EQUATION,gene = as.vector(Potato2$GENE.ASSOCIATION),objective_function = 2054,file = "Potato2")
model2 <- readTSVmod(prefix = "Potato2",fpath = "~/Dropbox/model_analysis/Results",quoteChar = "\"")
#cuando se cargo el modelo, la función objetivo aun no tenia flujo. Despues de hacer varios experimentos se identifico que se debia forzar el modelo a que tomara la siguiente función de biomasa
model2 <- addExchReact(model = model2,met = c("Biomass[c]"),lb = 0,ub = 1000)
optimizeProb(model2)
#En este punto la reacción de reaccion de biomasa ya tiene flujo
#En este punto se deben incluir nuevamente las gprs a Potato2.csv (script sot_g2f) ya que se desorganizan.


####EXPERIMENTO GENERAL5
###EXPERIMENTOS PARA VERIFICAR QUE DIF. FUN-OBJ POR SEPARADO TENGAN FLUJO CON EL MODELO COMPLETAMENTE REVERSIBLE
##1 Reacción de biomasa de biomasa(en convert2sybil se ha definido la función objetivo)
## Este experiemento tambien se hizo incialmente para darle flujo a la reacción de biomasa
#Verificación de funcion objetivo para el modelo totalmente reversible
#Tener en cuenta que debo cargar nuevamente el modelo inicial (potato-f1.2)
lowbnd(model)<-rep(-1000,model@react_num)
uppbnd(model)<-rep(1000,model@react_num)
#model1 <- rmReact(model = model1,react = "RBS01") Eliminar reacciones del modelo
optimizeProb(model)
model <- addExchReact(model = model,met = c("Biomass[c]"),ub = 1000)
model <- addReact(model,id = "B",met = c("ATP[c]","dTTP[c]","dCTP[c]","UTP[c]","CTP[c]","dGTP[c]","dATP[c]","GTP[c]","Hexadecanoic-acid[c]","L-Lysine[c]","L-Valine[c]","L-Tyrosine[c]","L-Tryptophan[c]","L-Threonine[c]","L-Serine[c]","L-Proline[c]","L-Phenylalanine[c]","L-Methionine[c]","L-Leucine[c]","L-Isoleucine[c]","Glycine[c]","L-Glutamine[c]","L-Glutamate[c]","L-Cysteine[c]","L-Aspartate[c]","L-Asparagine[c]","L-Alanine[c]","Sinapyl-alcohol[c]","Coniferyl-alcohol[c]","4-Coumaryl-alcohol[c]","D-Xylose[c]","Sucrose[c]","alpha-D-Glucose[c]","beta-D-Fructose[c]","Starch[c]","Cellulose[c]","Biomass[c]"),
                   Scoef = c(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1),obj = 1)
optimizeProb(model)
table(getFluxDist(optimizeProb(model))!=0)
plot(getFluxDist(optimizeProb(model)),cex=0.2)
findExchReact(model)

##2 Reacción de fijación de carbono 
##Verificación de funcion objetivo para el modelo totalmente reversible
model1@obj_coef[model1@react_id=="RBS01"]<-0 
model1@obj_coef[model1@react_id=="R00024"]<-1 #activar
optimizeProb(model1)
##Verificación de la función objetivo con la versión de potato (potato2.csv) a la cual se le cambio la reversibilidad de algunas reacciones para tener flujo en RBS01
#Esto se hace depues de que se genera potato2(el que tiene cambio de reversibilidades para tner flujo en la reacción de biomasa)
#Cambiar reversibilidades de reacciones particulares cuando sea necesario para validar
#Potato2$EQUATION[71]<- gsub("<=>", "=>", Potato2$EQUATION[71]## 71 es la posición de la reacción para fijación de carbono
convert2sybil(reactions = Potato2$EQUATION,abbreviation = as.vector(Potato2$ID),file = "Potato2")
model2 <- readTSVmod(prefix = "Potato2",fpath = "~/Dropbox/model_analysis/Results",quoteChar = "\"")
model2 <- addExchReact(model = model1,met = c("Biomass[c]"),ub = 1000)
optimizeProb(model2)

###INCLUIR REACCIÓNES DE INTERCAMBIO PARA LOS DEAD-ENDS
## Este experimento se hizo para intentar que diferentes funciones objetivo tengas las reacciones de Uptake/Excretion
#Antes de hacer este paso se verifican las reacciones de Uptake/Excretion con el script (Exper_simulation.R)
dem <- deadEndMetabolites(model2)
for(metabolite in dem$dem){
  model2 <- addExchReact(model = model2,met = c(metabolite),lb = 0,ub = 1000)
}
optimizeProb(model2) 

#Posteriormente se verfican nuevamente las reacciones de Uptake/Excretion con el script (Exper_simulation.R)
#No se evidencian los cambios esperados

###CAMBIAR REVERSIBILIDADES DE POTATO PARA OBTENER FLUJO EN EL MODELO CUADO SE LE EXIGE QUE CUMPLA CON VARIAS FUNCIONES OBJETIVO 
#las funciones objetivo de definen para evaluar las condiciones fisiologicas definidas en el proyecto
#La descripción de cada función objetivo se esta en el xls (review_recon-plantas)

# Cargar modelo 
setwd("~/Dropbox/model_analysis/Results/")
Potato1 <- read.xls("~/Dropbox/model_analysis/Data/model_sybil/Potato-F1.2.xlsx", sheet = 1)
Potato1 <- as.data.frame.array(Potato1)
convert2sybil(reactions = Potato1$EQUATION,abbreviation = as.vector(Potato1$ID),file = "Potato1")
model <- readTSVmod(prefix = "Potato1",fpath = "~/Dropbox/model_analysis/Results",quoteChar = "\"")
#Verificar que toma solo los metabolitos deseados uptake
warnings()
findExchReact(model)
#Establecer varias funciones objetivo
#model@obj_coef[model@react_id== c("RBS01","R00024","RK0004","RK0005","R01015","R02110","R05196","R00028")]<-1
model@obj_coef[model@react_id=="R00028"]<-1# serealiza para todas las funciones objetivo

#Cambiar reversibilidades
lowbnd(model)<-rep(-1000,model@react_num)
uppbnd(model)<-rep(1000,model@react_num)
model <- addExchReact(model = model,met = c("Biomass[c]"),lb = 0,ub = 1000)
model
Potato3 <- Potato1
rl <- model@react_id[getFluxDist(optimizeProb(model))<0] #(reversibles)
Potato3$EQUATION[Potato3$ID%in%rl] <- gsub("<?=>","<=>",Potato3$EQUATION[Potato3$ID%in%rl])
Potato3$LOWER.BOUND[Potato3$ID%in%rl]<--1000
convert2sybil(reactions = Potato3$EQUATION,abbreviation = as.vector(Potato3$ID),file = "Potato3")
model3 <- readTSVmod(prefix = "Potato3",fpath = "~/Dropbox/model_analysis/Results",quoteChar = "\"")
model3 <- addExchReact(model = model3,met = c("Biomass[c]"),ub = 1000)

#verificar entradas y salidas
#Da lo mismo





