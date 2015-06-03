rm(list = ls())
#définit le répertoire de travail
setwd("~/nemo_cluster/Sims2Conti")

#Charge les bibliothèques nécessaires
library("ggplot2", lib.loc="/home/pcubry/R/x86_64-pc-linux-gnu-library/3.1") ; library("Hmisc", lib.loc="/home/pcubry/R/x86_64-pc-linux-gnu-library/3.1") ; library("nlme", lib.loc="/home/pcubry/R/x86_64-pc-linux-gnu-library/3.1") ; library("lme4", lib.loc="/home/pcubry/R/x86_64-pc-linux-gnu-library/3.1") ; library(doBy)

#Définit une fonction basée sur les packages ggplot2 et Hmisc pour le plot de la moyenne et de l'écart-type
stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data=fun, geom=geom, width=0.2, ...)
}
stat_sum_single <- function(fun, geom="point", ...) {
  stat_summary(fun.y=fun, colour="red", geom=geom, size = 3, ...)
}
####Définit les variables utilisées dans les simulations et le nom des fichiers####

model <- c("2","3") #On ne considère plus le modèle 1 qui est un cas intermédiaire entre les autres et présente peu d'intérêt
cas <- c("3","4") #On ne considère que ces cas qui sont les plus vraisemblables
fecundity <- c("010","100","200") ; selfing <- c("0","03") ; qtls <- c("10","50") ; selection <- c("01","05","20","50")
repetition <- c("01","02","03","04","05","06","07","08","09","10")

#Création d'un objet vide qui va accueillir les résultats
results = NULL

########### Boucle de lecture des fichiers et calcul des variances ###########
for(m in model) {
  for(s in selfing) {
    for(q in qtls) {
      for(f in fecundity) {
        for(w in selection) {
          for(c in cas) {
            for(r in repetition){
            untar(tarfile = paste("Sim2_conti_Model",m,"_",s,"self_",q,"QTLs_fitperpop/quanti.tar.gz",sep=""),files = paste("quanti/Sim2_conti_Model",m,"_",s,"self_",q,"QTLs_Cas",c,"_SelInt",w,"_fec",f,"_200_",r,".quanti",sep=""),compressed = "gzip",exdir = paste("Sim2_conti_Model",m,"_",s,"self_",q,"QTLs_fitperpop",sep=""))
            temp <- read.table(paste("Sim2_conti_Model",m,"_",s,"self_",q,"QTLs_fitperpop/quanti/Sim2_conti_Model",m,"_",s,"self_",q,"QTLs_Cas",c,"_SelInt",w,"_fec",f,"_200_",r,".quanti",sep=""),header=TRUE)
            unlink(paste("Sim2_conti_Model",m,"_",s,"self_",q,"QTLs_fitperpop/quanti",sep=""),recursive = TRUE)
#Définition des noms de paramètres
            Model <- rep(paste(m,"_conti",sep=""),length(temp[,1])) ; Selfing <- rep(s,length(temp[,1])) ; QTLs <- rep(q,length(temp[,1])) ; Fecundity <- rep(f,length(temp[,1]))
            Selection <- rep(w,length(temp[,1])) ; Cas <- rep(c,length(temp[,1])) ; Repetition <- rep(r,length(temp[,1])) ; Environment <- rep(NA,length(temp[,1]))
            Population <- rep(NA,length(temp[,1]))
            temp_var<-NULL
#On ajoute les valeurs des paramètres à temp
          temp <- cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,Repetition,Population,Environment,temp) 
#On teste si toutes les pops sont occupées
          ifelse(length(unique(temp$pop))==24,
                 #Si oui, on met à jour les infos de Population et d'Environnement
                 {for (i in 1:3) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-1}
                 for (i in 4:6) {temp[temp$pop==i,]$Environment <-2 ; temp[temp$pop==i,]$Population <-1}
                 for (i in 7:9) {temp[temp$pop==i,]$Environment <-3 ; temp[temp$pop==i,]$Population <-2}
                 for (i in 10:12) {temp[temp$pop==i,]$Environment <-4 ; temp[temp$pop==i,]$Population <-2}
                 for (i in 13:15) {temp[temp$pop==i,]$Environment <-5 ; temp[temp$pop==i,]$Population <-3}
                 for (i in 16:18) {temp[temp$pop==i,]$Environment <-6 ; temp[temp$pop==i,]$Population <-3}
                 for (i in 19:21) {temp[temp$pop==i,]$Environment <-7 ; temp[temp$pop==i,]$Population <-4}
                 for (i in 22:24) {temp[temp$pop==i,]$Environment <-8 ; temp[temp$pop==i,]$Population <-4}
                 #Création de l'objet qui accueillera les résultats du calcul de variance
                 temp_var<-lme(G1~1,data = temp, random = ~1|Population/Environment/pop,na.action = na.fail)
                 #Création des composantes de la variance
                 #VarPopulation <- VarCorr(temp_var)$Population[1] #avec nlme intervals(temp_var)$reStruct$Population[,2]**2; names(VarPopulation) <- "VarPopulation"
                 VarPopulation<-VarCorr(temp_var)[2]
                 #VarEnvironment <- VarCorr(temp_var)$Environment[1] #avec nlme intervals(temp_var)$reStruct$Environment[,2]**2; names(VarEnvironment) <- "VarEnvironment"
                 VarEnvironment<-VarCorr(temp_var)[4]
                 #VarPatch <- VarCorr(temp_var)$pop[1] #avec nlme intervals(temp_var)$reStruct$pop[,2]**2; names(VarPatch) <- "VarPatch"
                 VarPatch<-VarCorr(temp_var)[6]
                 #VarInd <- sigma(temp_var)**2 #avec nlme temp_var$sigma**2 ; names(VarInd) <- "VarInd"
                 VarInd<-VarCorr(temp_var)[7]
                 },
                 #Sinon on mets des "NA" dans les colonnes
                 {VarPopulation <- NA ; VarEnvironment <- NA ; VarPatch <- NA ; VarInd <- NA})
    Model <- paste(m,"_conti",sep="") ; Selfing <- s ; QTLs <- q ; Fecundity <- f ; Selection <- w ; Cas <- c ; Repetition <- r
    res_temp<-as.data.frame(cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,Repetition,VarPopulation,VarEnvironment,VarPatch,VarInd))
    
    #Et on change la classe des variables variances de factor à numeric
    res_temp$VarPopulation<-as.numeric(levels(res_temp$VarPopulation))[res_temp$VarPopulation]
    res_temp$VarEnvironment<-as.numeric(levels(res_temp$VarEnvironment))[res_temp$VarEnvironment]
    res_temp$VarPatch<-as.numeric(levels(res_temp$VarPatch))[res_temp$VarPatch]
    res_temp$VarInd<-as.numeric(levels(res_temp$VarInd))[res_temp$VarInd]
    #On ajoute les colonnes qui vont recevoir le calcul des QSTs
    QSTpatch <- NA; QSTenv <-NA ; QSTpop <- NA
    res_temp<-cbind(res_temp,QSTpatch,QSTenv,QSTpop)
    #On fait le calcul des QSTs    
    res_temp$QSTpop<-(res_temp$VarPopulation/(res_temp$VarPopulation+res_temp$VarEnvironment+res_temp$VarPatch+2*res_temp$VarInd))
    res_temp$QSTenv<-(res_temp$VarEnvironment/(res_temp$VarEnvironment+res_temp$VarPatch+2*res_temp$VarInd))
    res_temp$QSTpatch<-(res_temp$VarPatch/(res_temp$VarPatch+2*res_temp$VarInd))
    results <- rbind(res_temp,results)
    #et on fait du ménage
    rm(list = c("QSTpop","QSTpatch","QSTenv","Model", "Selfing", "QTLs", "Fecundity", "Selection", "Cas", "Repetition", "Environment", "Population","temp_var","temp", "VarPopulation", "VarEnvironment", "VarPatch", "VarInd"))
            }
          }  
        }
      }
    }
  }
}

#####Et la même pour les modèles "nuls"#####
for(m in model) {
  for(s in selfing) {
    for(q in qtls) {
      for(f in fecundity) {
        for(w in selection) {
          for(c in cas) {
            for(r in repetition){
              untar(tarfile = paste("Sim2_conti_Model",m,"nul_",s,"self_",q,"QTLs_fitperpop/quanti.tar.gz",sep=""),files = paste("quanti/Sim2_conti_Model",m,"nul_",s,"self_",q,"QTLs_Cas",c,"_SelInt",w,"_fec",f,"_200_",r,".quanti",sep=""),compressed = "gzip",exdir = paste("Sim2_conti_Model",m,"nul_",s,"self_",q,"QTLs_fitperpop",sep=""))
              temp <- read.table(paste("Sim2_conti_Model",m,"nul_",s,"self_",q,"QTLs_fitperpop/quanti/Sim2_conti_Model",m,"nul_",s,"self_",q,"QTLs_Cas",c,"_SelInt",w,"_fec",f,"_200_",r,".quanti",sep=""),header=TRUE)
              unlink(paste("Sim2_conti_Model",m,"nul_",s,"self_",q,"QTLs_fitperpop/quanti",sep=""),recursive = TRUE)
              #Définition des noms de paramètres
              Model <- rep(paste(m,"_conti",sep=""),length(temp[,1])) ; Selfing <- rep(s,length(temp[,1])) ; QTLs <- rep(q,length(temp[,1])) ; Fecundity <- rep(f,length(temp[,1]))
              Selection <- rep(w,length(temp[,1])) ; Cas <- rep(c,length(temp[,1])) ; Repetition <- rep(r,length(temp[,1])) ; Environment <- rep(NA,length(temp[,1]))
              Population <- rep(NA,length(temp[,1]))
              temp_var<-NULL
              #On ajoute les valeurs des paramètres à temp
              temp <- cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,Repetition,Population,Environment,temp) 
              #On teste si toutes les pops sont occupées
              ifelse(length(unique(temp$pop))==24,
                     #Si oui, on met à jour les infos de Population et d'Environnement
              {for (i in 1:3) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-1}
               for (i in 4:6) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-1}
               for (i in 7:9) {temp[temp$pop==i,]$Environment <-2 ; temp[temp$pop==i,]$Population <-2}
               for (i in 10:12) {temp[temp$pop==i,]$Environment <-2 ; temp[temp$pop==i,]$Population <-2}
               for (i in 13:15) {temp[temp$pop==i,]$Environment <-3 ; temp[temp$pop==i,]$Population <-3}
               for (i in 16:18) {temp[temp$pop==i,]$Environment <-3 ; temp[temp$pop==i,]$Population <-3}
               for (i in 19:21) {temp[temp$pop==i,]$Environment <-4 ; temp[temp$pop==i,]$Population <-4}
               for (i in 22:24) {temp[temp$pop==i,]$Environment <-4 ; temp[temp$pop==i,]$Population <-4}
 #Création de l'objet qui accueillera les résultats du calcul de variance
 temp_var<-lme(G1~1,data = temp, random = ~1|Population/pop,na.action = na.fail)
 #Création des composantes de la variance
 
 #VarPopulation <- VarCorr(temp_var)$Population[1] #avec nlme intervals(temp_var)$reStruct$Population[,2]**2; names(VarPopulation) <- "VarPopulation"
 VarPopulation<-VarCorr(temp_var)[2]
 #VarEnvironment <- VarCorr(temp_var)$Environment[1] #avec nlme intervals(temp_var)$reStruct$Environment[,2]**2; names(VarEnvironment) <- "VarEnvironment"
 VarEnvironment<-NA
 #VarPatch <- VarCorr(temp_var)$pop[1] #avec nlme intervals(temp_var)$reStruct$pop[,2]**2; names(VarPatch) <- "VarPatch"
 VarPatch<-VarCorr(temp_var)[4]
 #VarInd <- sigma(temp_var)**2 #avec nlme temp_var$sigma**2 ; names(VarInd) <- "VarInd"
 VarInd<-VarCorr(temp_var)[5]
              },
#Sinon on mets des "NA" dans les colonnes
              {VarPopulation <- NA ; VarEnvironment <- NA; VarPatch <- NA ; VarInd <- NA})

              Model <- paste(m,"conti_nul") ; Selfing <- s ; QTLs <- q ; Fecundity <- f ; Selection <- w ; Cas <- c ; Repetition <- r
              res_temp<-as.data.frame(cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,Repetition,VarPopulation,VarEnvironment,VarPatch,VarInd))

#Et on change la classe des variables variances de factor à numeric
              res_temp$VarPopulation<-as.numeric(levels(res_temp$VarPopulation))[res_temp$VarPopulation]
              res_temp$VarEnvironment<-as.numeric(levels(res_temp$VarEnvironment))[res_temp$VarEnvironment]
              res_temp$VarPatch<-as.numeric(levels(res_temp$VarPatch))[res_temp$VarPatch]
              res_temp$VarInd<-as.numeric(levels(res_temp$VarInd))[res_temp$VarInd]
#On ajoute les colonnes qui vont recevoir le calcul des QSTs
              QSTpatch <- NA; QSTenv <-NA ; QSTpop <- NA
              res_temp<-cbind(res_temp,QSTpatch,QSTenv,QSTpop)
#On fait le calcul en itérant sur la liste
              for (i in 1:length(res_temp$Model)) {res_temp$QSTpop<-(res_temp$VarPopulation/(res_temp$VarPopulation+res_temp$VarPatch+2*res_temp$VarInd))}
              for (i in 1:length(res_temp$Model)) {res_temp$QSTpatch<-(res_temp$VarPatch/(res_temp$VarPatch+2*res_temp$VarInd))}

#et on ajoute temp aux résultats !!!
results <- rbind(res_temp,results)
#et on fait du ménage
rm(list = c("QSTpatch","QSTpop","QSTenv","res_temp","Model", "Selfing", "QTLs", "Fecundity", "Selection", "Cas", "Repetition", "Environment", "Population","temp_var","temp", "VarPopulation", "VarEnvironment", "VarPatch", "VarInd"))
            }
          }  
        }
      }
    }
  }
}

#####On lit aussi les résultats des initialisation pour comparaison#####
setwd("~/nemo_cluster")
for(s in selfing) {
  for(q in qtls) {
              untar(tarfile = "INIT2.tar.gz",files = paste("INIT2/conti_init2_carto_",s,"self_",q,"QTLs/quanti/conti_init2_carto_",s,"self_",q,"QTLs.bin_200000_1.quanti",sep=""),compressed = "gzip",exdir = "INIT2")
              temp <- read.table(paste("INIT2/INIT2/conti_init2_carto_",s,"self_",q,"QTLs/quanti/conti_init2_carto_",s,"self_",q,"QTLs.bin_200000_1.quanti",sep=""),header=TRUE)
              #On ne garde que les adultes
              temp<-temp[which(temp$age==2),]
              #Et on efface les fichiers extraits
              unlink("INIT2",recursive = TRUE)
              #Définition des noms de paramètres
              Model <- rep("equilibrium",length(temp[,1])) ; Selfing <- rep(s,length(temp[,1])) ; QTLs <- rep(q,length(temp[,1])) ; Fecundity <- rep("up to capacity",length(temp[,1]))
              Selection <- rep("minimal",length(temp[,1])) ; Cas <- rep("NA",length(temp[,1])) ; Repetition <- rep("NA",length(temp[,1])) ; Environment <- rep(NA,length(temp[,1]))
              Population <- rep(NA,length(temp[,1]))
              temp_var<-NULL
              #On ajoute les valeurs des paramètres à temp
              temp <- cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,Repetition,Population,Environment,temp) 
              #On teste si toutes les pops sont occupées
              ifelse(length(unique(temp$pop))==24,
                     #Si oui, on met à jour les infos de Population et d'Environnement
{for (i in 1:3) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-1}
 for (i in 4:6) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-1}
 for (i in 7:9) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-2}
 for (i in 10:12) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-2}
 for (i in 13:15) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-3}
 for (i in 16:18) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-3}
 for (i in 19:21) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-4}
 for (i in 22:24) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-4}
 #Création de l'objet qui accueillera les résultats du calcul de variance
 temp_var<-lme(G1~1,data = temp, random = ~1|Population/pop,na.action = na.fail)
 #Création des composantes de la variance
 
 #VarPopulation <- VarCorr(temp_var)$Population[1] #avec nlme intervals(temp_var)$reStruct$Population[,2]**2; names(VarPopulation) <- "VarPopulation"
 VarPopulation<-VarCorr(temp_var)[2]
 #VarEnvironment <- VarCorr(temp_var)$Environment[1] #avec nlme intervals(temp_var)$reStruct$Environment[,2]**2; names(VarEnvironment) <- "VarEnvironment"
 VarEnvironment<-NA
 #VarPatch <- VarCorr(temp_var)$pop[1] #avec nlme intervals(temp_var)$reStruct$pop[,2]**2; names(VarPatch) <- "VarPatch"
 VarPatch<-VarCorr(temp_var)[4]
 #VarInd <- sigma(temp_var)**2 #avec nlme temp_var$sigma**2 ; names(VarInd) <- "VarInd"
 VarInd<-VarCorr(temp_var)[5]#Création de l'objet qui accueillera les résultats du calcul de variance
},
#Sinon on mets des "NA" dans les colonnes
{VarPopulation <- NA ; VarEnvironment <- NA; VarPatch <- NA ; VarInd <- NA})

Model <- paste("conti_equilibrium") ; Selfing <- s ; QTLs <- q ; Fecundity <- "up to capacity" ; Selection <- "minimal" ; Cas <- NA ; Repetition <- NA
res_temp<-as.data.frame(cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,Repetition,VarPopulation,VarEnvironment,VarPatch,VarInd))

#Et on change la classe des variables variances de factor à numeric
res_temp$VarPopulation<-as.numeric(levels(res_temp$VarPopulation))[res_temp$VarPopulation]
res_temp$VarEnvironment<-as.numeric(levels(res_temp$VarEnvironment))[res_temp$VarEnvironment]
res_temp$VarPatch<-as.numeric(levels(res_temp$VarPatch))[res_temp$VarPatch]
res_temp$VarInd<-as.numeric(levels(res_temp$VarInd))[res_temp$VarInd]
#On ajoute les colonnes qui vont recevoir le calcul des QSTs
QSTpatch <- NA; QSTenv <-NA ; QSTpop <- NA
res_temp<-cbind(res_temp,QSTpatch,QSTenv,QSTpop)
#On fait le calcul en itérant sur la liste
for (i in 1:length(res_temp$Model)) {res_temp$QSTpop<-(res_temp$VarPopulation/(res_temp$VarPopulation+res_temp$VarPatch+2*res_temp$VarInd))}
for (i in 1:length(res_temp$Model)) {res_temp$QSTpatch<-(res_temp$VarPatch/(res_temp$VarPatch+2*res_temp$VarInd))}

#et on ajoute temp aux résultats !!!
results <- rbind(res_temp,results)
#et on fait du ménage
rm(list = c("QSTpatch","QSTpop","QSTenv","res_temp","Model", "Selfing", "QTLs", "Fecundity", "Selection", "Cas", "Repetition", "Environment", "Population","temp_var","temp", "VarPopulation", "VarEnvironment", "VarPatch", "VarInd"))
  }
}
setwd("~/nemo_cluster/Sims2Conti")

####On stocke les résultats dans un fichier
write.table(results, file="QSTs_Sims2Conti.txt")
####

####On représente les QSTs####
q<-ggplot(data = results)
q<-q+(aes(x=Selection,y=QSTpop))+ylim(0,1)+geom_point() +aes(colour = QTLs:Model, shape = Model,size=3)+guides(size=FALSE,alpha=FALSE,colour=guide_legend("Interactions of Model and QTLs"))+facet_grid(Selfing~Cas+Fecundity,labeller=label_both)
q+labs(title="QST between Populations")
q<-ggplot(data = results)
q<-q+(aes(x=Selection))+ylim(0,1)+geom_point(aes(y=QSTenv))+ aes(colour = QTLs:Model, shape = QTLs,size=4,alpha=0.5)+guides(size=FALSE,alpha=FALSE,colour=guide_legend("Interactions of Model and QTLs"))+facet_grid(Selfing~Cas+Fecundity,labeller=label_both)
q+labs(title="QST between Environments")
q<-ggplot(data = results)
q<-q+(aes(x=Selection))+ylim(0,1)+geom_point(aes(y=QSTpatch))+ aes(colour = QTLs:Model, shape = QTLs,size=4,alpha=0.5)+guides(size=FALSE,alpha=FALSE,colour=guide_legend("Interactions of Model and QTLs"))+facet_grid(Selfing~Cas+Fecundity,labeller=label_both)
q+labs(title="QST between Patches")
