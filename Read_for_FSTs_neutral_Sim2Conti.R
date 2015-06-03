rm(list=ls())
#définit le répertoire de travail
setwd("~/nemo_cluster/Sims2Conti")

#Charge les bibliothèques nécessaires
library("ggplot2", lib.loc="/home/pcubry/R/x86_64-pc-linux-gnu-library/3.1") ; library("Hmisc", lib.loc="/home/pcubry/R/x86_64-pc-linux-gnu-library/3.1") ; library("nlme", lib.loc="/home/pcubry/R/x86_64-pc-linux-gnu-library/3.1") ; library("lme4", lib.loc="/home/pcubry/R/x86_64-pc-linux-gnu-library/3.1") ; library(doBy)
library(adegenet) ; library(hierfstat)
####Définit une fonction basée sur les packages ggplot2 et Hmisc pour le plot de la moyenne et de l'écart-type####
stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data=fun, geom=geom, width=0.2, ...)
}
stat_sum_single <- function(fun, geom="point", ...) {
  stat_summary(fun.y=fun, colour="red", geom=geom, size = 3, ...)
}

####Définit les variables utilisées dans les simulations et le nom des fichiers####

model <- c("2","3") #On ne considère plus le modèle 1 qui est un cas intermédiaire entre les autres et présente peu d'intérêt
cas <- c("4") #On ne considère que ces cas qui sont les plus vraisemblables
fecundity <- c("010","100","200") ; selfing <- c("0","03") ; qtls <- c(10,50) ; selection <- c("01","05","20","50")
repetition <- c("01","02","03","04","05","06","07","08","09","10")

########### Boucle de lecture des fichiers ###########
for(m in model)
{
  for(s in selfing)
  {
    for(q in qtls)
    {
      untar(tarfile = paste("Sim2_conti_Model",m,"_",s,"self_",q,"QTLs_fitperpop/neutral.tar.gz",sep=""),compressed = "gzip",exdir = paste("Sim2_conti_Model",m,"_",s,"self_",q,"QTLs_fitperpop",sep=""))
      for(f in fecundity)
      {
        for(w in selection)
        {
          for(c in cas)
          {
            for(r in repetition)
            {
              if(file.exists(paste("VarCompForFSTs_neutralConti/var_comp_model",m,"_conti_cas",c,"_",f,"fec_",s,"self_selection",w,"_",q,"QTLs_rep",r,"_neutral.rda",sep=""))==FALSE)  
              {
                temp <- read.fstat(paste("Sim2_conti_Model",m,"_",s,"self_",q,"QTLs_fitperpop/neutral/Sim2_conti_Model",m,"_",s,"self_",q,"QTLs_Cas",c,"_SelInt",w,"_fec",f,"_200_",r,".dat",sep=""))
                #Définition des noms de paramètres
                Model <- rep(m,length(temp[,1])) ; Selfing <- rep(s,length(temp[,1])) ; QTLs <- rep(q,length(temp[,1])) ; Fecundity <- rep(f,length(temp[,1]))
                Selection <- rep(w,length(temp[,1])) ; Cas <- rep(c,length(temp[,1])) ; Repetition <- rep(r,length(temp[,1])) ; Environment <- rep(NA,length(temp[,1]))
                Population <- rep(NA,length(temp[,1]))
                
                #On vire les colonnes qui ne nous intéressent pas
                temp <- temp[c(seq(1,10001,20))]
                #On ajoute les valeurs des paramètres à temp
                temp <- cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,Repetition,Population,Environment,temp)
                #On teste si toutes les pops sont occupées
                if(length(unique(temp$Pop))==24)
                {
                  #Si oui, on met à jour les infos de Population et d'Environnement
                  for (i in 1:3) {temp[temp$Pop==i,]$Environment <-1 ; temp[temp$Pop==i,]$Population <-1}
                  for (i in 4:6) {temp[temp$Pop==i,]$Environment <-2 ; temp[temp$Pop==i,]$Population <-1}
                  for (i in 7:9) {temp[temp$Pop==i,]$Environment <-3 ; temp[temp$Pop==i,]$Population <-2}
                  for (i in 10:12) {temp[temp$Pop==i,]$Environment <-4 ; temp[temp$Pop==i,]$Population <-2}
                  for (i in 13:15) {temp[temp$Pop==i,]$Environment <-5 ; temp[temp$Pop==i,]$Population <-3}
                  for (i in 16:18) {temp[temp$Pop==i,]$Environment <-6 ; temp[temp$Pop==i,]$Population <-3}
                  for (i in 19:21) {temp[temp$Pop==i,]$Environment <-7 ; temp[temp$Pop==i,]$Population <-4}
                  for (i in 22:24) {temp[temp$Pop==i,]$Environment <-8 ; temp[temp$Pop==i,]$Population <-4}
                }
                
                if (unique(is.na(temp$Population)==FALSE))
                {
                  var_comp<-varcomp.glob(data.frame(temp$Population,temp$Environment,temp$Pop),(temp[,-c(1:10)])) ;
                  save(file=paste("VarCompForFSTs_neutralConti/var_comp_model",m,"_conti_cas",c,"_",f,"fec_",s,"self_selection",w,"_",q,"QTLs_rep",r,"_neutral.rda",sep=""),var_comp) 
                }
                Model <- m ; Selfing <- s ; QTLs <- q ; Fecundity <- f ; Selection <- w ; Cas <- c ; Repetition <- r
                #et on fait du ménage
                rm(list = c("var_comp","Model", "Selfing", "QTLs", "Fecundity", "Selection", "Cas", "Repetition", "Environment", "Population","temp"))
              }
            }
          }  
        }
      }
      unlink(paste("Sim2_conti_Model",m,"_",s,"self_",q,"QTLs_fitperpop/neutral",sep=""),recursive = TRUE)
    }
  }
}

########### Boucle de lecture des fichiers pour modèles "nuls" ###########
for(m in model)
{
  for(s in selfing)
  {
    for(q in qtls)
    {
      untar(tarfile = paste("Sim2_conti_Model",m,"nul_",s,"self_",q,"QTLs_fitperpop/neutral.tar.gz",sep=""),compressed = "gzip",exdir = paste("Sim2_conti_Model",m,"nul_",s,"self_",q,"QTLs_fitperpop",sep=""))
      for(f in fecundity)
      {
        for(w in selection)
        {
          for(c in cas)
          {
            for(r in repetition)
            {
  if(file.exists(paste("VarCompForFSTsNul_neutralConti/var_comp_model",m,"_conti_nul_cas",c,"_",f,"fec_",s,"self_selection",w,"_",q,"QTLs_rep",r,"_neutral.rda",sep=""))==FALSE)  
{
              temp <- read.fstat(paste("Sim2_conti_Model",m,"nul_",s,"self_",q,"QTLs_fitperpop/neutral/Sim2_conti_Model",m,"nul_",s,"self_",q,"QTLs_Cas",c,"_SelInt",w,"_fec",f,"_200_",r,".dat",sep=""))
              #Définition des noms de paramètres
              Model <- rep(paste(m,"nul",sep=""),length(temp[,1])) ; Selfing <- rep(s,length(temp[,1])) ; QTLs <- rep(q,length(temp[,1])) ; Fecundity <- rep(f,length(temp[,1]))
              Selection <- rep(w,length(temp[,1])) ; Cas <- rep(c,length(temp[,1])) ; Repetition <- rep(r,length(temp[,1])) ; Environment <- rep(NA,length(temp[,1]))
              Population <- rep(NA,length(temp[,1]))
              #On vire les colonnes qui ne nous intéressent pas
              temp <- temp[c(seq(1,10001,20))]
              #On ajoute les valeurs des paramètres à temp
              temp <- cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,Repetition,Population,Environment,temp)
              #On teste si toutes les pops sont occupées
              if(length(unique(temp$Pop))==24)
              {
                #Si oui, on met à jour les infos de Population et d'Environnement
                for (i in 1:3) {temp[temp$Pop==i,]$Environment <-1 ; temp[temp$Pop==i,]$Population <-1}
                for (i in 4:6) {temp[temp$Pop==i,]$Environment <-1 ; temp[temp$Pop==i,]$Population <-1}
                for (i in 7:9) {temp[temp$Pop==i,]$Environment <-2 ; temp[temp$Pop==i,]$Population <-2}
                for (i in 10:12) {temp[temp$Pop==i,]$Environment <-2 ; temp[temp$Pop==i,]$Population <-2}
                for (i in 13:15) {temp[temp$Pop==i,]$Environment <-3 ; temp[temp$Pop==i,]$Population <-3}
                for (i in 16:18) {temp[temp$Pop==i,]$Environment <-3 ; temp[temp$Pop==i,]$Population <-3}
                for (i in 19:21) {temp[temp$Pop==i,]$Environment <-4 ; temp[temp$Pop==i,]$Population <-4}
                for (i in 22:24) {temp[temp$Pop==i,]$Environment <-4 ; temp[temp$Pop==i,]$Population <-4}
              }
              
           if (unique(is.na(temp$Population)==FALSE))
              {
                var_comp<-varcomp.glob(data.frame(temp$Population,temp$Pop),(temp[,-c(1:10)])) ;
                save(file=paste("VarCompForFSTsNul_neutralConti/var_comp_model",m,"_conti_nul_cas",c,"_",f,"fec_",s,"self_selection",w,"_",q,"QTLs_rep",r,"_neutral.rda",sep=""),var_comp) 
              }
                 Model <- m ; Selfing <- s ; QTLs <- q ; Fecundity <- f ; Selection <- w ; Cas <- c ; Repetition <- r
                #et on fait du ménage
                 rm(list = c("var_comp","Model", "Selfing", "QTLs", "Fecundity", "Selection", "Cas", "Repetition", "Environment", "Population","temp"))
              }
            }
          }  
        }
      }
      unlink(paste("Sim2_conti_Model",m,"nul_",s,"self_",q,"QTLs_fitperpop/neutral",sep=""),recursive = TRUE)
    }
  }
}

########### Boucle de lecture des fichiers pour modèles "initiaux" ###########
setwd("~/nemo_cluster")
for(s in selfing)
{
  for(q in qtls)
  {
    untar(tarfile = "INIT2.tar.gz",files = paste("INIT2/conti_init2_carto_",s,"self_",q,"QTLs/neutral/conti_init2_carto_",s,"self_",q,"QTLs.bin_200000_1.dat",sep=""),compressed = "gzip",exdir = "INIT2")
    temp <- read.fstat(paste("INIT2/INIT2/conti_init2_carto_",s,"self_",q,"QTLs/neutral/conti_init2_carto_",s,"self_",q,"QTLs.bin_200000_1.sat",sep=""))
    #On ne garde que les adultes
    temp<-temp[which(temp$age==2),]
    #Et on efface les fichiers extraits
    unlink("INIT2",recursive = TRUE)
    #Définition des noms de paramètres
    Model <- rep(paste(m,"equilibrium",sep=""),length(temp[,1])) ; Selfing <- rep(s,length(temp[,1])) ; QTLs <- rep(q,length(temp[,1])) ; Fecundity <- rep("up to capacity",length(temp[,1]))
    Selection <- rep("minimal",length(temp[,1])) ; Cas <- rep(NA,length(temp[,1])) ; Repetition <- rep(NA,length(temp[,1])) ; Environment <- rep(NA,length(temp[,1]))
    Population <- rep(NA,length(temp[,1]))
    
    #On enlève les colonnes qui ne nous intéressent pas
    temp <- temp[c(seq(1,10001,20))]
    #On ajoute les valeurs des paramètres à temp
    temp <- cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,Repetition,Population,Environment,temp)
    #On teste si toutes les pops sont occupées
    if(length(unique(temp$pop))==24)
    {
      #Si oui, on met à jour les infos de Population et d'Environnement
      for (i in 1:3) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-1}
      for (i in 4:6) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-1}
      for (i in 7:9) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-2}
      for (i in 10:12) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-2}
      for (i in 13:15) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-3}
      for (i in 16:18) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-3}
      for (i in 19:21) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-4}
      for (i in 22:24) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-4}
    }
    
    if (unique(is.na(temp$Population)==FALSE))
    {
      setwd("~/nemo_cluster/Sims2Conti")
      var_comp<-varcomp.glob(data.frame(temp$Population,temp$pop),(temp[,-c(1,10)])) ;
      save(file=paste("VarCompForFSTsEq_neutralConti/var_comp_equilibrium_",s,"self_",q,"QTLs_neutralConti.rda",sep=""),var_comp) 
    }
    Model <- m ; Selfing <- s ; QTLs <- q ; Fecundity <- f ; Selection <- w ; Cas <- c ; Repetition <- r
    #et on fait du ménage
    rm(list = c("var_comp","Model", "Selfing", "QTLs", "Fecundity", "Selection", "Cas", "Repetition", "Environment", "Population","temp","v"))
    setwd("~/nemo_cluster")
  }
}
setwd("~/nemo_cluster/Sims2Conti")
