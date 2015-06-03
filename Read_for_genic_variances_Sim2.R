#définit le répertoire de travail
setwd("~/nemo_cluster/Sims2")

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
            untar(tarfile = paste("Sim2_Model",m,"_",s,"self_",q,"QTLs_fitperpop/quanti.tar.gz",sep=""),files = paste("quanti/Sim2_Model",m,"_",s,"self_",q,"QTLs_Cas",c,"_SelInt",w,"_fec",f,"_200_",r,".quanti",sep=""),compressed = "gzip",exdir = paste("Sim2_Model",m,"_",s,"self_",q,"QTLs_fitperpop",sep=""))
            temp <- read.table(paste("Sim2_Model",m,"_",s,"self_",q,"QTLs_fitperpop/quanti/Sim2_Model",m,"_",s,"self_",q,"QTLs_Cas",c,"_SelInt",w,"_fec",f,"_200_",r,".quanti",sep=""),header=TRUE)
            unlink(paste("Sim2_Model",m,"_",s,"self_",q,"QTLs_fitperpop/quanti",sep=""),recursive = TRUE)
#Définition des noms de paramètres
            Model <- rep(m,length(temp[,1])) ; Selfing <- rep(s,length(temp[,1])) ; QTLs <- rep(q,length(temp[,1])) ; Fecundity <- rep(f,length(temp[,1]))
            Selection <- rep(w,length(temp[,1])) ; Cas <- rep(c,length(temp[,1])) ; Repetition <- rep(r,length(temp[,1])) ; Environment <- rep(NA,length(temp[,1]))
            Population <- rep(NA,length(temp[,1]))
            temp_var<-NULL
#On ajoute les valeurs des paramètres à temp
          temp <- cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,Repetition,Population,Environment,temp) 

#On teste si toutes les pops sont occupées
          ifelse(length(unique(temp$pop))==24,
                 #Si oui, on met à jour les infos de Population et d'Environnement
                 {
                 for (i in 1:3) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-1}
                 for (i in 4:6) {temp[temp$pop==i,]$Environment <-2 ; temp[temp$pop==i,]$Population <-1}
                 for (i in 7:9) {temp[temp$pop==i,]$Environment <-3 ; temp[temp$pop==i,]$Population <-2}
                 for (i in 10:12) {temp[temp$pop==i,]$Environment <-4 ; temp[temp$pop==i,]$Population <-2}
                 for (i in 13:15) {temp[temp$pop==i,]$Environment <-5 ; temp[temp$pop==i,]$Population <-3}
                 for (i in 16:18) {temp[temp$pop==i,]$Environment <-6 ; temp[temp$pop==i,]$Population <-3}
                 for (i in 19:21) {temp[temp$pop==i,]$Environment <-7 ; temp[temp$pop==i,]$Population <-4}
                 for (i in 22:24) {temp[temp$pop==i,]$Environment <-8 ; temp[temp$pop==i,]$Population <-4}
                 
                 ##########Calcul des variances pour chaque marqueur
                 res_temp_var<-NULL
                 for (j in 1:q) {
                 #Création des valeurs géniques
                   temp_bis <- temp
                   temp_bis[["l"]] <- temp_bis[[paste("t1l",j,"1",sep="")]] + temp_bis[[paste("t1l",j,"2",sep="")]]
                   
                   #Vérifie si la variance n'est pas nulle
                   ifelse (var(temp_bis$l)!=0,{
                      
                #Création de l'objet qui accueillera les résultats du calcul de variance
                 temp_var<-try(lmer(l~1|Population/Environment/pop,data = temp_bis,na.action = na.fail),silent=TRUE)
                #Si temp_var ne retourne pas d'erreur
                ifelse(class(temp_var)!= "try-error",
                  {
                #Création des composantes de la variance
                Genic_VarPopulation <- as.numeric(VarCorr(temp_var)$Population[1])
                Genic_VarEnvironment <- as.numeric(VarCorr(temp_var)$Environment[1])
                Genic_VarPatch <- as.numeric( VarCorr(temp_var)$pop[1])
                Genic_VarInd <- as.numeric(sigma(temp_var)**2)
                  },
                  
                  {
                  temp_var<-try(lme(l~1,random=~1|Population/Environment/pop,data = temp_bis,na.action = na.fail),silent=TRUE)
                  ifelse (class(temp_var)!="try-error",
                  {
                    Genic_VarPopulation<-as.numeric(VarCorr(temp_var)[2])
                    Genic_VarEnvironment<-as.numeric(VarCorr(temp_var)[4])
                    Genic_VarPatch<-as.numeric(VarCorr(temp_var)[6])
                    Genic_VarInd<-as.numeric(VarCorr(temp_var)[7])
                  },
                  {
                    #On crée des objets avec NA pour le cas ni lme ni lmer ne fonctionnent
                    Genic_VarPopulation <- NA ; Genic_VarEnvironment <- NA ; Genic_VarPatch <- NA ; Genic_VarInd <- NA
                    
                  })
                })
                   }
                , #Sinon on mets des 0
                {Genic_VarPopulation <- 0
                 Genic_VarEnvironment <- 0
                 Genic_VarPatch <- 0
                 Genic_VarInd <- 0}
                )
                 res_temp_var_bis <-as.data.frame(cbind(Genic_VarPopulation,Genic_VarEnvironment,Genic_VarPatch,Genic_VarInd), row.names = as.character(paste("L",j,sep="")))
                 res_temp_var<-rbind(res_temp_var_bis,res_temp_var)
                  }
                 #On sauve cet objet
                write.table(x = res_temp_var,file = paste("Genic_variances_components/Model",m,"_Selfing",s,"_QTLs",q,"_Fecundity",f,"_Selection",w,"_Cas",c,"_repetition",r,".txt",sep=""))
                
                 #Et on calcule la somme des variances géniques pour chaque niveau
                Genic_VarPopulation <- sum(res_temp_var$Genic_VarPopulation)
                Genic_VarEnvironment <- sum(res_temp_var$Genic_VarEnvironment)
                Genic_VarPatch <- sum(res_temp_var$Genic_VarPatch)
                Genic_VarInd <- sum(res_temp_var$Genic_VarInd)
                 },
                 #Sinon on mets des "NA" dans les objets
                {Genic_VarPopulation <- NA ; Genic_VarEnvironment <- NA ; Genic_VarPatch <- NA ; Genic_VarInd <- NA})

                  #On combine le tout
                Model <- m ; Selfing <- s ; QTLs <- q ; Fecundity <- f ; Selection <- w ; Cas <- c ; Repetition <- r
                res_temp<-as.data.frame(cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,Repetition,Genic_VarPopulation,Genic_VarEnvironment,Genic_VarPatch,Genic_VarInd))
                results <- rbind(res_temp,results)
              }
    #et on fait du ménage
    rm(list = c("Model", "Selfing", "QTLs", "Fecundity", "Selection", "Cas", "Repetition", "Environment", "Population","temp_var","temp", "Genic_VarPopulation", "Genic_VarEnvironment", "Genic_VarPatch", "Genic_VarInd","res_temp_var","res_temp_var_bis","temp_bis"))
            }
          }
        }
      }
    }
  }

#On sauve le tout
write.table(x = results,file = paste("Genic_variances_components/main/normal",sep=""))

#####Et la même pour les modèles "nuls"#####.
rm(list=ls())
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
              untar(tarfile = paste("Sim2_Model",m,"nul_",s,"self_",q,"QTLs_fitperpop/quanti.tar.gz",sep=""),files = paste("quanti/Sim2_Model",m,"nul_",s,"self_",q,"QTLs_Cas",c,"_SelInt",w,"_fec",f,"_200_",r,".quanti",sep=""),compressed = "gzip",exdir = paste("Sim2_Model",m,"nul_",s,"self_",q,"QTLs_fitperpop",sep=""))
              temp <- read.table(paste("Sim2_Model",m,"nul_",s,"self_",q,"QTLs_fitperpop/quanti/Sim2_Model",m,"nul_",s,"self_",q,"QTLs_Cas",c,"_SelInt",w,"_fec",f,"_200_",r,".quanti",sep=""),header=TRUE)
              unlink(paste("Sim2_Model",m,"nul_",s,"self_",q,"QTLs_fitperpop/quanti",sep=""),recursive = TRUE)
              #Définition des noms de paramètres
              Model <- rep(paste(m,"nul",sep=""),length(temp[,1])) ; Selfing <- rep(s,length(temp[,1])) ; QTLs <- rep(q,length(temp[,1])) ; Fecundity <- rep(f,length(temp[,1]))
              Selection <- rep(w,length(temp[,1])) ; Cas <- rep(c,length(temp[,1])) ; Repetition <- rep(r,length(temp[,1])) ; Environment <- rep(NA,length(temp[,1]))
              Population <- rep(NA,length(temp[,1]))
              temp_var<-NULL
              #On ajoute les valeurs des paramètres à temp
              temp <- cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,Repetition,Population,Environment,temp) 
              
              #On teste si toutes les pops sont occupées
              ifelse(length(unique(temp$pop))==24,
                     #Si oui, on met à jour les infos de Population et d'Environnement
{
  for (i in 1:3) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-1}
  for (i in 4:6) {temp[temp$pop==i,]$Environment <-1 ; temp[temp$pop==i,]$Population <-1}
  for (i in 7:9) {temp[temp$pop==i,]$Environment <-2 ; temp[temp$pop==i,]$Population <-2}
  for (i in 10:12) {temp[temp$pop==i,]$Environment <-2 ; temp[temp$pop==i,]$Population <-2}
  for (i in 13:15) {temp[temp$pop==i,]$Environment <-3 ; temp[temp$pop==i,]$Population <-3}
  for (i in 16:18) {temp[temp$pop==i,]$Environment <-3 ; temp[temp$pop==i,]$Population <-3}
  for (i in 19:21) {temp[temp$pop==i,]$Environment <-4 ; temp[temp$pop==i,]$Population <-4}
  for (i in 22:24) {temp[temp$pop==i,]$Environment <-4 ; temp[temp$pop==i,]$Population <-4}
  
  ##########Calcul des variances pour chaque marqueur
  res_temp_var<-NULL
  for (j in 1:q) {
    #Création des valeurs géniques
    temp_bis <- temp
    temp_bis[["l"]] <- temp_bis[[paste("t1l",j,"1",sep="")]] + temp_bis[[paste("t1l",j,"2",sep="")]]
    
    #Vérifie si la variance n'est pas nulle
    ifelse (var(temp_bis$l)!=0,{
      
      #Création de l'objet qui accueillera les résultats du calcul de variance
      temp_var<-try(lmer(l~1|Population/pop,data = temp_bis,na.action = na.fail),silent=TRUE)
      #Si temp_var ne retourne pas d'erreur
      ifelse(class(temp_var)!= "try-error",
{
  #Création des composantes de la variance
  Genic_VarPopulation <- as.numeric(VarCorr(temp_var)$Population[1])
  Genic_VarEnvironment <- NA
  Genic_VarPatch <- as.numeric( VarCorr(temp_var)$pop[1])
  Genic_VarInd <- as.numeric(sigma(temp_var)**2)
},

{
  temp_var<-try(lme(l~1,random=~1|Population/pop,data = temp_bis,na.action = na.fail),silent=TRUE)
  ifelse (class(temp_var)!="try-error",
{
  Genic_VarPopulation<-as.numeric(VarCorr(temp_var)[2])
  Genic_VarEnvironment<-NA
  Genic_VarPatch<-as.numeric(VarCorr(temp_var)[6])
  Genic_VarInd<-as.numeric(VarCorr(temp_var)[7])
},
{
  #On crée des objets avec NA pour le cas ni lme ni lmer ne fonctionnent
  Genic_VarPopulation <- NA ; Genic_VarEnvironment <- NA ; Genic_VarPatch <- NA ; Genic_VarInd <- NA
  
})
})
    }
, #Sinon on mets des 0
{Genic_VarPopulation <- 0
 Genic_VarEnvironment <- NA
 Genic_VarPatch <- 0
 Genic_VarInd <- 0}
    )
res_temp_var_bis <-as.data.frame(cbind(Genic_VarPopulation,Genic_VarEnvironment,Genic_VarPatch,Genic_VarInd), row.names = as.character(paste("L",j,sep="")))
res_temp_var<-rbind(res_temp_var_bis,res_temp_var)
  }
#On sauve cet objet
write.table(x = res_temp_var,file = paste("Genic_variances_components/Model",m,"nul_Selfing",s,"_QTLs",q,"_Fecundity",f,"_Selection",w,"_Cas",c,"_repetition",r,".txt",sep=""))

#Et on calcule la somme des variances géniques pour chaque niveau
Genic_VarPopulation <- sum(res_temp_var$Genic_VarPopulation)
Genic_VarEnvironment <- NA
Genic_VarPatch <- sum(res_temp_var$Genic_VarPatch)
Genic_VarInd <- sum(res_temp_var$Genic_VarInd)
},
#Sinon on mets des "NA" dans les objets
{Genic_VarPopulation <- NA ; Genic_VarEnvironment <- NA ; Genic_VarPatch <- NA ; Genic_VarInd <- NA})

#On combine le tout
Model <- paste(m,"nul",sep="") ; Selfing <- s ; QTLs <- q ; Fecundity <- f ; Selection <- w ; Cas <- c ; Repetition <- r
res_temp<-as.data.frame(cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,Repetition,Genic_VarPopulation,Genic_VarEnvironment,Genic_VarPatch,Genic_VarInd))
results <- rbind(res_temp,results)
            }
#et on fait du ménage
rm(list = c("Model", "Selfing", "QTLs", "Fecundity", "Selection", "Cas", "Repetition", "Environment", "Population","temp_var","temp", "Genic_VarPopulation", "Genic_VarEnvironment", "Genic_VarPatch", "Genic_VarInd","res_temp_var","res_temp_var_bis","temp_bis"))
          }
        }
      }
    }
  }
}

#On sauve le tout
write.table(x = results,file = paste("Genic_variances_components/main/nul",sep=""))