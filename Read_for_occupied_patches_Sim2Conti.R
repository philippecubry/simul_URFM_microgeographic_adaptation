#définit le répertoire de travail
setwd("~/nemo_cluster/Sims2Conti")

##### Charge les bibliothèques nécessaires #####
library("ggplot2", lib.loc="/home/pcubry/R/x86_64-pc-linux-gnu-library/3.1")
library("Hmisc", lib.loc="/home/pcubry/R/x86_64-pc-linux-gnu-library/3.1")

##### Définit une fonction basée sur les packages ggplot2 et Hmisc pour le plot de la moyenne et de l'écart-type #####
stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data=fun, geom=geom, width=0.2, ...)
}

##### Définit les variables utilisées dans les simulations et le nom des fichiers #####

model <- c("2","3","2nul","3nul") #On ne considère plus le modèle 1 qui est un cas intermédiaire entre les autres et présente peu d'intérêt
cas <- c("1","2","3","4")
fecundity <- c("010","100","200")
selfing <- c("0","03")
qtls <- c("10","50")
selection <- c("01","05","20","50")
replicate <- c(1:10)
patch <- c(1:24)

#Création d'un objet vide qui va accueillir les résultats
results = NULL

##### Boucle de lecture des fichiers #####
for(m in model) {
  for(s in selfing) {
    for(q in qtls) {
      for(f in fecundity) {
        for(w in selection) {
          for(c in cas) {
            temp <- read.delim(file=(paste("Sim2_conti_Model",m,"_",s,"self_",q,"QTLs_fitperpop","/stats/Sim2_conti_Model",m,"_",s,"self_",q,"QTLs_Cas",c,"_SelInt",w,"_fec",f,".txt",sep = "")),header=T)
            temp<-temp[which(temp$generation==200),]
            #Définition des noms de paramètres
            Model <- rep(paste(m,"_conti",sep=""),length(temp[,1]))
            Selfing <- rep(s,length(temp[,1]))
            QTLs <- rep(q,length(temp[,1]))
            Fecundity <- rep(f,length(temp[,1]))
            Selection <- rep(w,length(temp[,1]))
            Cas <- rep(c,length(temp[,1]))
            #On ajoute les valeurs des paramètres à temp
            temp <- cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,temp)
            temp_occupied = NULL
            for (r in replicate){
              nb_occupied_patches <- 0
              temp_r <- temp[which(temp$replicate==r),]
              for (p in patch){
                test<-try(if (temp_r[,eval(parse(text=paste('"adlt.fem.p',p,'"',sep="")))] != 0){},silent=TRUE)
                if (class(test)!= "try-error"){
                if (temp_r[,eval(parse(text=paste('"adlt.fem.p',p,'"',sep="")))] != 0) {nb_occupied_patches <- nb_occupied_patches + 1}
                } else {
                  temp_r[1,1:length(temp_r)]<-NA ; temp_r$Model <- paste(m,"_conti",sep=""); temp_r$Selfing <- s ; temp_r$QTLs <- q; temp_r$Fecundity <- f ; temp_r$Selection <- w ; temp_r$Cas <- c; temp_r$replicate <- r
                }
            } 
            temp_occupied<-rbind(temp_occupied,cbind(temp_r,nb_occupied_patches))
            nb_occupied_patches <- 0
            }
                        
            #et on ajoute temp aux résultats !!!
            results <- rbind(temp_occupied,results)
          }  
        }
      }
    }
  }
}

results_sauve <- results

##### On réordonne les cas 1, 2, 3, 4 en 2, 1, 3, 4 et on les renomme, on réordonne les sélections en 50,20,5,1 et les QTLs en 10,50 #####
results$Cas<-factor(results$Cas,levels=c("2","1","3","4")) ; levels(results$Cas) <- c("A","B","C","D")
results$Selection<-factor(results$Selection,levels=c("50","20","05","01"))
results$QTLs<-factor(results$QTLs,levels=c("10","50"))
results$Selfing<-factor(results$Selfing,levels=c("0","03")) ; levels(results$Selfing) <- c("0","0.3")
results$Model<-factor(results$Model,levels=c("2_conti","2nul_conti","3_conti","3nul_conti")) ; levels(results$Model) <- c("Overlapping with within population environmental heterogeneity","Overlapping without within population environmental heterogeneity","Divided with within population environmental heterogeneity","Divided without within population environmental heterogeneity")

Case <- c('A'='Case A','B'='Case B','C'='Case C','D'='Case D'); Fecun <- c("010"="Mean Fecundity = 10","100"="Mean Fecundity = 100","200"="Mean Fecundity = 200")

##### Création des graphes d'occupation des patches #####
for (s in c("0","0.3")) {
  g <- ggplot(results[which(results$Selfing == s & results$Fecundity == "100"),])
  h <- g + facet_grid(Fecundity+QTLs~Cas, labeller = labeller(Fecundity=Fecun, Cas = Case,QTLs=label_both)) + aes(Selection,nb_occupied_patches,shape = Model)+ labs(y="Number of Occupied Patches", title = paste("Continuous Mutation Model, Selfing Rate of",s)) + theme_bw(base_size=15)
  h <- h +  scale_shape_manual(values = c(15,0,16,1)) +
    guides(shape=guide_legend(direction = "horizontal",title.position = "top", title.hjust=0.5, title.theme = element_text(size = 15,angle = 0),label.position = "bottom", label.hjust=0.5,label.theme = element_text(angle = 90))) + scale_y_continuous(breaks=c(seq(0,24,3))) + theme(panel.grid.minor=element_blank()) +
    stat_summary(geom="pointrange", fun.y = median, fun.ymin = min, fun.ymax=max, position = position_dodge(width = 0.75), size = 0.8) +
    scale_y_continuous(breaks=c(seq(0,24,3))) +
    theme(panel.grid.minor=element_blank())
 
  #Ecriture du pdf
  pdf(file = paste("ContinuousModel_Selfing",s,"_conti_nb_occupied_patches_fec100.pdf",sep=""), height = 12, width = 16, paper = "special")
  print(h)
  dev.off() 
    
}

for (s in c("0","0.3")) {
  g <- ggplot(results[which(results$Selfing == s & results$Fecundity != "100"),])
  h <- g + facet_grid(Fecundity+QTLs~Cas, labeller = labeller(Fecundity=Fecun, Cas = Case,QTLs=label_both)) + aes(Selection,nb_occupied_patches,shape = Model)+ labs(y="Number of Occupied Patches", title = paste("Continuous Mutation Model, Selfing Rate of",s)) + theme_bw(base_size=15)
  h <- h +  scale_shape_manual(values = c(15,0,16,1)) +
    guides(shape=guide_legend(direction = "horizontal",title.position = "top", title.hjust=0.5, title.theme = element_text(size = 15,angle = 0),label.position = "bottom", label.hjust=0.5,label.theme = element_text(angle = 90))) + scale_y_continuous(breaks=c(seq(0,24,3))) + theme(panel.grid.minor=element_blank()) +
    stat_summary(geom="pointrange", fun.y = median, fun.ymin = min, fun.ymax=max, position = position_dodge(width = 0.75), size = 0.8) +
    scale_y_continuous(breaks=c(seq(0,24,3))) +
    theme(panel.grid.minor=element_blank())
  
  #Ecriture du pdf
  pdf(file = paste("ContinuousModel_Selfing",s,"_conti_nb_occupied_patches_fec10-200.pdf",sep=""), height = 12, width = 16, paper = "special")
  print(h)
  dev.off() 
  
}

##### Definition des optimums locaux demandés #####
results$Zopt_Env1 <- NA;results$Zopt_Env2 <- NA;results$Zopt_Env3 <- NA;results$Zopt_Env4 <- NA;results$Zopt_Env5 <- NA;results$Zopt_Env6 <- NA;results$Zopt_Env7 <- NA;results$Zopt_Env8 <- NA

Z_Model2_CasB_10QTLs <- c(-13.49,-4.49666666666667,-7.4944444444,1.49888888888889,13.49,4.49666666666667,7.4944444444,-1.49888888888889)[c(1,2,3,4,8,7,6,5)]
Z_Model2_CasA_10QTLs <- c(-16.8625,-5.62083333333333,-9.36805555555556,1.87361111111111,16.8625,5.62083333333333,9.36805555555556,-1.87361111111111)[c(1,2,3,4,8,7,6,5)]
Z_Model2_CasC_10QTLs <- c(-6.745,-2.24833333333333,-3.74722222222222,0.749444444444444,6.745,2.24833333333333,3.74722222222222,-0.749444444444444)[c(1,2,3,4,8,7,6,5)]
Z_Model2_CasD_10QTLs <- c(-3.3725,-1.1241666667,-1.87361111111111,0.374722222222222,3.3725,1.1241666667,1.87361111111111,-0.374722222222222)[c(1,2,3,4,8,7,6,5)]
Z_Model2_CasB_50QTLs <- c(-67.45,-22.4833333333333,-37.4722222222222,7.49444444444444,67.45,22.4833333333333,37.4722222222222,-7.49444444444444)[c(1,2,3,4,8,7,6,5)]
Z_Model2_CasA_50QTLs <- c(-84.3125,-28.1041666666667,-46.8402777777778,9.36805555555556,84.3125,28.1041666666667,46.8402777777778,-9.36805555555556)[c(1,2,3,4,8,7,6,5)]
Z_Model2_CasC_50QTLs <- c(-33.725,-11.2416666666667,-18.7361111111111,3.74722222222222,33.725,11.2416666666667,18.7361111111111,-3.74722222222222)[c(1,2,3,4,8,7,6,5)]
Z_Model2_CasD_50QTLs <- c(-16.8625,-5.62083333333333,-9.36805555555556,1.87361111111111,16.8625,5.62083333333333,9.36805555555556,-1.87361111111111)[c(1,2,3,4,8,7,6,5)]

Z_Model3_CasB_10QTLs <- c(-13.49,-9.63571428571429,-5.78142857142857,-1.92714285714286,13.49,9.63571428571429,5.78142857142857,1.92714285714286)[c(1,2,3,4,8,7,6,5)]
Z_Model3_CasA_10QTLs <- c(-16.8625,-12.0446428571429,-7.22678571428572,-2.40892857142857,16.8625,12.0446428571429,7.22678571428572,2.40892857142857)[c(1,2,3,4,8,7,6,5)]
Z_Model3_CasC_10QTLs <- c(-6.745,-4.81785714285714,-2.89071428571429,-0.963571428571429,6.745,4.81785714285714,2.89071428571429,0.963571428571429)[c(1,2,3,4,8,7,6,5)]
Z_Model3_CasD_10QTLs <- c(-3.3725,-2.40892857142857,-1.44535714285714,-0.481785714285714,3.3725,2.40892857142857,1.44535714285714,0.481785714285714)[c(1,2,3,4,8,7,6,5)]
Z_Model3_CasB_50QTLs <- c(-67.45,-48.1785714285714,-28.9071428571429,-9.63571428571429,67.45,48.1785714285714,28.9071428571429,9.63571428571429)[c(1,2,3,4,8,7,6,5)]
Z_Model3_CasA_50QTLs <- c(-84.3125,-60.2232142857143,-36.1339285714286,-12.0446428571429,84.3125,60.2232142857143,36.1339285714286,12.0446428571429)[c(1,2,3,4,8,7,6,5)]
Z_Model3_CasC_50QTLs <- c(-33.725,-24.0892857142857,-14.4535714285714,-4.81785714285714,33.725,24.0892857142857,14.4535714285714,4.81785714285714)[c(1,2,3,4,8,7,6,5)]
Z_Model3_CasD_50QTLs <- c(-16.8625,-12.0446428571429,-7.22678571428572,-2.40892857142857,16.8625,12.0446428571429,7.22678571428572,2.40892857142857)[c(1,2,3,4,8,7,6,5)]

Z_Model2nul_CasB_10QTLs <- c(-8.99333333333333  , -8.99333333333333  ,  -2.99777777777778	, -2.99777777777778  ,	2.99777777777778	,  2.99777777777778	,	8.99333333333333,  8.99333333333333)
Z_Model2nul_CasA_10QTLs <- c(-11.2416666666667  ,	-11.2416666666667  ,  -3.74722222222222	, -3.74722222222222  ,	3.74722222222222	,  3.74722222222222	,	11.2416666666667,  11.2416666666667)
Z_Model2nul_CasC_10QTLs <- c(-4.49666666666667  ,-4.49666666666667  ,	-1.49888888888889	,  -1.49888888888889	,	1.49888888888889	,  1.49888888888889	,	4.49666666666667,  4.49666666666667)
Z_Model2nul_CasD_10QTLs <- c(-2.24833333333333  ,-2.24833333333333  ,	-0.749444444444444	,  -0.749444444444444	,	0.749444444444444	,  0.749444444444444	,	2.24833333333333,  2.24833333333333)
Z_Model2nul_CasB_50QTLs <- c(-44.9666666666667  ,	-44.9666666666667  ,  -14.9888888888889	,-14.9888888888889  ,	14.9888888888889	,  14.9888888888889	,	44.9666666666667,  44.9666666666667)
Z_Model2nul_CasA_50QTLs <- c(-56.2083333333333  ,	-56.2083333333333  ,  -18.7361111111111	,-18.7361111111111  ,	18.7361111111111	,  18.7361111111111	,	56.2083333333333,  56.2083333333333)
Z_Model2nul_CasC_50QTLs <- c(-22.4833333333333  ,-22.4833333333333  ,	-7.49444444444444	,  -7.49444444444444	,	7.49444444444444	,  7.49444444444444	,	22.4833333333333,  22.4833333333333)
Z_Model2nul_CasD_50QTLs <- c(-11.2416666666667  ,-11.2416666666667  ,	-3.74722222222222	,  -3.74722222222222	,	3.74722222222222	,  3.74722222222222	,	11.2416666666667,  11.2416666666667)

Z_Model3nul_CasB_10QTLs <- c(-11.5628571428571  ,  -11.5628571428571  ,  -3.85428571428571	,-3.85428571428571  ,	3.85428571428571	,3.85428571428571  ,	11.5628571428571,  11.5628571428571)
Z_Model3nul_CasA_10QTLs <- c(-14.4535714285714  ,	-14.4535714285714  ,  -4.81785714285714	,-4.81785714285714  ,	4.81785714285714	,  4.81785714285714	,	14.4535714285714,  14.4535714285714)
Z_Model3nul_CasC_10QTLs <- c(-5.78142857142857  ,	-5.78142857142857  ,  -1.92714285714286	,-1.92714285714286  ,	1.92714285714286	,  1.92714285714286	,	5.78142857142857,  5.78142857142857)
Z_Model3nul_CasD_10QTLs <- c(-2.89071428571429  ,	-2.89071428571429  ,  -0.963571428571428	,-0.963571428571428  ,	0.963571428571428	,  0.963571428571428	,	2.89071428571429,  2.89071428571429)
Z_Model3nul_CasB_50QTLs <- c(-57.8142857142857  ,-57.8142857142857  ,	-19.2714285714286	,  -19.2714285714286	,	19.2714285714286	,  19.2714285714286	,	57.8142857142857,  57.8142857142857)
Z_Model3nul_CasA_50QTLs <- c(-72.2678571428571  ,-72.2678571428571  ,	-24.0892857142857	,  -24.0892857142857	,	24.0892857142857	,  24.0892857142857	,	72.2678571428571,  72.2678571428571)
Z_Model3nul_CasC_50QTLs <- c(-28.9071428571429  ,-28.9071428571429  ,	-9.63571428571429	,  -9.63571428571429	,	9.63571428571429	,  9.63571428571429	,	28.9071428571429,  28.9071428571429)
Z_Model3nul_CasD_50QTLs <- c(-14.4535714285714  ,-14.4535714285714  ,	-4.81785714285714	,  -4.81785714285714	,	4.81785714285714	,  4.81785714285714	,	14.4535714285714,  14.4535714285714)

for (i in (1:8)){
  results[results$Cas=="B" & results$QTLs=="10" & results$Model == "Overlapping with within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model2_CasB_10QTLs[i]
  results[results$Cas=="A" & results$QTLs=="10" & results$Model == "Overlapping with within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model2_CasA_10QTLs[i]
  results[results$Cas=="C" & results$QTLs=="10" & results$Model == "Overlapping with within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model2_CasC_10QTLs[i]
  results[results$Cas=="D" & results$QTLs=="10" & results$Model == "Overlapping with within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model2_CasD_10QTLs[i]
  
  results[results$Cas=="B" & results$QTLs=="50" & results$Model == "Overlapping with within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model2_CasB_50QTLs[i]
  results[results$Cas=="A" & results$QTLs=="50" & results$Model == "Overlapping with within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model2_CasA_50QTLs[i]
  results[results$Cas=="C" & results$QTLs=="50" & results$Model == "Overlapping with within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model2_CasC_50QTLs[i]
  results[results$Cas=="D" & results$QTLs=="50" & results$Model == "Overlapping with within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model2_CasD_50QTLs[i]
  
  results[results$Cas=="B" & results$QTLs=="10" & results$Model == "Divided with within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model3_CasB_10QTLs[i]
  results[results$Cas=="A" & results$QTLs=="10" & results$Model == "Divided with within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model3_CasA_10QTLs[i]
  results[results$Cas=="C" & results$QTLs=="10" & results$Model == "Divided with within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model3_CasC_10QTLs[i]
  results[results$Cas=="D" & results$QTLs=="10" & results$Model == "Divided with within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model3_CasD_10QTLs[i]
  
  results[results$Cas=="B" & results$QTLs=="50" & results$Model == "Divided with within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model3_CasB_50QTLs[i]
  results[results$Cas=="A" & results$QTLs=="50" & results$Model == "Divided with within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model3_CasA_50QTLs[i]
  results[results$Cas=="C" & results$QTLs=="50" & results$Model == "Divided with within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model3_CasC_50QTLs[i]
  results[results$Cas=="D" & results$QTLs=="50" & results$Model == "Divided with within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model3_CasD_50QTLs[i]
  
  results[results$Cas=="B" & results$QTLs=="10" & results$Model == "Overlapping without within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model2nul_CasB_10QTLs[i]
  results[results$Cas=="A" & results$QTLs=="10" & results$Model == "Overlapping without within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model2nul_CasA_10QTLs[i]
  results[results$Cas=="C" & results$QTLs=="10" & results$Model == "Overlapping without within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model2nul_CasC_10QTLs[i]
  results[results$Cas=="D" & results$QTLs=="10" & results$Model == "Overlapping without within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model2nul_CasD_10QTLs[i]
  
  results[results$Cas=="B" & results$QTLs=="50" & results$Model == "Overlapping without within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model2nul_CasB_50QTLs[i]
  results[results$Cas=="A" & results$QTLs=="50" & results$Model == "Overlapping without within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model2nul_CasA_50QTLs[i]
  results[results$Cas=="C" & results$QTLs=="50" & results$Model == "Overlapping without within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model2nul_CasC_50QTLs[i]
  results[results$Cas=="D" & results$QTLs=="50" & results$Model == "Overlapping without within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model2nul_CasD_50QTLs[i]
  
  results[results$Cas=="B" & results$QTLs=="10" & results$Model == "Divided without within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model3nul_CasB_10QTLs[i]
  results[results$Cas=="A" & results$QTLs=="10" & results$Model == "Divided without within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model3nul_CasA_10QTLs[i]
  results[results$Cas=="C" & results$QTLs=="10" & results$Model == "Divided without within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model3nul_CasC_10QTLs[i]
  results[results$Cas=="D" & results$QTLs=="10" & results$Model == "Divided without within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model3nul_CasD_10QTLs[i]
  
  results[results$Cas=="B" & results$QTLs=="50" & results$Model == "Divided without within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model3nul_CasB_50QTLs[i]
  results[results$Cas=="A" & results$QTLs=="50" & results$Model == "Divided without within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model3nul_CasA_50QTLs[i]
  results[results$Cas=="C" & results$QTLs=="50" & results$Model == "Divided without within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model3nul_CasC_50QTLs[i]
  results[results$Cas=="D" & results$QTLs=="50" & results$Model == "Divided without within population environmental heterogeneity",][paste("Zopt_Env",i,sep="")] <- Z_Model3nul_CasD_50QTLs[i]
}

#Par populations
results$Zopt_Pop1 <- (results$Zopt_Env1+results$Zopt_Env2)/2
results$Zopt_Pop2 <- (results$Zopt_Env3+results$Zopt_Env4)/2
results$Zopt_Pop3 <- (results$Zopt_Env5+results$Zopt_Env6)/2
results$Zopt_Pop4 <- (results$Zopt_Env7+results$Zopt_Env8)/2

#On rajoutte les espérances/max de la valeur génotypique
results$MaxG <- NA
results$MaxG[(which(results$QTLs == 10))] <- 13.49*2
results$MaxG[(which(results$QTLs == 50))] <- 67.45*2

##### Création des statistiques de phénotypes réalisés par environnement/population #####
results$Env1 <- ((results$adlt.q1.p1*results$adlt.fem.p1+results$adlt.q1.p2*results$adlt.fem.p2+results$adlt.q1.p3*results$adlt.fem.p3)/(results$adlt.fem.p1+results$adlt.fem.p2+results$adlt.fem.p3))
results$Env2 <- ((results$adlt.q1.p4*results$adlt.fem.p4+results$adlt.q1.p5*results$adlt.fem.p5+results$adlt.q1.p6*results$adlt.fem.p6)/(results$adlt.fem.p4+results$adlt.fem.p5+results$adlt.fem.p6))
results$Env3 <- ((results$adlt.q1.p7*results$adlt.fem.p7+results$adlt.q1.p8*results$adlt.fem.p8+results$adlt.q1.p9*results$adlt.fem.p9)/(results$adlt.fem.p7+results$adlt.fem.p8+results$adlt.fem.p9))
results$Env4 <- ((results$adlt.q1.p10*results$adlt.fem.p10+results$adlt.q1.p11*results$adlt.fem.p11+results$adlt.q1.p12*results$adlt.fem.p12)/(results$adlt.fem.p10+results$adlt.fem.p11+results$adlt.fem.p12))
results$Env5 <- ((results$adlt.q1.p13*results$adlt.fem.p13+results$adlt.q1.p14*results$adlt.fem.p14+results$adlt.q1.p15*results$adlt.fem.p15)/(results$adlt.fem.p13+results$adlt.fem.p14+results$adlt.fem.p15))
results$Env6 <- ((results$adlt.q1.p16*results$adlt.fem.p16+results$adlt.q1.p17*results$adlt.fem.p17+results$adlt.q1.p18*results$adlt.fem.p18)/(results$adlt.fem.p16+results$adlt.fem.p17+results$adlt.fem.p18))
results$Env7 <- ((results$adlt.q1.p19*results$adlt.fem.p19+results$adlt.q1.p20*results$adlt.fem.p20+results$adlt.q1.p21*results$adlt.fem.p21)/(results$adlt.fem.p19+results$adlt.fem.p20+results$adlt.fem.p21))
results$Env8 <- ((results$adlt.q1.p22*results$adlt.fem.p22+results$adlt.q1.p23*results$adlt.fem.p23+results$adlt.q1.p24*results$adlt.fem.p24)/(results$adlt.fem.p22+results$adlt.fem.p23+results$adlt.fem.p24))

results$Pop1 <- ((results$adlt.q1.p1*results$adlt.fem.p1+results$adlt.q1.p2*results$adlt.fem.p2+results$adlt.q1.p3*results$adlt.fem.p3+results$adlt.q1.p4*results$adlt.fem.p4+results$adlt.q1.p5*results$adlt.fem.p5+results$adlt.q1.p6*results$adlt.fem.p6)/(results$adlt.fem.p1+results$adlt.fem.p2+results$adlt.fem.p3+results$adlt.fem.p4+results$adlt.fem.p5+results$adlt.fem.p6))
results$Pop2 <- ((results$adlt.q1.p7*results$adlt.fem.p7+results$adlt.q1.p8*results$adlt.fem.p8+results$adlt.q1.p9*results$adlt.fem.p9+results$adlt.q1.p10*results$adlt.fem.p10+results$adlt.q1.p11*results$adlt.fem.p11+results$adlt.q1.p12*results$adlt.fem.p12)/(results$adlt.fem.p7+results$adlt.fem.p8+results$adlt.fem.p9+results$adlt.fem.p10+results$adlt.fem.p11+results$adlt.fem.p12))
results$Pop3 <- ((results$adlt.q1.p13*results$adlt.fem.p13+results$adlt.q1.p14*results$adlt.fem.p14+results$adlt.q1.p15*results$adlt.fem.p15+results$adlt.q1.p16*results$adlt.fem.p16+results$adlt.q1.p17*results$adlt.fem.p17+results$adlt.q1.p18*results$adlt.fem.p18)/(results$adlt.fem.p13+results$adlt.fem.p14+results$adlt.fem.p15+results$adlt.fem.p16+results$adlt.fem.p17+results$adlt.fem.p18))
results$Pop4 <- ((results$adlt.q1.p19*results$adlt.fem.p19+results$adlt.q1.p20*results$adlt.fem.p20+results$adlt.q1.p21*results$adlt.fem.p21+results$adlt.q1.p22*results$adlt.fem.p22+results$adlt.q1.p23*results$adlt.fem.p23+results$adlt.q1.p24*results$adlt.fem.p24)/(results$adlt.fem.p19+results$adlt.fem.p20+results$adlt.fem.p21+results$adlt.fem.p22+results$adlt.fem.p23+results$adlt.fem.p24))

##### Calcul des écarts aux Zopt demandé par environnement ou par pop #####
results$LagEnv = (abs(results$adlt.q1.p1 - results$Zopt_Env1)/results$MaxG+abs(results$adlt.q1.p2 - results$Zopt_Env1)/results$MaxG+abs(results$adlt.q1.p3 - results$Zopt_Env1)/results$MaxG+
                    abs(results$adlt.q1.p4 - results$Zopt_Env2)/results$MaxG+abs(results$adlt.q1.p5 - results$Zopt_Env2)/results$MaxG+abs(results$adlt.q1.p6 - results$Zopt_Env2)/results$MaxG+
                    abs(results$adlt.q1.p7 - results$Zopt_Env3)/results$MaxG+abs(results$adlt.q1.p8 - results$Zopt_Env3)/results$MaxG+abs(results$adlt.q1.p9 - results$Zopt_Env3)/results$MaxG+
                    abs(results$adlt.q1.p10 - results$Zopt_Env4)/results$MaxG+abs(results$adlt.q1.p11 - results$Zopt_Env4)/results$MaxG+abs(results$adlt.q1.p12 - results$Zopt_Env4)/results$MaxG+
                    abs(results$adlt.q1.p13 - results$Zopt_Env5)/results$MaxG+abs(results$adlt.q1.p14 - results$Zopt_Env5)/results$MaxG+abs(results$adlt.q1.p15 - results$Zopt_Env5)/results$MaxG+
                    abs(results$adlt.q1.p16 - results$Zopt_Env6)/results$MaxG+abs(results$adlt.q1.p17 - results$Zopt_Env6)/results$MaxG+abs(results$adlt.q1.p18 - results$Zopt_Env6)/results$MaxG+
                    abs(results$adlt.q1.p19 - results$Zopt_Env7)/results$MaxG+abs(results$adlt.q1.p20 - results$Zopt_Env7)/results$MaxG+abs(results$adlt.q1.p21 - results$Zopt_Env7)/results$MaxG+
                    abs(results$adlt.q1.p22 - results$Zopt_Env8)/results$MaxG+abs(results$adlt.q1.p23 - results$Zopt_Env8)/results$MaxG+abs(results$adlt.q1.p24 - results$Zopt_Env8)/results$MaxG)/24
results$LagPop = (abs(results$adlt.q1.p1 - results$Zopt_Pop1)/results$MaxG+abs(results$adlt.q1.p2 - results$Zopt_Pop1)/results$MaxG+abs(results$adlt.q1.p3 - results$Zopt_Pop1)/results$MaxG+
                    abs(results$adlt.q1.p4 - results$Zopt_Pop1)/results$MaxG+abs(results$adlt.q1.p5 - results$Zopt_Pop1)/results$MaxG+abs(results$adlt.q1.p6 - results$Zopt_Pop1)/results$MaxG+
                    abs(results$adlt.q1.p7 - results$Zopt_Pop2)/results$MaxG+abs(results$adlt.q1.p8 - results$Zopt_Pop2)/results$MaxG+abs(results$adlt.q1.p9 - results$Zopt_Pop2)/results$MaxG+
                    abs(results$adlt.q1.p10 - results$Zopt_Pop2)/results$MaxG+abs(results$adlt.q1.p11 - results$Zopt_Pop2)/results$MaxG+abs(results$adlt.q1.p12 - results$Zopt_Pop2)/results$MaxG+
                    abs(results$adlt.q1.p13 - results$Zopt_Pop3)/results$MaxG+abs(results$adlt.q1.p14 - results$Zopt_Pop3)/results$MaxG+abs(results$adlt.q1.p15 - results$Zopt_Pop3)/results$MaxG+
                    abs(results$adlt.q1.p16 - results$Zopt_Pop3)/results$MaxG+abs(results$adlt.q1.p17 - results$Zopt_Pop3)/results$MaxG+abs(results$adlt.q1.p18 - results$Zopt_Pop3)/results$MaxG+
                    abs(results$adlt.q1.p19 - results$Zopt_Pop4)/results$MaxG+abs(results$adlt.q1.p20 - results$Zopt_Pop4)/results$MaxG+abs(results$adlt.q1.p21 - results$Zopt_Pop4)/results$MaxG+
                    abs(results$adlt.q1.p22 - results$Zopt_Pop4)/results$MaxG+abs(results$adlt.q1.p23 - results$Zopt_Pop4)/results$MaxG+abs(results$adlt.q1.p24 - results$Zopt_Pop4)/results$MaxG)/24

##### Création des graphiques #####
#Pour le modèle 2
for (q in qtls){
  if (q=="10"){
    for (c in names(Case)){
      g <- ggplot(results[which((results$Model=="Overlapping with within population environmental heterogeneity" &
                                   results$QTLs == q & results$Cas==c |
                                   results$Model== "Overlapping without within population environmental heterogeneity" &
                                   results$QTLs == q & results$Cas==c)),])
      h <- g + aes(colour = Selection, shape = Selection) +  scale_shape_manual(values = c(1,0,2,5)) +
      facet_grid(Fecundity~Selfing,labeller=label_both) +
        labs(x= "Asked phenotypic value",y="Difference between realized and asked phenotypic value normalized by the esperance of max genotypic value", title = paste("Continuous Mutation Model, Overlapping Model, Case ",c, ", ", q, "QTLs"))  +
        geom_point(aes(Zopt_Env1,(Env1-Zopt_Env1)/(13.49*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env2,(Env2-Zopt_Env2)/(13.49*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env3,(Env3-Zopt_Env3)/(13.49*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env4,(Env4-Zopt_Env4)/(13.49*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env5,(Env5-Zopt_Env5)/(13.49*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env6,(Env6-Zopt_Env6)/(13.49*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env7,(Env7-Zopt_Env7)/(13.49*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env8,(Env8-Zopt_Env8)/(13.49*2)),alpha=0.8) +
        scale_x_continuous(breaks=c(eval(parse(text=paste("Z_Model2_Cas", c,"_",q,"QTLs",sep=""))),eval(parse(text=paste("Z_Model2nul_Cas", c,"_",q,"QTLs",sep="")))[c(1,3,5,7)]),labels = c("E 1a", "E 1b","E 2a", "E 2b","E 3a",  "E 3b","E 4a", "E 4b", "E 1", "E 2", "E 3", "E 4")) +
        theme_bw() + theme(panel.grid=element_blank()) +
        theme(axis.text.x  = element_text(angle=-90, vjust=1, hjust=0, size = 8))
      
      #       for (i in (1:8)){
      #         h <-h + geom_vline(xintercept= eval(parse(text=paste("Z_Model2_Cas", c,"_",q,"QTLs",sep="")))[i],size = 0.5, linetype = "solid",colour = "purple")
      #         h <-h + geom_vline(xintercept= eval(parse(text=paste("Z_Model2nul_Cas", c,"_",q,"QTLs",sep="")))[i],size = 0.5, linetype = "dashed", colour = "orange")
      #       }
      h <-h + geom_hline(yintercept= 0,size = 0.5, linetype = "dashed", colour = "gray")  
      h <- h+ geom_point(aes(x=c(-13.49,13.49), y=0.3), colour="gray", fill = "gray", text="Optimum",  size = 5, shape = 25,alpĥa = 0.4)
      if (c == "A"|c == "B") {h <- h + ylim(-0.3,0.3)} else {h <- h + ylim (-0.3,0.3)}
      pdf(file = paste("ContinuousModel_Overlapping_PhenotValues_Case",c,"_",q,"QTLs",sep=""), height = 12, width = 16, paper = "special")        
      print(h)
      dev.off()
      
    }
  } else {
    for (c in names(Case)){
      g <- ggplot(results[which((results$Model=="Overlapping with within population environmental heterogeneity" &
                                   results$QTLs == q & results$Cas==c |
                                   results$Model== "Overlapping without within population environmental heterogeneity" &
                                   results$QTLs == q & results$Cas==c)),])
      h <- g + aes(colour = Selection, shape = Selection) +  scale_shape_manual(values = c(1,0,2,5)) +
        facet_grid(Fecundity~Selfing,labeller=label_both) +
        labs(x= "Asked phenotypic value",y="Difference between realized and asked phenotypic value normalized by the esperance of max genotypic value", title = paste("Continuous Mutation Model, Overlapping Model, Case ",c, ", ", q, "QTLs"))  +
        geom_point(aes(Zopt_Env1,(Env1-Zopt_Env1)/(67.45*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env2,(Env2-Zopt_Env2)/(67.45*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env3,(Env3-Zopt_Env3)/(67.45*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env4,(Env4-Zopt_Env4)/(67.45*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env5,(Env5-Zopt_Env5)/(67.45*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env6,(Env6-Zopt_Env6)/(67.45*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env7,(Env7-Zopt_Env7)/(67.45*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env8,(Env8-Zopt_Env8)/(67.45*2)),alpha=0.8) +
        scale_x_continuous(breaks=c(eval(parse(text=paste("Z_Model2_Cas", c,"_",q,"QTLs",sep=""))),eval(parse(text=paste("Z_Model2nul_Cas", c,"_",q,"QTLs",sep="")))[c(1,3,5,7)]),labels = c("E 1a", "E 1b","E 2a", "E 2b","E 3a",  "E 3b","E 4a", "E 4b", "E 1", "E 2", "E 3", "E 4")) +
        theme_bw() + theme(panel.grid=element_blank()) +
        theme(axis.text.x  = element_text(angle=-90, vjust= 1,hjust=0,size = 8))
      
      #       for (i in (1:8)){
      #         h <-h + geom_vline(xintercept= eval(parse(text=paste("Z_Model2_Cas", c,"_",q,"QTLs",sep="")))[i],size = 0.5, linetype = "solid",colour = "purple")
      #         h <-h + geom_vline(xintercept= eval(parse(text=paste("Z_Model2nul_Cas", c,"_",q,"QTLs",sep="")))[i],size = 0.5, linetype = "dashed", colour = "orange")
      #           } 
      h <-h + geom_hline(yintercept= 0,size = 0.5, linetype = "dashed", colour = "gray")  
      h <- h+ geom_point(aes(x=c(-67.45,67.45), y=0.3), colour="gray", fill = "gray", text="Optimum",  size = 5, shape = 25,alpĥa = 0.4)
      if (c == "A"|c == "B") {h <- h + ylim(-0.3,0.3)} else {h <- h + ylim (-0.3,0.3)}      
      pdf(file = paste("ContinuousModel_Overlapping_PhenotValues_Case",c,"_",q,"QTLs",sep=""), height = 12, width = 16, paper = "special")        
      print(h)
      dev.off()
      
    }
  }
}

#Pour le modèle 3
for (q in qtls){
  if (q=="10"){
    for (c in names(Case)){
      g <- ggplot(results[which((results$Model=="Divided with within population environmental heterogeneity" &
                                   results$QTLs == q & results$Cas==c |
                                   results$Model== "Divided without within population environmental heterogeneity" &
                                   results$QTLs == q & results$Cas==c)),])
      h <- g + aes(colour = Selection, shape = Selection) +  scale_shape_manual(values = c(1,0,2,5)) +
        facet_grid(Fecundity~Selfing,labeller=label_both) +
        labs(x= "Asked phenotypic value",y="Difference between realized and asked phenotypic value normalized by the esperance of max genotypic value", title = paste("Continuous Mutation Model, Divided Model, Case ",c, ", ", q, "QTLs"))  +
        geom_point(aes(Zopt_Env1,(Env1-Zopt_Env1)/(13.49*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env2,(Env2-Zopt_Env2)/(13.49*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env3,(Env3-Zopt_Env3)/(13.49*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env4,(Env4-Zopt_Env4)/(13.49*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env5,(Env5-Zopt_Env5)/(13.49*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env6,(Env6-Zopt_Env6)/(13.49*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env7,(Env7-Zopt_Env7)/(13.49*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env8,(Env8-Zopt_Env8)/(13.49*2)),alpha=0.8) +
        scale_x_continuous(breaks=c(eval(parse(text=paste("Z_Model3_Cas", c,"_",q,"QTLs",sep=""))),eval(parse(text=paste("Z_Model3nul_Cas", c,"_",q,"QTLs",sep="")))[c(1,3,5,7)]),labels = c("E 1a", "E 1b","E 2a", "E 2b","E 3a",  "E 3b","E 4a", "E 4b", "E 1", "E 2", "E 3", "E 4")) +
        theme_bw() + theme(panel.grid=element_blank()) + 
        theme(axis.text.x  = element_text(angle=-90, vjust=1, hjust=0))
      
      #       for (i in (1:8)){
      #         h <-h + geom_vline(xintercept= eval(parse(text=paste("Z_Model3_Cas", c,"_",q,"QTLs",sep="")))[i],size = 0.5, linetype = "solid",colour = "purple")
      #         h <-h + geom_vline(xintercept= eval(parse(text=paste("Z_Model3nul_Cas", c,"_",q,"QTLs",sep="")))[i],size = 0.5, linetype = "dashed", colour = "orange")
      #       }
      h <-h + geom_hline(yintercept= 0,size = 0.5, linetype = "dashed", colour = "gray")  
      h <- h+ geom_point(aes(x=c(-13.49,13.49), y=0.3), colour="gray", fill = "gray", text="Optimum",  size = 5, shape = 25,alpĥa = 0.4)
      if (c == "A"|c == "B") {h <- h + ylim(-0.3,0.3)} else {h <- h + ylim (-0.3,0.3)}
      pdf(file = paste("ContinuousModel_Divided_PhenotValues_Case",c,"_",q,"QTLs",sep=""), height = 12, width = 16, paper = "special")        
      print(h)
      dev.off()
      
    }
  } else {
    for (c in names(Case)){
      g <- ggplot(results[which((results$Model=="Divided with within population environmental heterogeneity" &
                                   results$QTLs == q & results$Cas==c |
                                   results$Model== "Divided without within population environmental heterogeneity" &
                                   results$QTLs == q & results$Cas==c)),])
      h <- g + aes(colour = Selection, shape = Selection) +  scale_shape_manual(values = c(1,0,2,5)) +
        facet_grid(Fecundity~Selfing,labeller=label_both) +
        labs(x= "Asked phenotypic value",y="Difference between realized and asked phenotypic value normalized by the esperance of max genotypic value", title = paste("Continuous Mutation Model, Divided Model, Case ",c, ", ", q, "QTLs"))  +
        geom_point(aes(Zopt_Env1,(Env1-Zopt_Env1)/(67.45*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env2,(Env2-Zopt_Env2)/(67.45*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env3,(Env3-Zopt_Env3)/(67.45*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env4,(Env4-Zopt_Env4)/(67.45*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env5,(Env5-Zopt_Env5)/(67.45*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env6,(Env6-Zopt_Env6)/(67.45*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env7,(Env7-Zopt_Env7)/(67.45*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env8,(Env8-Zopt_Env8)/(67.45*2)),alpha=0.8) +
        scale_x_continuous(breaks=c(eval(parse(text=paste("Z_Model3_Cas", c,"_",q,"QTLs",sep=""))),eval(parse(text=paste("Z_Model3nul_Cas", c,"_",q,"QTLs",sep="")))[c(1,3,5,7)]),labels = c("E 1a", "E 1b","E 2a", "E 2b","E 3a",  "E 3b","E 4a", "E 4b", "E 1", "E 2", "E 3", "E 4")) +
        theme_bw() + theme(panel.grid=element_blank()) +
        theme(axis.text.x  = element_text(angle=-90, vjust= 1,hjust=0))
      
      #       for (i in (1:8)){
      #         h <-h + geom_vline(xintercept= eval(parse(text=paste("Z_Model3_Cas", c,"_",q,"QTLs",sep="")))[i],size = 0.5, linetype = "solid",colour = "purple")
      #         h <-h + geom_vline(xintercept= eval(parse(text=paste("Z_Model3nul_Cas", c,"_",q,"QTLs",sep="")))[i],size = 0.5, linetype = "dashed", colour = "orange")
      #       } 
      h <-h + geom_hline(yintercept= 0,size = 0.5, linetype = "dashed", colour = "gray")  
      h <- h+ geom_point(aes(x=c(-67.45,67.45), y=0.3), colour="gray", fill = "gray", text="Optimum",  size = 5, shape = 25,alpĥa = 0.4)
      if (c == "A"|c == "B") {h <- h + ylim(-0.3,0.3)} else {h <- h + ylim (-0.3,0.3)}      
      pdf(file = paste("ContinuousModel_Divided_PhenotValues_Case",c,"_",q,"QTLs",sep=""), height = 12, width = 16, paper = "special")        
      print(h)
      dev.off()
      
    }
  }
}

###### Autres possibilités de graphiques plus résumés #######
for (s in c("0","0.3")) {
  g <- ggplot(results[which(results$Selfing == s & results$Fecundity == "100"),])
  h <- g + facet_grid(Fecundity+QTLs~Cas, labeller = labeller(Fecundity=Fecun, Cas = Case,QTLs=label_both)) + aes(Selection,LagEnv,shape = Model)+ labs(y="Phenotypic lag for environments", title = paste("Continuous Mutation Model, Selfing Rate of",s)) + theme_bw(base_size=15)
  h <- h +  scale_shape_manual(values = c(15,0,16,1)) +
    guides(shape=guide_legend(direction = "horizontal",title.position = "top", title.hjust=0.5, title.theme = element_text(size = 15,angle = 0),label.position = "bottom", label.hjust=0.5,label.theme = element_text(angle = 90))) + scale_y_continuous(breaks=c(seq(0,24,3))) + theme(panel.grid.minor=element_blank()) +
    stat_summary(geom="pointrange", fun.data = mean_cl_boot, position = position_dodge(width = 0.75), size = 0.8) +
    scale_y_continuous(breaks=c(seq(0,0.2,0.05))) + coord_cartesian(ylim = c(-0.05, 0.2)) +
    theme(panel.grid.minor=element_blank())
  print(h)
  #Ecriture du pdf
  pdf(file = paste("ContinuousModel_Selfing",s,"_conti_Lag_Env_fec100.pdf",sep=""), height = 12, width = 16, paper = "special")
  print(h)
  dev.off() 
  
}
for (s in c("0","0.3")) {
  g <- ggplot(results[which(results$Selfing == s & results$Fecundity != "100"),])
  h <- g + facet_grid(Fecundity+QTLs~Cas, labeller = labeller(Fecundity=Fecun, Cas = Case,QTLs=label_both)) + aes(Selection,LagEnv,shape = Model)+ labs(y="Phenotypic lag for environments", title = paste("Continuous Mutation Model, Selfing Rate of",s)) + theme_bw(base_size=15)
  h <- h +  scale_shape_manual(values = c(15,0,16,1)) +
    guides(shape=guide_legend(direction = "horizontal",title.position = "top", title.hjust=0.5, title.theme = element_text(size = 15,angle = 0),label.position = "bottom", label.hjust=0.5,label.theme = element_text(angle = 90))) + scale_y_continuous(breaks=c(seq(0,24,3))) + theme(panel.grid.minor=element_blank()) +
    stat_summary(geom="pointrange", fun.data = mean_cl_boot, position = position_dodge(width = 0.75), size = 0.8) +
    scale_y_continuous(breaks=c(seq(0,0.2,0.05))) + coord_cartesian(ylim = c(-0.05, 0.2)) +
    theme(panel.grid.minor=element_blank())
  print(h)
  #Ecriture du pdf
  pdf(file = paste("ContinuousModel_Selfing",s,"_conti_Lag_Env_fec10-200.pdf",sep=""), height = 12, width = 16, paper = "special")
  print(h)
  dev.off() 
  
}


for (s in c("0","0.3")) {
  g <- ggplot(results[which(results$Selfing == s & results$Fecundity == "100"),])
  h <- g + facet_grid(Fecundity+QTLs~Cas, labeller = labeller(Fecundity=Fecun, Cas = Case,QTLs=label_both)) + aes(Selection,LagPop,shape = Model)+ labs(y="Phenotypic lag for populations", title = paste("Continuous Mutation Model, Selfing Rate of",s)) + theme_bw(base_size=15)
  h <- h +  scale_shape_manual(values = c(15,0,16,1)) +
    guides(shape=guide_legend(direction = "horizontal",title.position = "top", title.hjust=0.5, title.theme = element_text(size = 15,angle = 0),label.position = "bottom", label.hjust=0.5,label.theme = element_text(angle = 90))) + scale_y_continuous(breaks=c(seq(0,24,3))) + theme(panel.grid.minor=element_blank()) +
    stat_summary(geom="pointrange", fun.data = mean_cl_boot, position = position_dodge(width = 0.75), size = 0.8) +
    scale_y_continuous(breaks=c(seq(0,0.2,0.05)))  + coord_cartesian(ylim = c(-0.05, 0.2)) +
    theme(panel.grid.minor=element_blank())
  print(h)
  #Ecriture du pdf
  pdf(file = paste("ContinuousModel_Selfing",s,"_conti_Lag_Pop_fec100.pdf",sep=""), height = 12, width = 16, paper = "special")
  print(h)
  dev.off() 
  
}
for (s in c("0","0.3")) {
  g <- ggplot(results[which(results$Selfing == s & results$Fecundity != "100"),])
  h <- g + facet_grid(Fecundity+QTLs~Cas, labeller = labeller(Fecundity=Fecun, Cas = Case,QTLs=label_both)) + aes(Selection,LagPop,shape = Model)+ labs(y="Phenotypic lag for populations", title = paste("Continuous Mutation Model, Selfing Rate of",s)) + theme_bw(base_size=15)
  h <- h +  scale_shape_manual(values = c(15,0,16,1)) +
    guides(shape=guide_legend(direction = "horizontal",title.position = "top", title.hjust=0.5, title.theme = element_text(size = 15,angle = 0),label.position = "bottom", label.hjust=0.5,label.theme = element_text(angle = 90))) + scale_y_continuous(breaks=c(seq(0,24,3))) + theme(panel.grid.minor=element_blank()) +
    stat_summary(geom="pointrange", fun.data = mean_cl_boot, position = position_dodge(width = 0.75), size = 0.8) +
    scale_y_continuous(breaks=c(seq(0,0.2,0.05)))  + coord_cartesian(ylim = c(-0.05, 0.2)) +
    theme(panel.grid.minor=element_blank())
  print(h)
  #Ecriture du pdf
  pdf(file = paste("ContinuousModel_Selfing",s,"_conti_Lag_Pop_fec10-200.pdf",sep=""), height = 12, width = 16, paper = "special")
  print(h)
  dev.off() 
  
}


##### Sauvegarde du fichier #####

write.table(results, file= "Number_of_Occupied_Patches_Conti.txt")

##### Lecture du fichier et test des effets ####

results_conti<- read.table(file= "Number_of_Occupied_Patches_Conti.txt")
results_conti$freq_occupied_patches <- results_conti$nb_occupied_patches/24 ; results_conti$nb_unoccupied_patches <- 24 - results_conti$nb_occupied_patches
results_conti$Fecundity <- as.factor(results_conti$Fecundity);results_conti$QTLs <- as.factor(results_conti$QTLs);results_conti$Selfing <- as.factor(results_conti$Selfing);results_conti$Selection <- as.factor(results_conti$Selection)


######Et on fait des GLM######

#On change le système de contraintes par défaut
options(contrasts=c("contr.sum","contr.sum"))

#On fait le GLM sur l'ensemble des données
occupied_glm_conti<-glm(data = results_conti,formula = (nb_occupied_patches/24) ~ (Model+QTLs+Selfing+Cas+Fecundity+Selection)^2, weights = rep(24,times=length(nb_occupied_patches)), family = binomial(link = "logit"),start=NULL)
summary(occupied_glm_conti)
aov_glm_occupied_patches_conti <- anova(occupied_glm_conti,test="Chisq")
drop_glm_occupied_patches_conti <- drop1(occupied_glm_conti, scope = occupied_glm_conti$call, test = "Chisq")
coeff_glm_occupied_patches_conti <- dummy.coef(occupied_glm_conti)

write.table(aov_glm_occupied_patches_conti,file = "aov_glm_occupied_patches_conti.txt")
write.table(drop_glm_occupied_patches_conti,file = "drop_glm_occupied_patches_conti.txt")

for (i in 1:length(coeff_glm_occupied_patches_conti)) {
  write.table(coeff_glm_occupied_patches_conti[i], file = paste("glm_coeff_conti_",names(coeff_glm_occupied_patches_conti[i]),".txt",sep=""))
}

#On fait un GLM sur les modèles "nuls"
occupied_glm_conti_nul<-glm(data = results_conti[which(results_conti$Model=="Overlapping without within population environmental heterogeneity"|results_conti$Model == "Divided without within population environmental heterogeneity"),],
                        formula = (nb_occupied_patches/24) ~ (Model+QTLs+Selfing+Cas+Fecundity+Selection)^2,
                        weights = rep(24,times=length(nb_occupied_patches)),
                        family = binomial(link = "logit"),start=NULL)
summary(occupied_glm_conti_nul)
aov_glm_occupied_patches_conti_nul <- anova(occupied_glm_conti_nul,test="Chisq")
drop_glm_occupied_patches_conti_nul <- drop1(occupied_glm_conti_nul, scope = occupied_glm_conti_nul$call, test = "Chisq")
coeff_glm_occupied_patches_conti_nul <- dummy.coef(occupied_glm_conti_nul)

write.table(aov_glm_occupied_patches_conti_nul,file = "aov_glm_occupied_patches_conti_nul.txt")
write.table(drop_glm_occupied_patches_conti_nul,file = "drop_glm_occupied_patches_conti_nul.txt")

for (i in 1:length(coeff_glm_occupied_patches_conti_nul)) {
  write.table(coeff_glm_occupied_patches_conti_nul[i], file = paste("glm_coeff_conti_nul_",names(coeff_glm_occupied_patches_conti_nul[i]),".txt",sep=""))
}

#On fait un GLM sur les modèles avec hétérogénéité
occupied_glm_conti_het<-glm(results_conti[which(results_conti$Model=="Overlapping with within population environmental heterogeneity"|results_conti$Model == "Divided with within population environmental heterogeneity"),],
                                   formula = (nb_occupied_patches/24) ~ (Model+QTLs+Selfing+Cas+Fecundity+Selection)^2,
                                   weights = rep(24,times=length(nb_occupied_patches)),
                                   family = binomial(link = "logit"),start=NULL)
summary(occupied_glm_conti_het)
aov_glm_occupied_patches_conti_het <- anova(occupied_glm_conti_het,test="Chisq")
drop_glm_occupied_patches_conti_het <- drop1(occupied_glm_conti_het, scope = occupied_glm_conti_het$call, test = "Chisq")
coeff_glm_occupied_patches_conti_het <- dummy.coef(occupied_glm_conti_het)

write.table(aov_glm_occupied_patches_conti_het,file = "aov_glm_occupied_patches_conti_het.txt")
write.table(drop_glm_occupied_patches_conti_het,file = "drop_glm_occupied_patches_conti_het.txt")

for (i in 1:length(coeff_glm_occupied_patches_conti_het)) {
  write.table(coeff_glm_occupied_patches_conti_het[i], file = paste("glm_coeff_conti_het_",names(coeff_glm_occupied_patches_conti_het[i]),".txt",sep=""))
}

#On fait un GLM sur les modèles avec divided
occupied_glm_conti_div<-glm(results_conti[which(results_conti$Model== "Divided without within population environmental heterogeneity"|results_conti$Model == "Divided with within population environmental heterogeneity"),],
                            formula = (nb_occupied_patches/24) ~ (Model+QTLs+Selfing+Cas+Fecundity+Selection)^2,
                            weights = rep(24,times=length(nb_occupied_patches)),
                            family = binomial(link = "logit"),start=NULL)
summary(occupied_glm_conti_div)
aov_glm_occupied_patches_conti_div <- anova(occupied_glm_conti_div,test="Chisq")
drop_glm_occupied_patches_conti_div <- drop1(occupied_glm_conti_div, scope = occupied_glm_conti_div$call, test = "Chisq")
coeff_glm_occupied_patches_conti_div <- dummy.coef(occupied_glm_conti_div)

write.table(aov_glm_occupied_patches_conti_div,file = "aov_glm_occupied_patches_conti_div.txt")
write.table(drop_glm_occupied_patches_conti_div,file = "drop_glm_occupied_patches_conti_div.txt")

for (i in 1:length(coeff_glm_occupied_patches_conti_div)) {
  write.table(coeff_glm_occupied_patches_conti_div[i], file = paste("glm_coeff_conti_div_",names(coeff_glm_occupied_patches_conti_div[i]),".txt",sep=""))
}

#On fait un GLM sur les modèles avec overlapping
occupied_glm_conti_over<-glm(results_conti[which(results_conti$Model== "Overlapping without within population environmental heterogeneity"|results_conti$Model == "Overlapping with within population environmental heterogeneity"),],
                            formula = (nb_occupied_patches/24) ~ (Model+QTLs+Selfing+Cas+Fecundity+Selection)^2,
                            weights = rep(24,times=length(nb_occupied_patches)),
                            family = binomial(link = "logit"),start=NULL)
summary(occupied_glm_conti_over)
aov_glm_occupied_patches_conti_over <- anova(occupied_glm_conti_over,test="Chisq")
drop_glm_occupied_patches_conti_over <- drop1(occupied_glm_conti_over, scope = occupied_glm_conti_over$call, test = "Chisq")
coeff_glm_occupied_patches_conti_over <- dummy.coef(occupied_glm_conti_over)

write.table(aov_glm_occupied_patches_conti_over,file = "aov_glm_occupied_patches_conti_over.txt")
write.table(drop_glm_occupied_patches_conti_over,file = "drop_glm_occupied_patches_conti_over.txt")

for (i in 1:length(coeff_glm_occupied_patches_conti_over)) {
  write.table(coeff_glm_occupied_patches_conti_over[i], file = paste("glm_coeff_conti_over_",names(coeff_glm_occupied_patches_conti_over[i]),".txt",sep=""))
}

interaction.plot(results_conti$Model,results_conti$Cas,results_conti$nb_occupied_patches)
interaction.plot(results_conti$Cas,results_conti$QTLs,results_conti$nb_occupied_patches)
interaction.plot(results_conti$Selection,results_conti$QTLs,results_conti$nb_occupied_patches)
interaction.plot(results_conti$Cas,results_conti$Fecundity,results_conti$nb_occupied_patches)
interaction.plot(results_conti$Selection,results_conti$Cas,results_conti$nb_occupied_patches)
