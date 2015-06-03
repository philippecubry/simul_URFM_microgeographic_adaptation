#définit le répertoire de travail
setwd("~/nemo_cluster/Sims2")

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

##### Création d'un objet vide qui va accueillir les résultats #####
results = NULL

##### Boucle de lecture des fichiers #####
for(m in model) {
  for(s in selfing) {
    for(q in qtls) {
      for(f in fecundity) {
        for(w in selection) {
          for(c in cas) {
            temp <- read.delim(file=(paste("Sim2_Model",m,"_",s,"self_",q,"QTLs_fitperpop","/stats/Sim2_Model",m,"_",s,"self_",q,"QTLs_Cas",c,"_SelInt",w,"_fec",f,".txt",sep = "")),header=T)
            temp<-temp[which(temp$generation==200),]
            #Définition des noms de paramètres
            Model <- rep(m,length(temp[,1]))
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
                  temp_r[1,1:length(temp_r)]<-NA ; temp_r$Model <- m; temp_r$Selfing <- s ; temp_r$QTLs <- q; temp_r$Fecundity <- f ; temp_r$Selection <- w ; temp_r$Cas <- c; temp_r$replicate <- r
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

results_sauve <-results

##### On réordonne les cas 1, 2, 3, 4 en 2, 1, 3, 4 et on les renomme, on réordonne les sélections en 50,20,5,1 et les QTLs en 10,50 #####
results$Cas<-factor(results$Cas,levels=c("2","1","3","4")) ; levels(results$Cas) <- c("A","B","C","D")
results$Selection<-factor(results$Selection,levels=c("50","20","05","01"))
results$QTLs<-factor(results$QTLs,levels=c("10","50"))
results$Selfing<-factor(results$Selfing,levels=c("0","03")) ; levels(results$Selfing) <- c("0","0.3")
results$Model<-factor(results$Model,levels=c("2","2nul","3","3nul")) ; levels(results$Model) <- c("Overlapping with within population environmental heterogeneity","Overlapping without within population environmental heterogeneity","Divided with within population environmental heterogeneity","Divided without within population environmental heterogeneity")

Case <- c('A'='Case A','B'='Case B','C'='Case C','D'='Case D'); Fecun <- c("010"="Mean Fecundity = 10","100"="Mean Fecundity = 100","200"="Mean Fecundity = 200")

##### Création des graphes d'occupation des patches #####
for (s in c("0","0.3")) {
  g <- ggplot(results[which(results$Selfing == s & results$Fecundity == "100"),])
  h <- g + facet_grid(Fecundity+QTLs~Cas, labeller = labeller(Fecundity=Fecun, Cas = Case,QTLs=label_both)) + aes(Selection,nb_occupied_patches,shape = Model)+ labs(y="Number of Occupied Patches", title = paste("Diallelic Mutation Model, Selfing Rate of",s)) + theme_bw(base_size=15)
  h <- h +  scale_shape_manual(values = c(15,0,16,1)) +
    guides(shape=guide_legend(direction = "horizontal",title.position = "top", title.hjust=0.5, title.theme = element_text(size = 15,angle = 0),label.position = "bottom", label.hjust=0.5,label.theme = element_text(angle = 90))) + scale_y_continuous(breaks=c(seq(0,24,3))) + theme(panel.grid.minor=element_blank()) +
    stat_summary(geom="pointrange", fun.y = median, fun.ymin = min, fun.ymax=max, position = position_dodge(width = 0.75), size = 0.8) +
    scale_y_continuous(breaks=c(seq(0,24,3))) +
    theme(panel.grid.minor=element_blank())
#Ecriture du pdf
  pdf(file = paste("DiallelicModel_Selfing",s,"_nb_occupied_patches_fec100.pdf",sep=""), height = 12, width = 16, paper = "special")
  print(h)
  dev.off()  
  
}
for (s in c("0","0.3")) {
  g <- ggplot(results[which(results$Selfing == s & results$Fecundity != "100"),])
  h <- g + facet_grid(Fecundity+QTLs~Cas, labeller = labeller(Fecundity=Fecun, Cas = Case,QTLs=label_both)) + aes(Selection,nb_occupied_patches,shape = Model)+ labs(y="Number of Occupied Patches", title = paste("Diallelic Mutation Model, Selfing Rate of",s)) + theme_bw(base_size=15)
  h <- h +  scale_shape_manual(values = c(15,0,16,1)) +
    guides(shape=guide_legend(direction = "horizontal",title.position = "top", title.hjust=0.5, title.theme = element_text(size = 15,angle = 0),label.position = "bottom", label.hjust=0.5,label.theme = element_text(angle = 90))) + scale_y_continuous(breaks=c(seq(0,24,3))) + theme(panel.grid.minor=element_blank()) +
    stat_summary(geom="pointrange", fun.y = median, fun.ymin = min, fun.ymax=max, position = position_dodge(width = 0.75), size = 0.8) +
    scale_y_continuous(breaks=c(seq(0,24,3))) +
    theme(panel.grid.minor=element_blank())
  #Ecriture du pdf
  pdf(file = paste("DiallelicModel_Selfing",s,"_nb_occupied_patches_fec10-200.pdf",sep=""), height = 12, width = 16, paper = "special")
  print(h)
  dev.off()  
  
}


##### Definition des optimums locaux demandés #####
results$Zopt_Env1 <- NA;results$Zopt_Env2 <- NA;results$Zopt_Env3 <- NA;results$Zopt_Env4 <- NA;results$Zopt_Env5 <- NA;results$Zopt_Env6 <- NA;results$Zopt_Env7 <- NA;results$Zopt_Env8 <- NA

Z_Model2_CasB_10QTLs <- c(-13.557282933,-4.519094311,-7.53182385166667,1.5063647703,-1.5063647703,7.53182385166667,4.519094311,13.557282933)
Z_Model2_CasA_10QTLs <- c(-16.9466036663,-5.6488678888,-9.4147798146,1.8829559629,-1.8829559629,9.4147798146,5.6488678888,16.9466036663)
Z_Model2_CasC_10QTLs <- c(-6.7786414665,-2.2595471555,-3.76591192583333,0.7531823852,-0.7531823852,3.76591192583333,2.2595471555,6.7786414665)
Z_Model2_CasD_10QTLs <- c(-3.38932073325,-1.12977357775,-1.88295596291667,0.3765911926,-0.3765911926,1.88295596291667,1.12977357775,3.38932073325)
Z_Model2_CasB_50QTLs <- c(-80.325348195,-26.775116065,-44.6251934417,8.92503868833333,-8.92503868833333,44.6251934417,26.775116065,80.325348195)
Z_Model2_CasA_50QTLs <- c(-100.40668524375,-33.46889508125,-55.7814918020833,11.1562983604167,-11.1562983604167,55.7814918020833,33.46889508125,100.40668524375)
Z_Model2_CasC_50QTLs <- c(-40.1626740975,-13.3875580325,-22.3125967208333,4.46251934416667,-4.46251934416667,22.3125967208333,13.3875580325,40.1626740975)
Z_Model2_CasD_50QTLs <- c(-20.08133704875,-6.69377901625,-11.1562983604167,2.23125967208333,-2.23125967208333,11.1562983604167,6.69377901625,20.08133704875)

Z_Model3_CasB_10QTLs <- c(-13.557282933,-9.68377352357143,-5.81026411414286,-1.93675470471429,1.93675470471429,5.81026411414286,9.68377352357143,13.557282933)
Z_Model3_CasA_10QTLs <- c(-16.94660366625,-12.1047169044643,-7.26283014267857,-2.42094338089286,2.42094338089286,7.26283014267857,12.1047169044643,16.94660366625)
Z_Model3_CasC_10QTLs <- c(-6.7786414665,-4.84188676178571,-2.90513205707143,-0.968377352357143,0.968377352357143,2.90513205707143,4.84188676178571,6.7786414665)
Z_Model3_CasD_10QTLs <- c(-3.38932073325,-2.42094338089286,-1.45256602853571,-0.484188676178571,0.484188676178571,1.45256602853571,2.42094338089286,3.38932073325)
Z_Model3_CasB_50QTLs <- c(-80.325348195,-57.3752487107143,-34.4251492264286,-11.4750497421429,11.4750497421429,34.4251492264286,57.3752487107143,80.325348195)
Z_Model3_CasA_50QTLs <- c(-100.40668524375,-71.7190608883929,-43.0314365330357,-14.3438121776786,14.3438121776786,43.0314365330357,71.7190608883929,100.40668524375)
Z_Model3_CasC_50QTLs <- c(-40.1626740975,-28.6876243553571,-17.2125746132143,-5.73752487107143,5.73752487107143,17.2125746132143,28.6876243553571,40.1626740975)
Z_Model3_CasD_50QTLs <- c(-20.08133704875,-14.3438121776786,-8.60628730660714,-2.86876243553571,2.86876243553571,8.60628730660714,14.3438121776786,20.08133704875)

Z_Model2nul_CasB_10QTLs <- c(-9.038188622  ,-9.038188622  ,  -3.01272954066667  ,  -3.01272954066667	,	3.01272954066667	,  3.01272954066667	,	9.038188622,  9.038188622)
Z_Model2nul_CasA_10QTLs <- c(-11.2977357775  ,-11.2977357775  ,  -3.76591192583333	,  -3.76591192583333	,	3.76591192583333	,  3.76591192583333	,	11.2977357775,  11.2977357775)
Z_Model2nul_CasC_10QTLs <- c(-4.519094311  ,-4.519094311  ,  -1.50636477033333	,  -1.50636477033333	,	1.50636477033333	,1.50636477033333  ,	4.519094311,  4.519094311)
Z_Model2nul_CasD_10QTLs <- c(-2.2595471555  ,-2.2595471555  ,  -0.753182385166667	,-0.753182385166667  ,	0.753182385166667	,  0.753182385166667	,	2.2595471555,  2.2595471555)
Z_Model2nul_CasB_50QTLs <- c(-53.55023213  ,-53.55023213  ,  -17.8500773766667  ,  -17.8500773766667	,	17.8500773766667	,	  17.8500773766667	,	53.55023213,  53.55023213)
Z_Model2nul_CasA_50QTLs <- c(-66.9377901625  ,-66.9377901625  ,  -22.3125967208333	,  -22.3125967208333	,	22.3125967208333	,22.3125967208333  ,	66.9377901625,  66.9377901625)
Z_Model2nul_CasC_50QTLs <- c(-26.775116065  ,-26.775116065  ,  -8.92503868833333	,  -8.92503868833333	,	8.92503868833333	,  8.92503868833333	,	26.775116065,  26.775116065)
Z_Model2nul_CasD_50QTLs <- c(-13.3875580325  ,-13.3875580325  ,  -4.46251934416667	,  -4.46251934416667	,	4.46251934416667	,  4.46251934416667	,	13.3875580325,  13.3875580325)

Z_Model3nul_CasB_10QTLs <- c(-11.6205282282857  ,-11.6205282282857  ,  -3.87350940942857  ,  -3.87350940942857	,	3.87350940942857	,  3.87350940942857	,	11.6205282282857,  11.6205282282857)
Z_Model3nul_CasA_10QTLs <- c(-14.5256602853571  ,-14.5256602853571  ,  -4.84188676178571	,  -4.84188676178571	,	4.84188676178571	,  4.84188676178571	,	14.5256602853571,  14.5256602853571)
Z_Model3nul_CasC_10QTLs <- c(-5.81026411414286  ,-5.81026411414286  ,  -1.93675470471429	,  -1.93675470471429	,	1.93675470471429	,  1.93675470471429	,	5.81026411414286,  5.81026411414286)
Z_Model3nul_CasD_10QTLs <- c(-2.90513205707143  ,-2.90513205707143  ,  -0.968377352357143	,  -0.968377352357143	,	0.968377352357143	,  0.968377352357143	,	2.90513205707143,  2.90513205707143)
Z_Model3nul_CasB_50QTLs <- c(-68.8502984528571  ,-68.8502984528571  ,  -22.9500994842857  ,  -22.9500994842857	,	22.9500994842857	,  22.9500994842857	,	68.8502984528571,  68.8502984528571)
Z_Model3nul_CasA_50QTLs <- c(-86.0628730660714  ,-86.0628730660714  ,  -28.6876243553571	,  -28.6876243553571	,	28.6876243553571	,  28.6876243553571	,	86.0628730660714,  86.0628730660714)
Z_Model3nul_CasC_50QTLs <- c(-34.4251492264286  ,-34.4251492264286  ,  -11.4750497421429	,  -11.4750497421429	,	11.4750497421429	,  11.4750497421429	,	34.4251492264286,  34.4251492264286)
Z_Model3nul_CasD_50QTLs <- c(-17.2125746132143  ,-17.2125746132143  ,  -5.73752487107143	,  -5.73752487107143	,	5.73752487107143	,  5.73752487107143	,	17.2125746132143,  17.2125746132143)

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
        labs(x= "Asked phenotypic value",y="Difference between realized and asked phenotypic value normalized by theoretical max genotypic range", title = paste("Diallelic Mutation Model, Overlapping Model, Case ",c, ", ", q, "QTLs"))  +
        geom_point(aes(Zopt_Env1,(Env1-Zopt_Env1)/(13.557282933*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env2,(Env2-Zopt_Env2)/(13.557282933*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env3,(Env3-Zopt_Env3)/(13.557282933*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env4,(Env4-Zopt_Env4)/(13.557282933*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env5,(Env5-Zopt_Env5)/(13.557282933*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env6,(Env6-Zopt_Env6)/(13.557282933*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env7,(Env7-Zopt_Env7)/(13.557282933*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env8,(Env8-Zopt_Env8)/(13.557282933*2)),alpha=0.8) +
        scale_x_continuous(breaks=c(eval(parse(text=paste("Z_Model2_Cas", c,"_",q,"QTLs",sep=""))),eval(parse(text=paste("Z_Model2nul_Cas", c,"_",q,"QTLs",sep="")))[c(1,3,5,7)]),labels = c("E 1a", "E 1b","E 2a", "E 2b","E 3a",  "E 3b","E 4a", "E 4b", "E 1", "E 2", "E 3", "E 4")) +
        theme_bw() + theme(panel.grid=element_blank()) +
        theme(axis.text.x  = element_text(angle=-90, vjust=1, hjust=0, size = 8))
      
#       for (i in (1:8)){
#         h <-h + geom_vline(xintercept= eval(parse(text=paste("Z_Model2_Cas", c,"_",q,"QTLs",sep="")))[i],size = 0.5, linetype = "solid",colour = "purple")
#         h <-h + geom_vline(xintercept= eval(parse(text=paste("Z_Model2nul_Cas", c,"_",q,"QTLs",sep="")))[i],size = 0.5, linetype = "dashed", colour = "orange")
#       }
      h <-h + geom_hline(yintercept= 0,size = 0.5, linetype = "dashed", colour = "gray")  
      h <- h+ geom_point(aes(x=c(-13.557282933,13.557282933), y=0.3), colour="gray", fill = "gray", text="Optimum",  size = 5, shape = 25,alpha = 0.6)
      if (c == "A"|c == "B") {h <- h + ylim(-0.3,0.3)} else {h <- h + ylim (-0.3,0.3)}
      pdf(file = paste("DiallelicModel_Overlapping_PhenotValues_Case",c,"_",q,"QTLs",sep=""), height = 12, width = 16, paper = "special")        
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
        labs(x= "Asked phenotypic value",y="Difference between realized and asked phenotypic value normalized by theoretical max genotypic range", title = paste("Diallelic Mutation Model, Overlapping Model, Case ",c, ", ", q, "QTLs"))  +
        geom_point(aes(Zopt_Env1,(Env1-Zopt_Env1)/(80.325348195*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env2,(Env2-Zopt_Env2)/(80.325348195*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env3,(Env3-Zopt_Env3)/(80.325348195*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env4,(Env4-Zopt_Env4)/(80.325348195*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env5,(Env5-Zopt_Env5)/(80.325348195*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env6,(Env6-Zopt_Env6)/(80.325348195*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env7,(Env7-Zopt_Env7)/(80.325348195*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env8,(Env8-Zopt_Env8)/(80.325348195*2)),alpha=0.8) +
        scale_x_continuous(breaks=c(eval(parse(text=paste("Z_Model2_Cas", c,"_",q,"QTLs",sep=""))),eval(parse(text=paste("Z_Model2nul_Cas", c,"_",q,"QTLs",sep="")))[c(1,3,5,7)]),labels = c("E 1a", "E 1b","E 2a", "E 2b","E 3a",  "E 3b","E 4a", "E 4b", "E 1", "E 2", "E 3", "E 4")) +
        theme_bw() + theme(panel.grid=element_blank()) +
        theme(axis.text.x  = element_text(angle=-90, vjust= 1,hjust=0,size = 8))
      
#       for (i in (1:8)){
#         h <-h + geom_vline(xintercept= eval(parse(text=paste("Z_Model2_Cas", c,"_",q,"QTLs",sep="")))[i],size = 0.5, linetype = "solid",colour = "purple")
#         h <-h + geom_vline(xintercept= eval(parse(text=paste("Z_Model2nul_Cas", c,"_",q,"QTLs",sep="")))[i],size = 0.5, linetype = "dashed", colour = "orange")
#           } 
    h <-h + geom_hline(yintercept= 0,size = 0.5, linetype = "dashed", colour = "gray")  
    h <- h+ geom_point(aes(x=c(-80.325348195,80.325348195), y=0.3), colour="gray", fill = "gray", text="Optimum",  size = 5, shape = 25,alpha = 0.6)
    if (c == "A"|c == "B") {h <- h + ylim(-0.3,0.3)} else {h <- h + ylim (-0.3,0.3)}    
    pdf(file = paste("DiallelicModel_Overlapping_PhenotValues_Case",c,"_",q,"QTLs",sep=""), height = 12, width = 16, paper = "special")        
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
        labs(x= "Asked phenotypic value",y="Difference between realized and asked phenotypic value normalized by theoretical max genotypic range", title = paste("Diallelic Mutation Model, Divided Model, Case ",c, ", ", q, "QTLs"))  +
        geom_point(aes(Zopt_Env1,(Env1-Zopt_Env1)/(13.557282933*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env2,(Env2-Zopt_Env2)/(13.557282933*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env3,(Env3-Zopt_Env3)/(13.557282933*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env4,(Env4-Zopt_Env4)/(13.557282933*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env5,(Env5-Zopt_Env5)/(13.557282933*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env6,(Env6-Zopt_Env6)/(13.557282933*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env7,(Env7-Zopt_Env7)/(13.557282933*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env8,(Env8-Zopt_Env8)/(13.557282933*2)),alpha=0.8) +
        scale_x_continuous(breaks=c(eval(parse(text=paste("Z_Model3_Cas", c,"_",q,"QTLs",sep=""))),eval(parse(text=paste("Z_Model3nul_Cas", c,"_",q,"QTLs",sep="")))[c(1,3,5,7)]),labels = c("E 1a", "E 1b","E 2a", "E 2b","E 3a",  "E 3b","E 4a", "E 4b", "E 1", "E 2", "E 3", "E 4")) +
        theme_bw() + theme(panel.grid=element_blank()) + 
        theme(axis.text.x  = element_text(angle=-90, vjust=1, hjust=0))
      
#       for (i in (1:8)){
#         h <-h + geom_vline(xintercept= eval(parse(text=paste("Z_Model3_Cas", c,"_",q,"QTLs",sep="")))[i],size = 0.5, linetype = "solid",colour = "purple")
#         h <-h + geom_vline(xintercept= eval(parse(text=paste("Z_Model3nul_Cas", c,"_",q,"QTLs",sep="")))[i],size = 0.5, linetype = "dashed", colour = "orange")
#       }
      h <-h + geom_hline(yintercept= 0,size = 0.5, linetype = "dashed", colour = "gray")  
      h <- h+ geom_point(aes(x=c(-13.557282933,13.557282933), y=0.3), colour="gray", fill = "gray", text="Optimum",  size = 5, shape = 25,alpha = 0.6)
      if (c == "A"|c == "B") {h <- h + ylim(-0.3,0.3)} else {h <- h + ylim (-0.3,0.3)}
      pdf(file = paste("DiallelicModel_Divided_PhenotValues_Case",c,"_",q,"QTLs",sep=""), height = 12, width = 16, paper = "special")        
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
        labs(x= "Asked phenotypic value",y="Difference between realized and asked phenotypic value normalized by theoretical max genotypic range", title = paste("Diallelic Mutation Model, Divided Model, Case ",c, ", ", q, "QTLs"))  +
        geom_point(aes(Zopt_Env1,(Env1-Zopt_Env1)/(80.325348195*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env2,(Env2-Zopt_Env2)/(80.325348195*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env3,(Env3-Zopt_Env3)/(80.325348195*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env4,(Env4-Zopt_Env4)/(80.325348195*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env5,(Env5-Zopt_Env5)/(80.325348195*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env6,(Env6-Zopt_Env6)/(80.325348195*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env7,(Env7-Zopt_Env7)/(80.325348195*2)),alpha=0.8) +
        geom_point(aes(Zopt_Env8,(Env8-Zopt_Env8)/(80.325348195*2)),alpha=0.8) +
        scale_x_continuous(breaks=c(eval(parse(text=paste("Z_Model3_Cas", c,"_",q,"QTLs",sep=""))),eval(parse(text=paste("Z_Model3nul_Cas", c,"_",q,"QTLs",sep="")))[c(1,3,5,7)]),labels = c("E 1a", "E 1b","E 2a", "E 2b","E 3a",  "E 3b","E 4a", "E 4b", "E 1", "E 2", "E 3", "E 4")) +
        theme_bw() + theme(panel.grid=element_blank()) +
        theme(axis.text.x  = element_text(angle=-90, vjust= 1,hjust=0))
      
#       for (i in (1:8)){
#         h <-h + geom_vline(xintercept= eval(parse(text=paste("Z_Model3_Cas", c,"_",q,"QTLs",sep="")))[i],size = 0.5, linetype = "solid",colour = "purple")
#         h <-h + geom_vline(xintercept= eval(parse(text=paste("Z_Model3nul_Cas", c,"_",q,"QTLs",sep="")))[i],size = 0.5, linetype = "dashed", colour = "orange")
#       } 
      h <-h + geom_hline(yintercept= 0,size = 0.5, linetype = "dashed", colour = "gray")  
      h <- h+ geom_point(aes(x=c(-80.325348195,80.325348195), y=0.3), colour="gray", fill = "gray", text="Optimum",  size = 5, shape = 25,alpha = 0.6)
      if (c == "A"|c == "B") {h <- h + ylim(-0.3,0.3)} else {h <- h + ylim (-0.3,0.3)}      
      pdf(file = paste("DiallelicModel_Divided_PhenotValues_Case",c,"_",q,"QTLs",sep=""), height = 12, width = 16, paper = "special")        
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
  pdf(file = paste("DiallelicModel_Selfing",s,"_conti_Lag_Env_fec100.pdf",sep=""), height = 12, width = 16, paper = "special")
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
  pdf(file = paste("DiallelicModel_Selfing",s,"_conti_Lag_Env_fec10-200.pdf",sep=""), height = 12, width = 16, paper = "special")
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
  pdf(file = paste("DiallelicModel_Selfing",s,"_conti_Lag_Pop_fec100.pdf",sep=""), height = 12, width = 16, paper = "special")
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
  pdf(file = paste("DiallelicModel_Selfing",s,"_conti_Lag_Pop_fec10-200.pdf",sep=""), height = 12, width = 16, paper = "special")
  print(h)
  dev.off() 
  
}

##### Sauvegarde des fichiers #####  
write.table(results, file= "Number_of_Occupied_Patches.txt")

##### Lecture des fichiers #####
results<- read.table(file= "Number_of_Occupied_Patches.txt")
results$freq_occupied_patches <- results$nb_occupied_patches/24 ; results$nb_unoccupied_patches <- 24 - results$nb_occupied_patches
results$Fecundity <- as.factor(results$Fecundity);results$QTLs <- as.factor(results$QTLs);results$Selfing <- as.factor(results$Selfing);results$Selection <- as.factor(results$Selection)

######Et on fait des GLM######

#On change le système de contraintes par défaut
options(contrasts=c("contr.sum","contr.sum"))

#On fait le GLM sur l'ensemble des données
occupied_glm<-glm(data = results,formula = (nb_occupied_patches/24) ~ (Model+ Selfing+QTLs +Fecundity+ Model/Cas +Selection), weights = rep(24,times=length(nb_occupied_patches)), family = binomial(link = "logit"),start=NULL)
summary(occupied_glm)
aov_glm_occupied_patches <- anova(occupied_glm,test="Chisq")
drop_glm_occupied_patches <- drop1(occupied_glm, scope = occupied_glm$call, test = "Chisq")
coeff_glm_occupied_patches <- dummy.coef(occupied_glm)

write.table(aov_glm_occupied_patches,file = "aov_glm_occupied_patches.txt")
write.table(drop_glm_occupied_patches,file = "drop_glm_occupied_patches.txt")

for (i in 1:length(coeff_glm_occupied_patches_conti)) {
  write.table(coeff_glm_occupied_patches[i], file = paste("glm_coeff_",names(coeff_glm_occupied_patches[i]),".txt",sep=""))
}

#On fait un GLM sur les modèles "nuls"
occupied_glm_nul<-glm(data = results[which(results$Model=="Overlapping without within population environmental heterogeneity"|results$Model == "Divided without within population environmental heterogeneity"),],
                            formula = (nb_occupied_patches/24) ~ (Model+QTLs+Selfing+Cas+Fecundity+Selection)^2,
                            weights = rep(24,times=length(nb_occupied_patches)),
                            family = binomial(link = "logit"),start=NULL)
summary(occupied_glm_nul)
aov_glm_occupied_patches_nul <- anova(occupied_glm_nul,test="Chisq")
drop_glm_occupied_patches_nul <- drop1(occupied_glm_nul, scope = occupied_glm_nul$call, test = "Chisq")
coeff_glm_occupied_patches_nul <- dummy.coef(occupied_glm_nul)

write.table(aov_glm_occupied_patches_nul,file = "aov_glm_occupied_patches_nul.txt")
write.table(drop_glm_occupied_patches_nul,file = "drop_glm_occupied_patches_nul.txt")

for (i in 1:length(coeff_glm_occupied_patches_nul)) {
  write.table(coeff_glm_occupied_patches_nul[i], file = paste("glm_coeff_nul_",names(coeff_glm_occupied_patches_nul[i]),".txt",sep=""))
}

#On fait un GLM sur les modèles avec hétérogénéité
occupied_glm_het<-glm(results[which(results$Model=="Overlapping with within population environmental heterogeneity"|results$Model == "Divided with within population environmental heterogeneity"),],
                            formula = (nb_occupied_patches/24) ~ (Model+Model/Cas+Selfing+QTLs+Fecundity+Selection),
                            weights = rep(24,times=length(nb_occupied_patches)),
                            family = binomial(link = "logit"),start=NULL)
summary(occupied_glm_het)
aov_glm_occupied_patches_het <- anova(occupied_glm_het,test="Chisq")
drop_glm_occupied_patches_het <- drop1(occupied_glm_het, scope = occupied_glm_het$call, test = "Chisq")
coeff_glm_occupied_patches_het <- dummy.coef(occupied_glm_het)

write.table(aov_glm_occupied_patches_het,file = "aov_glm_occupied_patches_het.txt")
write.table(drop_glm_occupied_patches_het,file = "drop_glm_occupied_patches_het.txt")

for (i in 1:length(coeff_glm_occupied_patches_het)) {
  write.table(coeff_glm_occupied_patches_het[i], file = paste("glm_coeff_het_",names(coeff_glm_occupied_patches_het[i]),".txt",sep=""))
}

#On fait un GLM sur les modèles avec divided
occupied_glm_div<-glm(results[which(results$Model== "Divided without within population environmental heterogeneity"|results$Model == "Divided with within population environmental heterogeneity"),],
                            formula = (nb_occupied_patches/24) ~ (Model+QTLs+Fecundity+Model/Cas+Selection+Selfing),
                            weights = rep(24,times=length(nb_occupied_patches)),
                            family = binomial(link = "logit"),start=NULL)
summary(occupied_glm_div)
aov_glm_occupied_patches_div <- anova(occupied_glm_div,test="Chisq")
drop_glm_occupied_patches_div <- drop1(occupied_glm_div, scope= ~(Model + Cas + QTLs + Selfing + Fecundity + Selection)^2, test = "Chisq")
coeff_glm_occupied_patches_div <- dummy.coef(occupied_glm_div)

write.table(aov_glm_occupied_patches_div,file = "aov_glm_occupied_patches_div.txt")
write.table(drop_glm_occupied_patches_div,file = "drop_glm_occupied_patches_div.txt")

for (i in 1:length(coeff_glm_occupied_patches_div)) {
  write.table(coeff_glm_occupied_patches_div[i], file = paste("glm_coeff_div_",names(coeff_glm_occupied_patches_div[i]),".txt",sep=""))
}

#On fait un GLM sur les modèles avec overlapping
occupied_glm_over<-glm(results[which(results$Model== "Overlapping without within population environmental heterogeneity"|results$Model == "Overlapping with within population environmental heterogeneity"),],
                             formula = (nb_occupied_patches/24) ~ (Model+QTLs+Selfing+Cas+Fecundity+Selection)^2,
                             weights = rep(24,times=length(nb_occupied_patches)),
                             family = binomial(link = "logit"),start=NULL)
summary(occupied_glm_over)
aov_glm_occupied_patches_over <- anova(occupied_glm_over,test="Chisq")
drop_glm_occupied_patches_over <- drop1(occupied_glm_over, scope = occupied_glm_over$call, test = "Chisq")
coeff_glm_occupied_patches_over <- dummy.coef(occupied_glm_over)

write.table(aov_glm_occupied_patches_over,file = "aov_glm_occupied_patches_over.txt")
write.table(drop_glm_occupied_patches_over,file = "drop_glm_occupied_patches_over.txt")

for (i in 1:length(coeff_glm_occupied_patches_over)) {
  write.table(coeff_glm_occupied_patches_over[i], file = paste("glm_coeff_over_",names(coeff_glm_occupied_patches_over[i]),".txt",sep=""))
}

interaction.plot(results$Model,results$Cas,results$nb_occupied_patches)
interaction.plot(results$Cas,results$QTLs,results$nb_occupied_patches)
interaction.plot(results$Selection,results$QTLs,results$nb_occupied_patches)
interaction.plot(results$Cas,results$Fecundity,results$nb_occupied_patches)
interaction.plot(results$Selection,results$Cas,results$nb_occupied_patches)

g <- ggplot(results)
g + aes(colour = Selection) + facet_grid(QTLs~Model+Cas) +
  geom_boxplot(aes(1,(adlt.W.avg.p1+adlt.W.avg.p2+adlt.W.avg.p3)/3)) +
  geom_boxplot(aes(2,(adlt.W.avg.p4+adlt.W.avg.p5+adlt.W.avg.p6)/3)) +
  geom_boxplot(aes(3,(adlt.W.avg.p7+adlt.W.avg.p8+adlt.W.avg.p9)/3)) +
  geom_boxplot(aes(4,(adlt.W.avg.p10+adlt.W.avg.p11+adlt.W.avg.p12)/3)) +
  geom_boxplot(aes(5,(adlt.W.avg.p13+adlt.W.avg.p14+adlt.W.avg.p15)/3)) +
  geom_boxplot(aes(6,(adlt.W.avg.p16+adlt.W.avg.p17+adlt.W.avg.p18)/3)) +
  geom_boxplot(aes(7,(adlt.W.avg.p19+adlt.W.avg.p20+adlt.W.avg.p21)/3)) +
  geom_boxplot(aes(8,(adlt.W.avg.p22+adlt.W.avg.p23+adlt.W.avg.p24)/3))