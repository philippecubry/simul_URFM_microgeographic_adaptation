setwd("~/nemo_cluster")
FstTable<-NULL
FstTableConti<-NULL
####Definition des variables####
model <- c(2,3) ; cas <- c(3,4) ; fecundity <- c("010","100","200") ; selfing <- c(0,"03") ; selection <- c("01","05",20,50) ; qtls <- c(10,50)
repetition <- c("01","02","03","04","05","06","07","08","09","10")

####Lecture des modèles "normaux" bialléliques####
for (m in model) {
  for (c in cas) {
    for (f in fecundity) {
      for (s in selfing) {
        for (w in selection) {
          for (q in qtls) {
            for (r in repetition){
              Model <- m ; Selfing <- s ; QTLs <- q ; Fecundity <- f ; Selection <- w ; Cas <- c ; Repetition <- r
              if (file.exists(paste("Sims2/VarCompForFSTqs/var_comp_model",m,"_cas",c,"_",f,"fec_",s,"self_selection",w,"_",q,"QTLs_rep",r,".rda",sep="")))
              {load(paste("Sims2/VarCompForFSTqs/var_comp_model",m,"_cas",c,"_",f,"fec_",s,"self_selection",w,"_",q,"QTLs_rep",r,".rda",sep=""))
              FstPop <- diag(var_comp$F)[1] ; FstEnv <- diag(var_comp$F)[2] ; FstPatch <- diag(var_comp$F)[3] ; Fis <- diag(var_comp$F)[4]
              } else
              {
                FstPop <- NA ; FstEnv <- NA ; FstPatch <- NA ; Fis <- NA  
              }
            res_temp<-as.data.frame(cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,Repetition,FstPop,FstEnv,FstPatch,Fis))
            FstTable<-rbind(FstTable,res_temp)
            rm("var_comp","FstPop","FstEnv","FstPatch","Fis","Model","Cas","Selfing","QTLs","Fecundity","Selection","Repetition","res_temp")
            }
          }
        }
      }
    }
  }
}

####Lecture des modèles "nuls" bialléliques####
for (m in model) {
  for (c in cas) {
    for (f in fecundity) {
      for (s in selfing) {
        for (w in selection) {
          for (q in qtls) {
            for (r in repetition){
              Model <- paste(m," nul",sep="") ; Selfing <- s ; QTLs <- q ; Fecundity <- f ; Selection <- w ; Cas <- c ; Repetition <- r
              if (file.exists(paste("Sims2/VarCompForFSTqsNul/var_comp_model",m,"nul_cas",c,"_",f,"fec_",s,"self_selection",w,"_",q,"QTLs_rep",r,".rda",sep="")))
              {load(paste("Sims2/VarCompForFSTqsNul/var_comp_model",m,"nul_cas",c,"_",f,"fec_",s,"self_selection",w,"_",q,"QTLs_rep",r,".rda",sep=""))
               FstPop <- diag(var_comp$F)[1] ; FstEnv <- NA ; FstPatch <- diag(var_comp$F)[2] ; Fis <- diag(var_comp$F)[3]
              } else
              {
                FstPop <- NA ; FstEnv <- NA ; FstPatch <- NA ; Fis <- NA  
              }
              res_temp<-as.data.frame(cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,Repetition,FstPop,FstEnv,FstPatch,Fis))
              FstTable<-rbind(FstTable,res_temp)
              rm("var_comp","FstPop","FstEnv","FstPatch","Fis","Model","Cas","Selfing","QTLs","Fecundity","Selection","Repetition","res_temp")
            }
          }
        }
      }
    }
  }
}

####Lecture des modèles bialléliques à l'équilibre####
for (s in selfing) {
  for (q in qtls) {
    Model <- "equilibrium" ; Selfing <- s ; QTLs <- q ; Fecundity <- "up to capacity" ; Selection <- "minimal" ; Cas <- NA ; Repetition <- NA
    if (file.exists(paste("Sims2/VarCompForFSTqsEq/var_comp_equilibrium_",s,"self_",q,"QTLs.rda",sep="")))
    {load(paste("Sims2/VarCompForFSTqsEq/var_comp_equilibrium_",s,"self_",q,"QTLs.rda",sep=""))
     FstPop <- diag(var_comp$F)[1] ; FstEnv <- NA ; FstPatch <- diag(var_comp$F)[2] ; Fis <- diag(var_comp$F)[3]
    } else
    {
      FstPop <- NA ; FstEnv <- NA ; FstPatch <- NA ; Fis <- NA  
    }
    res_temp<-as.data.frame(cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,Repetition,FstPop,FstEnv,FstPatch,Fis))
    FstTable<-rbind(FstTable,res_temp)
    rm("var_comp","FstPop","FstEnv","FstPatch","Fis","Model","Cas","Selfing","QTLs","Fecundity","Selection","Repetition","res_temp")
  }
}


####Lecture des modèles "normaux" continus####
for (m in model) {
  for (c in cas) {
    for (f in fecundity) {
      for (s in selfing) {
        for (w in selection) {
          for (q in qtls) {
            for (r in repetition){
              Model <- paste(m,"_conti",sep="") ; Selfing <- s ; QTLs <- q ; Fecundity <- f ; Selection <- w ; Cas <- c ; Repetition <- r
              if (file.exists(paste("Sims2Conti/VarCompForFSTqsConti/var_comp_conti_model",m,"_cas",c,"_",f,"fec_",s,"self_selection",w,"_",q,"QTLs_rep",r,".rda",sep="")))
              {load(paste("Sims2Conti/VarCompForFSTqsConti/var_comp_conti_model",m,"_cas",c,"_",f,"fec_",s,"self_selection",w,"_",q,"QTLs_rep",r,".rda",sep=""))
               FstPop <- diag(var_comp$F)[1] ; FstEnv <- diag(var_comp$F)[2] ; FstPatch <- diag(var_comp$F)[3] ; Fis <- diag(var_comp$F)[4]
              } else
              {
                FstPop <- NA ; FstEnv <- NA ; FstPatch <- NA ; Fis <- NA  
              }
              res_temp<-as.data.frame(cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,Repetition,FstPop,FstEnv,FstPatch,Fis))
              FstTableConti<-rbind(FstTableConti,res_temp)
              rm("var_comp","FstPop","FstEnv","FstPatch","Fis","Model","Cas","Selfing","QTLs","Fecundity","Selection","Repetition","res_temp")
            }
          }
        }
      }
    }
  }
}

####Lecture des modèles "nuls" continus####
for (m in model) {
  for (c in cas) {
    for (f in fecundity) {
      for (s in selfing) {
        for (w in selection) {
          for (q in qtls) {
            for (r in repetition){
              Model <- paste(m," conti_nul",sep="") ; Selfing <- s ; QTLs <- q ; Fecundity <- f ; Selection <- w ; Cas <- c ; Repetition <- r
              if (file.exists(paste("Sims2Conti/VarCompForFSTqsNulConti/var_comp_conti_model",m,"nul_cas",c,"_",f,"fec_",s,"self_selection",w,"_",q,"QTLs_rep",r,".rda",sep="")))
              {load(paste("Sims2Conti/VarCompForFSTqsNulConti/var_comp_conti_model",m,"nul_cas",c,"_",f,"fec_",s,"self_selection",w,"_",q,"QTLs_rep",r,".rda",sep=""))
               FstPop <- diag(var_comp$F)[1] ; FstEnv <- NA ; FstPatch <- diag(var_comp$F)[2] ; Fis <- diag(var_comp$F)[3]
              } else
              {
                FstPop <- NA ; FstEnv <- NA ; FstPatch <- NA ; Fis <- NA  
              }
              res_temp<-as.data.frame(cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,Repetition,FstPop,FstEnv,FstPatch,Fis))
              FstTableConti<-rbind(FstTableConti,res_temp)
              rm("var_comp","FstPop","FstEnv","FstPatch","Fis","Model","Cas","Selfing","QTLs","Fecundity","Selection","Repetition","res_temp")
            }
          }
        }
      }
    }
  }
}

####Lecture des modèles continus à l'équilibre####
for (s in selfing) {
  for (q in qtls) {
    Model <- "conti_equilibrium" ; Selfing <- s ; QTLs <- q ; Fecundity <- "up to capacity" ; Selection <- "minimal" ; Cas <- NA ; Repetition <- NA
    if (file.exists(paste("Sims2Conti/VarCompForFSTqsEq/var_comp_equilibrium_",s,"self_",q,"QTLs.rda",sep="")))
    {load(paste("Sims2Conti/VarCompForFSTqsEq/var_comp_equilibrium_",s,"self_",q,"QTLs.rda",sep=""))
     FstPop <- diag(var_comp$F)[1] ; FstEnv <- NA ; FstPatch <- diag(var_comp$F)[2] ; Fis <- diag(var_comp$F)[3]
    } else
    {
      FstPop <- NA ; FstEnv <- NA ; FstPatch <- NA ; Fis <- NA  
    }
    res_temp<-as.data.frame(cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,Repetition,FstPop,FstEnv,FstPatch,Fis))
    FstTableConti<-rbind(FstTableConti,res_temp)
    rm("var_comp","FstPop","FstEnv","FstPatch","Fis","Model","Cas","Selfing","QTLs","Fecundity","Selection","Repetition","res_temp")
  }
}

####Sauvegarde des fichiers####
write.table(FstTable,file="FSTq_Sims2.txt")
write.table(FstTableConti,file="FSTq_Sims2Conti.txt")