a=commandArgs()[4]
setwd("/mnt/travail/1_Base_de_donnees/Vautour/3_hist_de_vie/4_modeles_etats_baguage_mort/2_all_birds/1_model_S/chat_236")
library(RMark)
s.models=function()
{
  vaut=read.table("hdv_2012b.csv",sep="\t",h=T,colClasses="character")
  vaut$fond_nl=as.factor(vaut$fond_nl)
  vaut$age_lach=as.factor(vaut$age)
  vaut.process.age=process.data(vaut,model="MSLiveDead",groups=c("fond_nl","age_lach"),age.var=2,initial.ages=c(0,1,2,3,4,5),begin.time=1980)#
  vaut.ddl=make.design.data(vaut.process.age)


psi1980.bG=as.numeric(row.names(vaut.ddl$Psi[vaut.ddl$Psi$time==c(1981)&vaut.ddl$Psi$stratum=="b"&vaut.ddl$Psi$tostratum=="G",]))
psi1981.bG=as.numeric(row.names(vaut.ddl$Psi[vaut.ddl$Psi$time==c(1982)&vaut.ddl$Psi$stratum=="b"&vaut.ddl$Psi$tostratum=="G",]))
psi1982.bG=as.numeric(row.names(vaut.ddl$Psi[vaut.ddl$Psi$time==c(1983)&vaut.ddl$Psi$stratum=="b"&vaut.ddl$Psi$tostratum=="G",]))
psi1983.bG=as.numeric(row.names(vaut.ddl$Psi[vaut.ddl$Psi$time==c(1984)&vaut.ddl$Psi$stratum=="b"&vaut.ddl$Psi$tostratum=="G",]))
psi1984.bG=as.numeric(row.names(vaut.ddl$Psi[vaut.ddl$Psi$time==c(1985)&vaut.ddl$Psi$stratum=="b"&vaut.ddl$Psi$tostratum=="G",]))
psi1985.bG=as.numeric(row.names(vaut.ddl$Psi[vaut.ddl$Psi$time==c(1986)&vaut.ddl$Psi$stratum=="b"&vaut.ddl$Psi$tostratum=="G",]))
psi1986.bG=as.numeric(row.names(vaut.ddl$Psi[vaut.ddl$Psi$time==c(1987)&vaut.ddl$Psi$stratum=="b"&vaut.ddl$Psi$tostratum=="G",]))
psi1989.bG=as.numeric(row.names(vaut.ddl$Psi[vaut.ddl$Psi$time==c(1990)&vaut.ddl$Psi$stratum=="b"&vaut.ddl$Psi$tostratum=="G",]))
psi1990.bG=as.numeric(row.names(vaut.ddl$Psi[vaut.ddl$Psi$time==c(1991)&vaut.ddl$Psi$stratum=="b"&vaut.ddl$Psi$tostratum=="G",]))
psi1991.bG=as.numeric(row.names(vaut.ddl$Psi[vaut.ddl$Psi$time==c(1992)&vaut.ddl$Psi$stratum=="b"&vaut.ddl$Psi$tostratum=="G",]))
psi1992.bG=as.numeric(row.names(vaut.ddl$Psi[vaut.ddl$Psi$time==c(1993)&vaut.ddl$Psi$stratum=="b"&vaut.ddl$Psi$tostratum=="G",]))
psi1993.bG=as.numeric(row.names(vaut.ddl$Psi[vaut.ddl$Psi$time==c(1994)&vaut.ddl$Psi$stratum=="b"&vaut.ddl$Psi$tostratum=="G",]))
psi1994.bG=as.numeric(row.names(vaut.ddl$Psi[vaut.ddl$Psi$time==c(1995)&vaut.ddl$Psi$stratum=="b"&vaut.ddl$Psi$tostratum=="G",]))
psi1995.bG=as.numeric(row.names(vaut.ddl$Psi[vaut.ddl$Psi$time==c(1996)&vaut.ddl$Psi$stratum=="b"&vaut.ddl$Psi$tostratum=="G",]))
psi1996.bG=as.numeric(row.names(vaut.ddl$Psi[vaut.ddl$Psi$time==c(1997)&vaut.ddl$Psi$stratum=="b"&vaut.ddl$Psi$tostratum=="G",]))
psi1999.bG=as.numeric(row.names(vaut.ddl$Psi[vaut.ddl$Psi$time==c(2000)&vaut.ddl$Psi$stratum=="b"&vaut.ddl$Psi$tostratum=="G",]))
psi2001.bG=as.numeric(row.names(vaut.ddl$Psi[vaut.ddl$Psi$time==c(2002)&vaut.ddl$Psi$stratum=="b"&vaut.ddl$Psi$tostratum=="G",]))
psi2004.bG=as.numeric(row.names(vaut.ddl$Psi[vaut.ddl$Psi$time==c(2005)&vaut.ddl$Psi$stratum=="b"&vaut.ddl$Psi$tostratum=="G",]))
psi2008.bG=as.numeric(row.names(vaut.ddl$Psi[vaut.ddl$Psi$time==c(2008)&vaut.ddl$Psi$stratum=="b"&vaut.ddl$Psi$tostratum=="G",]))

 vaut.ddl.age=vaut.ddl
  vaut.ddl.age=add.design.data(vaut.process.age,vaut.ddl.age, parameter="S", type="age", bins=c(0,1,28,37),name="age_covar",right=FALSE,replace=TRUE)
  vaut.ddl.age$S$a0=0
  vaut.ddl.age$S$a0[vaut.ddl.age$S$Age==0]=1
  vaut.ddl.age$S$a1_27=0
  vaut.ddl.age$S$a1_27[vaut.ddl.age$S$age_covar=="[1,28)"]=1
  vaut.ddl.age$S$a28_37=0
  vaut.ddl.age$S$a28_37[vaut.ddl.age$S$age_covar=="[28,37]"]=1
  vaut.ddl.age$S$a1_37=0
  vaut.ddl.age$S$a1_37[vaut.ddl.age$S$age_covar=="[28,37]"|vaut.ddl.age$S$age_covar=="[1,28)"]=1
  vaut.ddl.age$S$a0_37=0
  vaut.ddl.age$S$a0_37[vaut.ddl.age$S$Age==0|vaut.ddl.age$S$age_covar=="[28,37]"|vaut.ddl.age$S$age_covar=="[1,28)"]=1

################## relevel intercept ###################
vaut.ddl.age$p$time=relevel(vaut.ddl.age$p$time,"2007")
vaut.ddl.age$p$stratum=relevel(vaut.ddl.age$p$stratum,"G")
vaut.ddl.age$S$time=relevel(vaut.ddl.age$S$time,"2007")
vaut.ddl.age$r$time=relevel(vaut.ddl.age$r$time,"2007")
#####################################################"

# if (a==1){S.time=list(formula=~-1+a0+a1_37:time)}
# if (a==2){S.age=list(formula=~age)}
# if (a==3){S.timeXgpe=list(formula=~-1+a0+a1_37:time:fond_nl)}
# if (a==4){S.dot=list(formula=~1)}
# if (a==5){S.group=list(formula=~-1+a0+a1_37:fond_nl)}
# if (a==6){S.AgePLUStime=list(formula=~age+time)}
# if (a==7){S.ageXfond=list(formula=~age*fond_nl)}
# if (a==8){S.agePLUSfond=list(formula=~age+fond_nl)}
# if (a==9){S.timePLUSfond=list(formula=~-1+a0+a1_37:time+fond_nl)}
# if (a==10){S.ageplusfondplustime=list(formula=~age+fond_nl+time)}
# if (a==11){S.ageplusfondXtime=list(formula=~(age*fond_nl)+time)}

S.ageXtime=list(formula=~-1+a0+a1_37)
S.ageXfondXtime=list(formula=~age)


p.stratumXtime=list(formula=~stratum+time)

Psi.stratumXtime=list(formula=~-1+stratum:tostratum:fond_nl,fixed=list(index=c(psi1980.bG,psi1981.bG,psi1982.bG,psi1983.bG,psi1984.bG,psi1985.bG,psi1986.bG,psi1989.bG,psi1990.bG,psi1991.bG,psi1992.bG,psi1993.bG,psi1994.bG,psi1995.bG,psi1996.bG,psi1999.bG,psi2001.bG,psi2004.bG,psi2008.bG),value=0))

r.time=list(formula=~1)

  model.list=create.model.list("MSLiveDead")
  return(results=mark.wrapper.parallel(model.list,data=vaut.process.age,ddl=vaut.ddl.age,parallel=T,cpus=2,adjust=TRUE))
}
s=s.models()

save.image(paste(a,"_s.Rdata",sep=""))
#save.image("s_age.Rdata")
#adjust.chat(s,chat=2.36)
q()
n
