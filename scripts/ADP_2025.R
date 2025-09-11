#from Cayman_2018
#also, Cayman_Main.R (not Rmd)

# Install libraries if necessary and load them into the environment
if (suppressWarnings(require("velociraptr"))==FALSE) {
    install.packages("velociraptr",repos="https://cloud.r-project.org");
    library("velociraptr");
    }

# Install libraries if necessary and load them into the environment
if (suppressWarnings(require("vegan"))==FALSE) {
    install.packages("vegan",repos="https://cloud.r-project.org");
    library("vegan");
    }

# Install libraries if necessary and load them into the environment
if (suppressWarnings(require("betapart"))==FALSE) {
    install.packages("betapart",repos="https://cloud.r-project.org");
    library("betapart");
    }

# Install libraries if necessary and load them into the environment
# Gross, deeply disappointing to have to install it... - AZ
if (suppressWarnings(require("tidyverse"))==FALSE) {
    install.packages("tidyverse",repos="https://cloud.r-project.org");
    library("tidyverse");
    }

# Install libraries if necessary and load them into the environment
if (suppressWarnings(require("lattice"))==FALSE) {
    install.packages("lattice",repos="https://cloud.r-project.org");
    library("lattice");
    }

#importing shell data####
#setwd("~/Documents/2_Research/1_Grand_Cayman_2018/Data")
#Cay.raw.preT<-read.csv(file="Cayman_2018.csv",row.names=1)
# Changed to using github links so that working directory is no longer an issue - AZ


#dead
  Cay.raw.preT<-read.csv(file="https://raw.githubusercontent.com/aazaff/antilles_biodiversity/refs/heads/master/data/Cayman_Data_18.csv",row.names=1)
  Cay.raw.preMNI<-as.data.frame(t(Cay.raw.preT))#raw values, transformed so species are columns and samples are rows
  colnames(Cay.raw.preMNI)<-str_replace(colnames(Cay.raw.preMNI),"_"," ")#replace "_" with " " in colnames; names don't match perfectly with Cay.fam, but they are the same (spelling differences only)
#live
  CayL.raw.preT<-read.csv(file="https://raw.githubusercontent.com/aazaff/antilles_biodiversity/refs/heads/master/data/Cayman_Data_Live_18.csv",row.names=1)
  CayL.raw.preMNI<-as.data.frame(t(CayL.raw.preT))#raw values, transformed so species are columns and samples are rows
  colnames(CayL.raw.preMNI)<-str_replace(colnames(Cay.raw.preMNI),"_"," ")#replace "_" with " " in colnames; names don't match perfectly with Cay.fam, but they are the same (spelling differences only)

#MNI and standardized dataset
#dead
  Cay.raw<-cbind(Cay.raw.preMNI[,1:67]/2,Cay.raw.preMNI[,68:ncol(Cay.raw.preMNI)])#MNI estimate for bivalves, dividing by 2
  Cay.s<-decostand(Cay.raw,method="total",na.rm=T)#relative abundances (percent transformed) so samples = 1
#live
  CayL.raw<-CayL.raw.preMNI#don't need MNI for live; could correct above, but more consistent across live/dead to keep separate
  CayL.s<-decostand(CayL.raw,method="total",na.rm=T)#relative abundances (percent transformed) so samples = 1

#Grouping
Grass<-which(substr(rownames(Cay.raw),3,3)=="G")
Sand<-which(substr(rownames(Cay.raw),3,3)=="S")
RP<-which(substr(rownames(Cay.raw),1,2)=="RP")
SS<-which(substr(rownames(Cay.raw),1,2)=="SS")
EB<-which(substr(rownames(Cay.raw),1,2)=="EB")
RPG<-which(substr(rownames(Cay.raw),1,3)=="RPG")
RPS<-which(substr(rownames(Cay.raw),1,3)=="RPS")
SSG<-which(substr(rownames(Cay.raw),1,3)=="SSG")
SSS<-which(substr(rownames(Cay.raw),1,3)=="SSS")
EBG<-which(substr(rownames(Cay.raw),1,3)=="EBG")
EBS<-which(substr(rownames(Cay.raw),1,3)=="EBS")
#sample.colors<-c(rep(71,length(RPG)),rep(68,length(RPS)),rep(35,length(SSG)),rep(33,length(SSS)),rep(89,length(EBG)),rep(86,length(EBS)))
#sample.col<-colors()[sample.colors]

#Dataframe columns combined by site
#dead
  Cay.raw.site<-t(data.frame(RPG=colSums(Cay.raw[RPG,]),RPS=colSums(Cay.raw[RPS,]),SSG=colSums(Cay.raw[SSG,]),SSS=colSums(Cay.raw[SSS,]),EBG=colSums(Cay.raw[EBG,]),EBS=colSums(Cay.raw[EBS,])))
#live
  CayL.raw.site<-t(data.frame(RPG=colSums(CayL.raw[RPG,]),RPS=colSums(CayL.raw[RPS,]),SSG=colSums(CayL.raw[SSG,]),SSS=colSums(CayL.raw[SSS,]),EBG=colSums(CayL.raw[EBG,]),EBS=colSums(CayL.raw[EBS,])))

#combining dead and live by site, with family names (Table 2)  
#Cay.fam<-read.csv("Cayman_Data_18_WithFamilies_new.csv",header=T,row.names=3)
Cay.fam<-read.csv(file="https://raw.githubusercontent.com/aazaff/antilles_biodiversity/refs/heads/master/data/Fam_Cay_newest.csv",header=T,row.names=3)
DL.bysite<-t(rbind(Cay.raw.site,CayL.raw.site))
DL.fam<-data.frame(Class=Cay.fam$Class,Family=Cay.fam$Family,DL.bysite) #MNI bivalves, live and dead, integers
# write.csv(DL.fam,"DL.fam.new.csv")

#Basic Diversity Statistics####
#summary
  sum(Cay.raw.preMNI)#total elements
  sum(Cay.raw.preMNI[,68:ncol(Cay.raw.preMNI)])#total gastropods
  sum(Cay.raw.preMNI[,1:67])#total bivalves
  biv.fam<-unique(Cay.fam$Family[1:67])#total bivalve families
  gast.fam<-unique(Cay.fam$Family[68:168])#total bivalve families
#dead diversity metrics
  Richness_all  <- specnumber(Cay.raw)
  Richness_grass <- specnumber(Cay.raw[Grass, ])
  Richness_sand  <- specnumber(Cay.raw[Sand, ])
  
  TrueD_all  <- exp(diversity(Cay.raw))
  TrueD_grass <- exp(diversity(Cay.raw[Grass, ]))
  TrueD_sand  <- exp(diversity(Cay.raw[Sand, ]))
  
  Evenness_all  <- TrueD_all / Richness_all
  Evenness_grass <- TrueD_grass / Richness_grass
  Evenness_sand  <- TrueD_sand / Richness_sand
  
  Div.df <- data.frame(
    substrate = c(rep("All", nrow(Cay.raw)),
                  rep("Grass", nrow(Cay.raw[Grass, ])),
                  rep("Sand", nrow(Cay.raw[Sand, ]))),
    Richness = c(Richness_all, Richness_grass, Richness_sand),
    TrueDiv = c(TrueD_all,TrueD_grass,TrueD_sand),
    Evenness = c(Evenness_all, Evenness_grass, Evenness_sand)
  )
  
#live diversity metrics
  RichnessL_all  <- specnumber(CayL.raw)
  RichnessL_grass <- specnumber(CayL.raw[Grass, ])
  RichnessL_sand  <- specnumber(CayL.raw[Sand, ])
  
  TrueDL_all  <- exp(diversity(CayL.raw))
  TrueDL_grass <- exp(diversity(CayL.raw[Grass, ]))
  TrueDL_sand  <- exp(diversity(CayL.raw[Sand, ]))
  
  EvennessL_all  <- TrueDL_all / RichnessL_all
  EvennessL_grass <- TrueDL_grass / RichnessL_grass
  EvennessL_sand  <- TrueDL_sand / RichnessL_sand
  
  DivL.df <- data.frame(
    substrate = c(rep("All", nrow(CayL.raw)),
                  rep("Grass", nrow(CayL.raw[Grass, ])),
                  rep("Sand", nrow(CayL.raw[Sand, ]))),
    Richness = c(RichnessL_all, RichnessL_grass, RichnessL_sand),
    TrueDiv = c(TrueDL_all,TrueDL_grass,TrueDL_sand),
    Evenness = c(EvennessL_all, EvennessL_grass, EvennessL_sand)
  )
  
  #combined dead and live
    Diversity.df<-data.frame(D_L=c(rep("Dead",nrow(Div.df)),rep("Live",nrow(DivL.df))),site=rep(c(rownames(Cay.raw),rownames(Cay.raw[Grass,]),rownames(Cay.raw[Sand,])),2),locality=rep(c(substr(rownames(Cay.raw),1,3),substr(rownames(Cay.raw[Grass,]),1,3),substr(rownames(Cay.raw[Sand,]),1,3)),2),rbind(Div.df,DivL.df))
    Diversity.df.long <- gather(Diversity.df, key = "measure", value = "value", c("Richness","TrueDiv", "Evenness"))
    Diversity.df.long$locality <- factor(Diversity.df.long$locality, levels = c("RPG", "RPS", "EBG", "EBS", "SSG", "SSS"))
    Diversity.df.long$measure <- factor(Diversity.df.long$measure, levels = c("Richness","TrueDiv", "Evenness"),labels = c("Richness","True Diversity", "Evenness"))
    
    ggplot(Diversity.df.long, aes(x=locality, y=value)) + 
      geom_boxplot(aes(fill=D_L),size=0.25) +
      labs(fill = "Assemblage Type",cex=0.5) +
      theme(
        axis.text = element_text(size = 8),         # x and y axis labels
        axis.title = element_text(size = 8),        # x and y axis titles
        strip.text = element_text(size = 9),        # facet labels
        legend.text = element_text(size = 9),       # legend labels
        legend.title = element_text(size = 9)       # legend title
      ) +
      facet_wrap(~ measure,nrow=3,ncol=1,scale="free")
  
#sample completeness
#dead
  collcurve<-estaccumR(round(Cay.raw,0))
  collcurve.RPG<-estaccumR(round(Cay.raw[RPG,],0))
  collcurve.RPS<-estaccumR(round(Cay.raw[RPS,],0))
  collcurve.SSG<-estaccumR(round(Cay.raw[SSG,],0))
  collcurve.SSS<-estaccumR(round(Cay.raw[SSS,],0))
  collcurve.EBG<-estaccumR(round(Cay.raw[EBG,],0))
  collcurve.EBS<-estaccumR(round(Cay.raw[EBS,],0))
  #completeness of sampling
  collcurve$means[7,2]/collcurve$means[7,3]
  collcurve.RPG$means[7,2]/collcurve.RPG$means[7,3]#completeness %
  collcurve.RPS$means[7,2]/collcurve.RPS$means[7,3]#completeness %
  collcurve.SSG$means[7,2]/collcurve.SSG$means[7,3]#completeness %
  collcurve.SSS$means[7,2]/collcurve.SSS$means[7,3]#completeness %
  collcurve.EBG$means[7,2]/collcurve.EBG$means[7,3]#completeness %
  collcurve.EBS$means[7,2]/collcurve.EBS$means[7,3]#completeness %
  
  
  #plots
  par(mfrow=c(3,2))
  #plot(collcurve)
  collRPG<-plot(collcurve.RPG,display=c("S","chao"),ylim=c(20,155))
  collRPS<-plot(collcurve.RPS,display=c("S","chao"),ylim=c(20,155))
  collSSG<-plot(collcurve.SSG,display=c("S","chao"),ylim=c(20,155))
  collSSS<-plot(collcurve.SSS,display=c("S","chao"),ylim=c(20,155))
  collEBG<-plot(collcurve.EBG,display=c("S","chao"),ylim=c(20,155))
  collEBS<-plot(collcurve.EBS,display=c("S","chao"),ylim=c(20,155))
  

  print(collRPG,position=c(0.01,0.7,0.49,0.99),more=T)#topleft
  print(collSSG,position=c(0.01,0.4,0.49,0.69),more=T)#midleft
  print(collEBG,position=c(0.01,0.1,0.49,0.39),more=T)#lowleft
  print(collRPS,position=c(0.51,0.7,0.99,0.99),more=T)#topright
  print(collSSS,position=c(0.51,0.4,0.99,0.69),more=T)#midright
  print(collEBS,position=c(0.51,0.1,0.99,0.39))#lowright
  
#live
  collcurveL<-estaccumR(round(CayL.raw,0))
  collcurveL.RPG<-estaccumR(round(CayL.raw[RPG,],0))
  collcurveL.RPS<-estaccumR(round(CayL.raw[RPS,],0))
  collcurveL.SSG<-estaccumR(round(CayL.raw[SSG,],0))
  collcurveL.SSS<-estaccumR(round(CayL.raw[SSS,],0))
  collcurveL.EBG<-estaccumR(round(CayL.raw[EBG,],0))
  collcurveL.EBS<-estaccumR(round(CayL.raw[EBS,],0))
  #completeness of sampling
  collcurveL$means[7,2]/collcurveL$means[7,3]
  collcurveL.RPG$means[7,2]/collcurveL.RPG$means[7,3]#completeness %
  collcurveL.RPS$means[7,2]/collcurveL.RPS$means[7,3]#completeness %
  collcurveL.SSG$means[7,2]/collcurveL.SSG$means[7,3]#completeness %
  collcurveL.SSS$means[7,2]/collcurveL.SSS$means[7,3]#completeness %
  collcurveL.EBG$means[7,2]/collcurveL.EBG$means[7,3]#completeness %
  collcurveL.EBS$means[7,2]/collcurveL.EBS$means[7,3]#completeness %
  #plots - don't work because I can't get rid of ACE...
  #plot(collcurveL,display=c("S","chao"))
  par(mfrow=c(3,2))
  collRPGL<-plot(collcurveL.RPG,display=c("S","chao"),ylim=c(0,85))
  collRPSL<-plot(collcurveL.RPS,display=c("S","chao"),ylim=c(0,85))
  collSSGL<-plot(collcurveL.SSG,display=c("S","chao"),ylim=c(0,85))
  collSSSL<-plot(collcurveL.SSS,display=c("S","chao"),ylim=c(0,85))
  collEBGL<-plot(collcurveL.EBG,display=c("S","chao"),ylim=c(0,85))
  collEBSL<-plot(collcurveL.EBS,display=c("S","chao"),ylim=c(0,85))
  
  print(collRPGL,position=c(0.01,0.7,0.49,0.99),more=T)#topleft
  print(collSSGL,position=c(0.01,0.4,0.49,0.69),more=T)#midleft
  print(collEBGL,position=c(0.01,0.1,0.49,0.39),more=T)#lowleft
  print(collRPSL,position=c(0.51,0.7,0.99,0.99),more=T)#topright
  print(collSSSL,position=c(0.51,0.4,0.99,0.69),more=T)#midright
  print(collEBSL,position=c(0.51,0.1,0.99,0.39))#lowright
##somehow the summary() function in Cayman_main estimates completeness
##plots were made here and then combined in Illustrator
### lol - AZ

######################################### ADP RICHNESS #########################################
# load data
ADP.groups<-read.csv(file="https://raw.githubusercontent.com/aazaff/antilles_biodiversity/refs/heads/master/data/ADP_groups.csv",header=T,row.names=1)

#dead
Cay.raw.nosingles.int<-round(Cay.raw,0)#integer values,singletons removed for ADP analysis
Cay.raw.nosingles<-Cay.raw.nosingles.int[,colSums(Cay.raw.nosingles.int)>1]#integer values,singletons removed for ADP analysis
Cay.ADP.S<-adipart(Cay.raw.nosingles ~ .,ADP.groups,index="richness",nsimul=999,method="r2dtable",na.rm=T)#ADP Richness

#live
  CayL.raw.int<-round(CayL.raw,0)#integer values,singletons NOT excluded
  CayL.ADP.S<-adipart(CayL.raw.int ~ .,ADP.groups,index="richness",nsimul=999,method="r2dtable",na.rm=T)#ADP Richness

#weighted by sample abundance (Don't want/need these when using Hill numbers)
#dead
#  Cay.ADP.S.weight<-adipart(Cay.raw.nosingles ~ .,ADP.groups,index="richness",nsimul=999,method="r2dtable",na.rm=T,weight="prop")#ADP Richness
#live
# CayL.ADP.S.weight<-adipart(CayL.raw.int ~ .,ADP.groups,index="richness",nsimul=999,method="r2dtable",na.rm=T,weight="prop")#ADP Richness

#split by grass and sand
#dead
  Cay.ADP.S.grass<-adipart(Cay.raw.nosingles[c(RPG,SSG,EBG),] ~ .,ADP.groups[c(RPG,SSG,EBG),],index="richness",nsimul=999,method="r2dtable",na.rm=T)
  Cay.ADP.S.sand<-adipart(Cay.raw.nosingles[c(RPS,SSS,EBS),] ~ .,ADP.groups[c(RPS,SSS,EBS),],index="richness",nsimul=999,method="r2dtable",na.rm=T)
  #weighted by sample abundance (Don't want/need these when using Hill numbers)
#  Cay.ADP.S.grass.weight<-adipart(Cay.raw.nosingles[c(RPG,SSG,EBG),] ~ .,ADP.groups[c(RPG,SSG,EBG),],index="richness",nsimul=999,method="r2dtable",na.rm=T,weight="prop")
#  Cay.ADP.S.sand.weight<-adipart(Cay.raw.nosingles[c(RPS,SSS,EBS),] ~ .,ADP.groups[c(RPS,SSS,EBS),],index="richness",nsimul=999,method="r2dtable",na.rm=T,weight="prop")
#live
  CayL.ADP.S.grass<-adipart(CayL.raw.int[c(RPG,SSG,EBG),] ~ .,ADP.groups[c(RPG,SSG,EBG),],index="richness",nsimul=999,method="r2dtable",na.rm=T)
  CayL.ADP.S.sand<-adipart(CayL.raw.int[c(RPS,SSS,EBS),] ~ .,ADP.groups[c(RPS,SSS,EBS),],index="richness",nsimul=999,method="r2dtable",na.rm=T)
  #weighted by sample abundance (Don't want/need these when using Hill numbers)
#  CayL.ADP.S.grass.weight<-adipart(CayL.raw.int[c(RPG,SSG,EBG),] ~ .,ADP.groups[c(RPG,SSG,EBG),],index="richness",nsimul=999,method="r2dtable",na.rm=T,weight="prop")
#  CayL.ADP.S.sand.weight<-adipart(CayL.raw.int[c(RPS,SSS,EBS),] ~ .,ADP.groups[c(RPS,SSS,EBS),],index="richness",nsimul=999,method="r2dtable",na.rm=T,weight="prop")
  

#DEAD: Master ADP table with all analyses (REMOVE WEIGHTED)
ADP.Table.S<-data.frame(D_L=rep("Dead",18),
                      sites=c(rep("all",6),rep("grass",6),rep("sand",6)),#36
                      #sites=c(rep("all",12),rep("grass",12),rep("sand",12)),#36
                      method=rep("Richness",18),
                      scale=rep(c("a1","B1","B2","B3","B4","g"),3), #36
                      #scale=rep(c("a1","B1","B2","B3","B4","g"),6), #36
                      #weight=rep(c(rep("uniform",6),rep("weighted",6)),3),#36
                      observed=round(c(Cay.ADP.S$statistic[c(c(1,6:9,5))],
                                       #Cay.ADP.S.weight$statistic[c(1,6:9,5)],
                                       Cay.ADP.S.grass$statistic[c(1,6:9,5)],
                                       #Cay.ADP.S.grass.weight$statistic[c(1,6:9,5)],
                                       Cay.ADP.S.sand$statistic[c(1,6:9,5)]#,
                                       #Cay.ADP.S.sand.weight$statistic[c(1,6:9,5)]
                                       )
                                    ),
                      percent=c(round(Cay.ADP.S$statistic[c(1,6:9,5)]/Cay.ADP.S$statistic[5]*100,2),
                                    #round(Cay.ADP.S.weight$statistic[c(1,6:9,5)]/Cay.ADP.S.weight$statistic[5]*100,2),
                                    round(Cay.ADP.S.grass$statistic[c(1,6:9,5)]/Cay.ADP.S.grass$statistic[5]*100,2),
                                    #round(Cay.ADP.S.grass.weight$statistic[c(1,6:9,5)]/Cay.ADP.S.grass.weight$statistic[5]*100,2),
                                    round(Cay.ADP.S.sand$statistic[c(1,6:9,5)]/Cay.ADP.S.sand$statistic[5]*100,2)#,
                                    #round(Cay.ADP.S.sand.weight$statistic[c(1,6:9,5)]/Cay.ADP.S.sand.weight$statistic[5]*100,2)
                                ),
                      expected=round(c(Cay.ADP.S$oecosimu$means[c(1,6:9,5)],
                                    #Cay.ADP.S.weight$oecosimu$means[c(1,6:9,5)],
                                    Cay.ADP.S.grass$oecosimu$means[c(1,6:9,5)],
                                    #Cay.ADP.S.grass.weight$oecosimu$means[c(1,6:9,5)],
                                    Cay.ADP.S.sand$oecosimu$means[c(1,6:9,5)]#,
                                    #Cay.ADP.S.sand.weight$oecosimu$means[c(1,6:9,5)]
                                      )
                                    ),
                      CI2.5=round(c(apply(Cay.ADP.S$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],
                                    #apply(Cay.ADP.S.weight$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],
                                    apply(Cay.ADP.S.grass$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],
                                    #apply(Cay.ADP.S.grass.weight$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],
                                    apply(Cay.ADP.S.sand$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)]#,
                                    #apply(Cay.ADP.S.sand.weight$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)]
                                    )
                                  ),
                      CI97.5=round(c(apply(Cay.ADP.S$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],
                                    #apply(Cay.ADP.S.weight$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],
                                    apply(Cay.ADP.S.grass$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],
                                    #apply(Cay.ADP.S.grass.weight$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],
                                    apply(Cay.ADP.S.sand$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)]#,
                                    #apply(Cay.ADP.S.sand.weight$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)]
                                    )
                                  ),
                      pval=c(Cay.ADP.S$oecosimu$pval[c(1,6:9,NA)],
                                    #Cay.ADP.S.weight$oecosimu$pval[c(1,6:9,NA)],
                                    Cay.ADP.S.grass$oecosimu$pval[c(1,6:9,NA)],
                                    #Cay.ADP.S.grass.weight$oecosimu$pval[c(1,6:9,NA)],
                                    Cay.ADP.S.sand$oecosimu$pval[c(1,6:9,NA)]#,
                                    #Cay.ADP.S.sand.weight$oecosimu$pval[c(1,6:9,NA)]
                             )
)

# write.csv(ADP.Table.S,"Full_Dead_ADP_Table_new.csv")

#Master ADP table in long form for faceting in ggplot only
ADP.Table.long<-gather(ADP.Table,key="measure",value="value",c("observed","expected"))

# LIVE: Master ADP table with alldead analyses
ADPL.Table.S<-data.frame(D_L=rep("Live",18),
                      sites=c(rep("all",6),rep("grass",6),rep("sand",6)),
                      method=rep("Richness",18),
                      scale=rep(c("a1","B1","B2","B3","B4","g"),3),
                      observed=round(c(CayL.ADP.S$statistic[c(c(1,6:9,5))],CayL.ADP.S.grass$statistic[c(1,6:9,5)],CayL.ADP.S.sand$statistic[c(1,6:9,5)])),#36
                      percent=c(round(CayL.ADP.S$statistic[c(1,6:9,5)]/CayL.ADP.S$statistic[5]*100,2),round(CayL.ADP.S.grass$statistic[c(1,6:9,5)]/CayL.ADP.S.grass$statistic[5]*100,2),round(CayL.ADP.S.sand$statistic[c(1,6:9,5)]/CayL.ADP.S.sand$statistic[5]*100,2)),#36
                      expected=round(c(CayL.ADP.S$oecosimu$means[c(1,6:9,5)],CayL.ADP.S.grass$oecosimu$means[c(1,6:9,5)],CayL.ADP.S.sand$oecosimu$means[c(1,6:9,5)])),#36
                      CI2.5=round(c(apply(CayL.ADP.S$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],apply(CayL.ADP.S.grass$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],apply(CayL.ADP.S.sand$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)])),#36  
                      CI97.5=round(c(apply(CayL.ADP.S$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],apply(CayL.ADP.S.grass$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],apply(CayL.ADP.S.sand$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)])),
                      pval=c(CayL.ADP.S$oecosimu$pval[c(1,6:9,NA)],CayL.ADP.S.grass$oecosimu$pval[c(1,6:9,NA)],CayL.ADP.S.sand$oecosimu$pval[c(1,6:9,NA)])#36
)
# write.csv(ADPL.Table.S,"Full_Live_ADP_Table_new.csv")

#Master ADPL table in long form for faceting in ggplot only
ADPL.Table.long<-gather(ADPL.Table.S,key="measure",value="value",c("observed","expected"))

Rich.ADP <- rbind(ADP.Table.S,ADPL.Table.S)#all dead and live ADP values, RICHNESS

#plotting ADP####

#ALL PLOTS IN ONE FIGURE!
#dead
ggplot(ADP.Table.long[!ADP.Table.long$scale=="g",], aes(fill=factor(ADP.Table.long[!ADP.Table.long$scale=="g",]$scale,levels=c("B4","B3","B2","B1","a1")),x=factor(measure,levels=c("observed","expected")), y=value)) +
  #geom_col() + #raw values
  geom_bar(position="fill",stat="identity") + #percentages
  labs(fill="partition level",x="",title="") +
  #scale_y_continuous(name="value", breaks=seq(0,1,.2)) + #uncomment if using percentages
  scale_x_discrete(labels=c("obs","exp")) +
  #facet_wrap(. ~ sites + weight,nrow=1,ncol=6)
  facet_wrap(. ~ sites,nrow=1,ncol=3)
#live
ggplot(ADPL.Table.long[!ADPL.Table.long$scale=="g",], aes(fill=factor(ADPL.Table.long[!ADPL.Table.long$scale=="g",]$scale,levels=c("B5","B4","B3","B2","B1","a1")),x=factor(measure,levels=c("observed","expected")), y=value)) +
  #geom_col() + #raw values
  geom_bar(position="fill",stat="identity") + #percentages
  labs(fill="partition level",x="",title="") +
  #scale_y_continuous(name="value", breaks=seq(0,1,.2)) + #uncomment if using percentages
  scale_x_discrete(labels=c("obs","exp")) +
  facet_wrap(. ~ sites + weight,nrow=1,ncol=6)

#############################################################################################################
############################################ MANUAL CHECK OF ADP  ###########################################
#############################################################################################################
# Calculate Alpha of jth taxon
taxonAlpha = function(CommunityMatrix) {
	Nj = apply(CommunityMatrix, MARGIN=2, FUN=function(x) length(which(x>0))) # this works, but is malformed and brittle... come back and do it right!
	Samples = nrow(CommunityMatrix)
	Alphaj = Nj/Samples
	names(Alphaj) = colnames(CommunityMatrix)
	return(Alphaj)
	}

# Calculate total alpha diversity
meanAlpha = function (CommunityMatrix) {
	return(sum(taxonAlpha(CommunityMatrix)))
	}

# Calculate beta of jth taxon
taxonBeta = function(CommunityMatrix) {
	Nj = apply(CommunityMatrix, MARGIN=2, FUN=function(x) length(which(x>0))) # this works, but is malformed and brittle... come back and do it right!
	Samples = nrow(CommunityMatrix)
	Betaj = (Samples - Nj) / Samples
	names(Betaj) = colnames(CommunityMatrix)
    return(Betaj)
    }

# Calculate sample beta
sampleBeta = function(CommunityMatrix) {
	Betaj = taxonBeta(CommunityMatrix)
	Nj = apply(CommunityMatrix, MARGIN=2, FUN=function(x) length(which(x>0))) # this works, but is malformed and brittle... come back and do it right!
	NSamples = nrow(CommunityMatrix)
	Betai = vector(mode="numeric", length=NSamples)
	for (i in seq_len(NSamples)) {
		Betai[i] = sum(Betaj / Nj * CommunityMatrix[i, ])
		}
	names(Betai) = rownames(CommunityMatrix)
	return(Betai)
	}

# Calculate beta diversity
totalBeta = function (CommunityMatrix) {
	return(sum(taxonBeta(CommunityMatrix)))
	}

# Calculate gamma diversity
totalGamma = function(CommunityMatrix) {
	Gamma = length(which(apply(CommunityMatrix,2,sum)>0))
    return(Gamma)
	}

# general formula for alpha (equal weights)
genAlpha = function(CommunityMatrix,q=0) {
    N = nrow(CommunityMatrix)
    Output = vector("numeric",length=N)
    for (j in seq_len(N)) {
        Nonzero = CommunityMatrix[j,CommunityMatrix[j,]>0]
        Relative = Nonzero/sum(Nonzero)
        Sum = sum(Relative^q)
        Output[j] = (1/N)*Sum
        }
    Mean = sum(Output)
    Result = Mean^(1/(1-q))
    return(Result)
    }

# load the breakdown of data into groups
ADP.groups<-read.csv(file="https://raw.githubusercontent.com/aazaff/antilles_biodiversity/refs/heads/master/data/ADP_groups.csv",header=T,row.names=1)

# round(x,0) is used elsewhere, but since you cannot have a half an individual it should always be rounded up - i.e., so that taxa with values of 0.5 don't become zero
Dead = ceiling(Cay.raw)
Live = ceiling(CayL.raw)

# Quick sanity check of raw calculation vs. adipart for alpha and beta
meanAlpha(Dead)
totalGamma(Dead)
adipart(Dead ~ .,ADP.groups,index="richness",nsimul=999,method="r2dtable",na.rm=T)
# Quick sanity check of raw caclulation for a4 and beta4
RP = subset(Dead,substr(rownames(Dead),1,2)=="RP")
SS = subset(Dead,substr(rownames(Dead),1,2)=="SS")
EB = subset(Dead,substr(rownames(Dead),1,2)=="EB")
mean(c(totalGamma(RP),totalGamma(SS),totalGamma(EB)))
B4 = rbind(apply(RP,2,sum),apply(SS,2,sum),apply(EB,2,sum))
totalBeta(B4) # the beta is off by a value of 1 (see above note about malformed count)
# Quick sanity check of Shannon and Hill's
mean(apply(B4,1,diversity)) # manual a4
diversity(apply(Dead,2,sum)) # manual gamma
adipart(Dead ~ .,ADP.groups,index="shannon",nsimul=999,method="r2dtable",na.rm=T)
# alpha using Jost alpha formula
genAlpha(Dead,0)
genAlpha(Dead,0.99)
genAlpha(Dead,2)
genAlpha(B4,0)


#dead
  Cay.ADP.S.grass<-adipart(Cay.raw.nosingles[c(RPG,SSG,EBG),] ~ .,ADP.groups[c(RPG,SSG,EBG),],index="richness",nsimul=999,method="r2dtable",na.rm=T)
  Cay.ADP.S.sand<-adipart(Cay.raw.nosingles[c(RPS,SSS,EBS),] ~ .,ADP.groups[c(RPS,SSS,EBS),],index="richness",nsimul=999,method="r2dtable",na.rm=T)
  #weighted by sample abundance (Don't want/need these when using Hill numbers)
#  Cay.ADP.S.grass.weight<-adipart(Cay.raw.nosingles[c(RPG,SSG,EBG),] ~ .,ADP.groups[c(RPG,SSG,EBG),],index="richness",nsimul=999,method="r2dtable",na.rm=T,weight="prop")
#  Cay.ADP.S.sand.weight<-adipart(Cay.raw.nosingles[c(RPS,SSS,EBS),] ~ .,ADP.groups[c(RPS,SSS,EBS),],index="richness",nsimul=999,method="r2dtable",na.rm=T,weight="prop")
#live
  CayL.ADP.S.grass<-adipart(CayL.raw.int[c(RPG,SSG,EBG),] ~ .,ADP.groups[c(RPG,SSG,EBG),],index="richness",nsimul=999,method="r2dtable",na.rm=T)
  CayL.ADP.S.sand<-adipart(CayL.raw.int[c(RPS,SSS,EBS),] ~ .,ADP.groups[c(RPS,SSS,EBS),],index="richness",nsimul=999,method="r2dtable",na.rm=T)
  #weighted by sample abundance (Don't want/need these when using Hill numbers)
#  CayL.ADP.S.grass.weight<-adipart(CayL.raw.int[c(RPG,SSG,EBG),] ~ .,ADP.groups[c(RPG,SSG,EBG),],index="richness",nsimul=999,method="r2dtable",na.rm=T,weight="prop")
#  CayL.ADP.S.sand.weight<-adipart(CayL.raw.int[c(RPS,SSS,EBS),] ~ .,ADP.groups[c(RPS,SSS,EBS),],index="richness",nsimul=999,method="r2dtable",na.rm=T,weight="prop")
  

#DEAD: Master ADP table with all analyses (REMOVE WEIGHTED)
ADP.Table.S<-data.frame(D_L=rep("Dead",18),
                      sites=c(rep("all",6),rep("grass",6),rep("sand",6)),#36
                      #sites=c(rep("all",12),rep("grass",12),rep("sand",12)),#36
                      method=rep("Richness",18),
                      scale=rep(c("a1","B1","B2","B3","B4","g"),3), #36
                      #scale=rep(c("a1","B1","B2","B3","B4","g"),6), #36
                      #weight=rep(c(rep("uniform",6),rep("weighted",6)),3),#36
                      observed=round(c(Cay.ADP.S$statistic[c(c(1,6:9,5))],
                                       #Cay.ADP.S.weight$statistic[c(1,6:9,5)],
                                       Cay.ADP.S.grass$statistic[c(1,6:9,5)],
                                       #Cay.ADP.S.grass.weight$statistic[c(1,6:9,5)],
                                       Cay.ADP.S.sand$statistic[c(1,6:9,5)]#,
                                       #Cay.ADP.S.sand.weight$statistic[c(1,6:9,5)]
                                       )
                                    ),
                      percent=c(round(Cay.ADP.S$statistic[c(1,6:9,5)]/Cay.ADP.S$statistic[5]*100,2),
                                    #round(Cay.ADP.S.weight$statistic[c(1,6:9,5)]/Cay.ADP.S.weight$statistic[5]*100,2),
                                    round(Cay.ADP.S.grass$statistic[c(1,6:9,5)]/Cay.ADP.S.grass$statistic[5]*100,2),
                                    #round(Cay.ADP.S.grass.weight$statistic[c(1,6:9,5)]/Cay.ADP.S.grass.weight$statistic[5]*100,2),
                                    round(Cay.ADP.S.sand$statistic[c(1,6:9,5)]/Cay.ADP.S.sand$statistic[5]*100,2)#,
                                    #round(Cay.ADP.S.sand.weight$statistic[c(1,6:9,5)]/Cay.ADP.S.sand.weight$statistic[5]*100,2)
                                ),
                      expected=round(c(Cay.ADP.S$oecosimu$means[c(1,6:9,5)],
                                    #Cay.ADP.S.weight$oecosimu$means[c(1,6:9,5)],
                                    Cay.ADP.S.grass$oecosimu$means[c(1,6:9,5)],
                                    #Cay.ADP.S.grass.weight$oecosimu$means[c(1,6:9,5)],
                                    Cay.ADP.S.sand$oecosimu$means[c(1,6:9,5)]#,
                                    #Cay.ADP.S.sand.weight$oecosimu$means[c(1,6:9,5)]
                                      )
                                    ),
                      CI2.5=round(c(apply(Cay.ADP.S$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],
                                    #apply(Cay.ADP.S.weight$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],
                                    apply(Cay.ADP.S.grass$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],
                                    #apply(Cay.ADP.S.grass.weight$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],
                                    apply(Cay.ADP.S.sand$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)]#,
                                    #apply(Cay.ADP.S.sand.weight$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)]
                                    )
                                  ),
                      CI97.5=round(c(apply(Cay.ADP.S$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],
                                    #apply(Cay.ADP.S.weight$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],
                                    apply(Cay.ADP.S.grass$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],
                                    #apply(Cay.ADP.S.grass.weight$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],
                                    apply(Cay.ADP.S.sand$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)]#,
                                    #apply(Cay.ADP.S.sand.weight$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)]
                                    )
                                  ),
                      pval=c(Cay.ADP.S$oecosimu$pval[c(1,6:9,NA)],
                                    #Cay.ADP.S.weight$oecosimu$pval[c(1,6:9,NA)],
                                    Cay.ADP.S.grass$oecosimu$pval[c(1,6:9,NA)],
                                    #Cay.ADP.S.grass.weight$oecosimu$pval[c(1,6:9,NA)],
                                    Cay.ADP.S.sand$oecosimu$pval[c(1,6:9,NA)]#,
                                    #Cay.ADP.S.sand.weight$oecosimu$pval[c(1,6:9,NA)]
                             )
)

# write.csv(ADP.Table.S,"Full_Dead_ADP_Table_new.csv")

#Master ADP table in long form for faceting in ggplot only
ADP.Table.long<-gather(ADP.Table,key="measure",value="value",c("observed","expected"))

# LIVE: Master ADP table with alldead analyses
ADPL.Table.S<-data.frame(D_L=rep("Live",18),
                      sites=c(rep("all",6),rep("grass",6),rep("sand",6)),
                      method=rep("Richness",18),
                      scale=rep(c("a1","B1","B2","B3","B4","g"),3),
                      observed=round(c(CayL.ADP.S$statistic[c(c(1,6:9,5))],CayL.ADP.S.grass$statistic[c(1,6:9,5)],CayL.ADP.S.sand$statistic[c(1,6:9,5)])),#36
                      percent=c(round(CayL.ADP.S$statistic[c(1,6:9,5)]/CayL.ADP.S$statistic[5]*100,2),round(CayL.ADP.S.grass$statistic[c(1,6:9,5)]/CayL.ADP.S.grass$statistic[5]*100,2),round(CayL.ADP.S.sand$statistic[c(1,6:9,5)]/CayL.ADP.S.sand$statistic[5]*100,2)),#36
                      expected=round(c(CayL.ADP.S$oecosimu$means[c(1,6:9,5)],CayL.ADP.S.grass$oecosimu$means[c(1,6:9,5)],CayL.ADP.S.sand$oecosimu$means[c(1,6:9,5)])),#36
                      CI2.5=round(c(apply(CayL.ADP.S$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],apply(CayL.ADP.S.grass$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],apply(CayL.ADP.S.sand$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)])),#36  
                      CI97.5=round(c(apply(CayL.ADP.S$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],apply(CayL.ADP.S.grass$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],apply(CayL.ADP.S.sand$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)])),
                      pval=c(CayL.ADP.S$oecosimu$pval[c(1,6:9,NA)],CayL.ADP.S.grass$oecosimu$pval[c(1,6:9,NA)],CayL.ADP.S.sand$oecosimu$pval[c(1,6:9,NA)])#36
)
# write.csv(ADPL.Table.S,"Full_Live_ADP_Table_new.csv")

#Master ADPL table in long form for faceting in ggplot only
ADPL.Table.long<-gather(ADPL.Table.S,key="measure",value="value",c("observed","expected"))

Rich.ADP <- rbind(ADP.Table.S,ADPL.Table.S)#all dead and live ADP values, RICHNESS

#plotting ADP####

#ALL PLOTS IN ONE FIGURE!
#dead
ggplot(ADP.Table.long[!ADP.Table.long$scale=="g",], aes(fill=factor(ADP.Table.long[!ADP.Table.long$scale=="g",]$scale,levels=c("B4","B3","B2","B1","a1")),x=factor(measure,levels=c("observed","expected")), y=value)) +
  #geom_col() + #raw values
  geom_bar(position="fill",stat="identity") + #percentages
  labs(fill="partition level",x="",title="") +
  #scale_y_continuous(name="value", breaks=seq(0,1,.2)) + #uncomment if using percentages
  scale_x_discrete(labels=c("obs","exp")) +
  #facet_wrap(. ~ sites + weight,nrow=1,ncol=6)
  facet_wrap(. ~ sites,nrow=1,ncol=3)
#live
ggplot(ADPL.Table.long[!ADPL.Table.long$scale=="g",], aes(fill=factor(ADPL.Table.long[!ADPL.Table.long$scale=="g",]$scale,levels=c("B5","B4","B3","B2","B1","a1")),x=factor(measure,levels=c("observed","expected")), y=value)) +
  #geom_col() + #raw values
  geom_bar(position="fill",stat="identity") + #percentages
  labs(fill="partition level",x="",title="") +
  #scale_y_continuous(name="value", breaks=seq(0,1,.2)) + #uncomment if using percentages
  scale_x_discrete(labels=c("obs","exp")) +
  facet_wrap(. ~ sites + weight,nrow=1,ncol=6)

######################################### ADP HILL1, Using vegan >2.8 #########################################
#change nsimul back to 999 for final iteration
#r2dtable and quasiswap_count are producing negative values in Cay.ADP.S.grass...
#dead
Cay.ADP.H1<-adipart(Cay.raw.nosingles ~ .,ADP.groups,index="hill1",nsimul=999,method="r2dtable",na.rm=T)#ADP Richness
#live
CayL.ADP.H1<-adipart(CayL.raw.int ~ .,ADP.groups,index="hill1",nsimul=999,method="r2dtable",na.rm=T)#ADP Richness

#split by grass and sand
#dead
Cay.ADP.H1.grass<-adipart(Cay.raw.nosingles[c(RPG,SSG,EBG),] ~ .,ADP.groups[c(RPG,SSG,EBG),],index="hill1",nsimul=999,method="r2dtable",na.rm=T)
Cay.ADP.H1.sand<-adipart(Cay.raw.nosingles[c(RPS,SSS,EBS),] ~ .,ADP.groups[c(RPS,SSS,EBS),],index="hill1",nsimul=999,method="r2dtable",na.rm=T)
#live
CayL.ADP.H1.grass<-adipart(CayL.raw.int[c(RPG,SSG,EBG),] ~ .,ADP.groups[c(RPG,SSG,EBG),],index="hill1",nsimul=999,method="r2dtable",na.rm=T)
CayL.ADP.H1.sand<-adipart(CayL.raw.int[c(RPS,SSS,EBS),] ~ .,ADP.groups[c(RPS,SSS,EBS),],index="hill1",nsimul=999,method="r2dtable",na.rm=T)

#DEAD: Master ADP table with all analyses (Hill=1)
ADP.Table.H1<-data.frame(D_L=c(rep("Dead",18)),
                       sites=c(rep("all",6),rep("grass",6),rep("sand",6)),
                       method=rep("Hill=1",18),
                       scale=rep(c("a1","B1","B2","B3","B4","g"),3),
                       observed=round(c(Cay.ADP.H1$statistic[c(c(1,6:9,5))],Cay.ADP.H1.grass$statistic[c(1,6:9,5)],Cay.ADP.H1.sand$statistic[c(1,6:9,5)])),#36
                       percent=c(round(Cay.ADP.H1$statistic[c(1,6:9,5)]/Cay.ADP.H1$statistic[5]*100,2),round(Cay.ADP.H1.grass$statistic[c(1,6:9,5)]/Cay.ADP.H1.grass$statistic[5]*100,2),round(Cay.ADP.H1.sand$statistic[c(1,6:9,5)]/Cay.ADP.H1.sand$statistic[5]*100,2)),#36
                       expected=round(c(Cay.ADP.H1$oecosimu$means[c(1,6:9,5)],Cay.ADP.H1.grass$oecosimu$means[c(1,6:9,5)],Cay.ADP.H1.sand$oecosimu$means[c(1,6:9,5)])),#36
                       CI2.5=round(c(apply(Cay.ADP.H1$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],apply(Cay.ADP.H1.grass$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],apply(Cay.ADP.H1.sand$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)])),#36  
                       CI97.5=round(c(apply(Cay.ADP.H1$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],apply(Cay.ADP.H1.grass$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],apply(Cay.ADP.H1.sand$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)])),
                       pval=c(Cay.ADP.H1$oecosimu$pval[c(1,6:9,NA)],Cay.ADP.H1.grass$oecosimu$pval[c(1,6:9,NA)],Cay.ADP.H1.sand$oecosimu$pval[c(1,6:9,NA)])#36
)
#write.csv(ADP.Table,"Full_Dead_ADP_Table_new.csv")#use if need to export table

#Master ADP table in long form for faceting in ggplot only
ADP.Table.long.H1<-gather(ADP.Table.H1,key="measure",value="value",c("observed","expected"))

# LIVE: Master ADP table with alldead analyses (REMOVE WEIGHTED)
ADPL.Table.H1<-data.frame(D_L=rep("Live",18),
                       sites=c(rep("all",6),rep("grass",6),rep("sand",6)),
                       method=rep("Hill=1",18),
                       scale=rep(c("a1","B1","B2","B3","B4","g"),3),
                       observed=round(c(CayL.ADP.H1$statistic[c(c(1,6:9,5))],CayL.ADP.H1.grass$statistic[c(1,6:9,5)],CayL.ADP.H1.sand$statistic[c(1,6:9,5)])),#36
                       percent=c(round(CayL.ADP.H1$statistic[c(1,6:9,5)]/CayL.ADP.H1$statistic[5]*100,2),round(CayL.ADP.H1.grass$statistic[c(1,6:9,5)]/CayL.ADP.H1.grass$statistic[5]*100,2),round(CayL.ADP.H1.sand$statistic[c(1,6:9,5)]/CayL.ADP.H1.sand$statistic[5]*100,2)),#36
                       expected=round(c(CayL.ADP.H1$oecosimu$means[c(1,6:9,5)],CayL.ADP.H1.grass$oecosimu$means[c(1,6:9,5)],CayL.ADP.H1.sand$oecosimu$means[c(1,6:9,5)])),#36
                       CI2.5=round(c(apply(CayL.ADP.H1$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],apply(CayL.ADP.H1.grass$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],apply(CayL.ADP.H1.sand$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)])),#36  
                       CI97.5=round(c(apply(CayL.ADP.H1$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],apply(CayL.ADP.H1.grass$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],apply(CayL.ADP.H1.sand$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)])),
                       pval=c(CayL.ADP.H1$oecosimu$pval[c(1,6:9,NA)],CayL.ADP.H1.grass$oecosimu$pval[c(1,6:9,NA)],CayL.ADP.H1.sand$oecosimu$pval[c(1,6:9,NA)])#36
)
#write.csv(ADPL.Table,"Full_Live_ADP_Table_new.csv")#to export table

#Master ADP table in long form for faceting in ggplot only
ADPL.Table.long.H1<-gather(ADPL.Table.H1,key="measure",value="value",c("observed","expected"))

H1.ADP <- rbind(ADP.Table.H1,ADPL.Table.H1)#all dead and live ADP values, RICHNESS

ADP.Table <- rbind(Rich.ADP,H1.ADP)
ADP.Table.long<-gather(ADP.Table,key="measure",value="value",c("observed","expected"))
#ADP.Table.long$measure <- factor(ADP.Table.long$measure, levels = c("observed","expected"),labels = c("Observed","Expected"))
ADP.Table.long$method <- factor(ADP.Table.long$method, levels = c("Richness","Hill=1"),labels = c("Richness","True Diversity"))
ADP.Table.long$sites <- factor(ADP.Table.long$sites, levels = c("all","grass","sand"),labels = c("all","grass","sand"))

#plotting ADP ****MASTER****####
#ALL PLOTS IN ONE FIGURE

ggplot(ADP.Table.long[!ADP.Table.long$scale=="g",], aes(fill=factor(ADP.Table.long[!ADP.Table.long$scale=="g",]$scale,levels=c("B4","B3","B2","B1","a1")),x=factor(measure,levels=c("observed","expected")), y=value)) +
  #geom_col() + #raw values
  geom_bar(position="fill",stat="identity") + #percentages
  labs(fill="partition level",x="",title="") +
  #scale_y_continuous(name="value", breaks=seq(0,1,.2)) + #uncomment if using percentages
  scale_x_discrete(labels=c("obs","exp")) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    axis.text = element_text(size = 8),         # x and y axis labels
    axis.title = element_text(size = 8),        # x and y axis titles
    #strip.text = element_text(size = 9),        # facet labels
    strip.text.x = element_text(size = 8),
    legend.text = element_text(size = 9),       # legend labels
    legend.title = element_text(size = 9)       # legend title
  ) +
  facet_wrap(. ~ method + D_L + sites,nrow=2,ncol=6)

#Spatial Decay of Assemblage Similarity####
Cay.Site.GPS<-read.csv(file="SampleSites.csv",header=T,row.names=2)
UTM.Easting<-Cay.Site.GPS$UTM_East
UTM.Northing<-Cay.Site.GPS$UTM_North
Cay.spat.dist<-dist(data.frame(UTM.Easting,UTM.Northing,row.names=rownames(Cay.Site.GPS)),method="euclidean")#actual distance between samples
Cay.spat.dist.mat<-as.matrix(Cay.spat.dist)
#log.Cay.spat.dist<-log(Cay.spat.dist)
Cay.G.spat.dist<-dist(data.frame(UTM.Easting[Grass],UTM.Northing[Grass],row.names=rownames(Cay.Site.GPS[Grass,])))
Cay.S.spat.dist<-dist(data.frame(UTM.Easting[Sand],UTM.Northing[Sand],row.names=rownames(Cay.Site.GPS[Sand,])))

#percent transform biotic data
#dead
  Cay.std<-decostand(Cay.raw.nosingles,"total")
#live
  CayL.std<-decostand(CayL.raw.int,"total")
  CayL.std.nonzero<-CayL.std[,c(colSums(CayL.std)>0)]#getting rid of empty species columns
#calculate bray-curtis distance
#dead
  Cay.std.dist<-vegdist(Cay.std,"bray") #dissimilarity
  Cay.sim<-1-Cay.std.dist #similarity
  Cay.G.sim<-1-(vegdist(Cay.std[Grass,],"bray"))#similarity of seagrass samples
  Cay.S.sim<-1-(vegdist(Cay.std[Sand,]))#similarity of sand samples
#live
  CayL.std.dist<-vegdist(CayL.std.nonzero,"bray") #dissimilarity
  CayL.sim<-1-CayL.std.dist #similarity
  CayL.G.sim<-1-(vegdist(CayL.std.nonzero[Grass,]))#similarity of grass samples
  CayL.S.sim<-1-(vegdist(CayL.std.nonzero[Sand,]))#similarity of sand samples
  
#Start here - this works for modeling distance decay-->
  #dead, grass and sand separately
    dist.sim.df<-data.frame(G_S=c(rep("G",length(Cay.G.spat.dist)),rep("S",length(Cay.G.spat.dist))),dist=c(Cay.G.spat.dist,Cay.S.spat.dist),sim=c(Cay.G.sim,Cay.S.sim))
    grass.lm.slope<-round(coef(lm(dist.sim.df$sim[which(dist.sim.df$G_S=="G")]~dist.sim.df$dist[which(dist.sim.df$G_S=="G")]))[2],7)#slope of the lm dead grass
    sand.lm.slope<-round(coef(lm(dist.sim.df$sim[which(dist.sim.df$G_S=="S")]~dist.sim.df$dist[which(dist.sim.df$G_S=="S")]))[2],6)#slope of the lm dead sand
    
    ggplot(data = dist.sim.df, aes(x=dist,y=sim)) +
      #geom_point(col="red") +
      geom_point(aes(color=factor(G_S))) +
      geom_smooth(method="lm", aes(color=factor(G_S)), fill = "gray") + #doesn't work yet
      ylim(0,1) +
      geom_label(aes(x=22000,y=0.95),label=paste("slope vegetated =",grass.lm.slope),show.legend=F) +
      geom_label(aes(x=22000,y=0.85),label=paste("slope unvegetated =",sand.lm.slope),show.legend=F) +
      scale_color_manual(values=c("G" = "red",
                                  "S" = "blue"))
      
  #live, grass and sand separate
    distL.sim.df<-data.frame(G_S=c(rep("G",length(Cay.G.spat.dist)),rep("S",length(Cay.G.spat.dist))),dist=c(Cay.G.spat.dist,Cay.S.spat.dist),sim=c(CayL.G.sim,CayL.S.sim))
    grassL.lm.slope<-round(coef(lm(distL.sim.df$sim[which(distL.sim.df$G_S=="G")]~distL.sim.df$dist[which(distL.sim.df$G_S=="G")]))[2],7)#slope of the lm dead grass
    sandL.lm.slope<-round(coef(lm(distL.sim.df$sim[which(distL.sim.df$G_S=="S")]~distL.sim.df$dist[which(distL.sim.df$G_S=="S")]))[2],6)#slope of the lm dead sand
    
    ggplot(data = distL.sim.df, aes(x=dist,y=sim)) +
      #geom_point(col="red") +
      geom_point(aes(color=factor(G_S))) +
      geom_smooth(method="lm", aes(color=factor(G_S)), fill = "gray") + #doesn't work yet
      ylim(0,1) +
      geom_label(aes(x=22000,y=0.95),label=paste("slope vegetated =",grassL.lm.slope),show.legend=F) +
      geom_label(aes(x=22000,y=0.85),label=paste("slope unvegetated =",sandL.lm.slope),show.legend=F) +
      scale_color_manual(values=c("G" = "red",
                                  "S" = "blue"))
    
#Decomposition of Beta diversity using Baselga's Bray-Curtis indices####
    #NOTE: cannot use percent transformed data!!!!!! - currently using log (can use raw by commenting out each of the decostand lines)
#dead
  #combined
    
  b1<-decostand(Cay.raw.nosingles,"log") #corresponds to alpha in ADP.groups
    
  b2<-as.data.frame(rbind(#b2 level - corresponds to b1 in ADP.groups
    RPG1=colSums(Cay.raw.nosingles[RPG[1:4],]),Cay.raw.nosingles[RPG[5],],Cay.raw.nosingles[RPG[6],],Cay.raw.nosingles[RPG[7],],
    RPS1=colSums(Cay.raw.nosingles[RPS[1:4],]),Cay.raw.nosingles[RPS[5],],Cay.raw.nosingles[RPS[6],],Cay.raw.nosingles[RPS[7],],
    SSG1=colSums(Cay.raw.nosingles[SSG[1:4],]),Cay.raw.nosingles[SSG[5],],Cay.raw.nosingles[SSG[6],],Cay.raw.nosingles[SSG[7],],
    SSS1=colSums(Cay.raw.nosingles[SSS[1:4],]),Cay.raw.nosingles[SSS[5],],Cay.raw.nosingles[SSS[6],],Cay.raw.nosingles[SSS[7],],
    EBG1=colSums(Cay.raw.nosingles[EBG[1:4],]),Cay.raw.nosingles[EBG[5],],Cay.raw.nosingles[EBG[6],],Cay.raw.nosingles[EBG[7],],
    EBS1=colSums(Cay.raw.nosingles[EBS[1:4],]),Cay.raw.nosingles[EBS[5],],Cay.raw.nosingles[EBS[6],],Cay.raw.nosingles[EBS[7],]
    ))
    b2<-decostand(b2,"log")
    
  b3<-as.data.frame(rbind(#b3 level - corresponds to b2 in ADP.groups
    RPG=colSums(Cay.raw.nosingles[RPG,]),
    RPS=colSums(Cay.raw.nosingles[RPS,]),
    SSG=colSums(Cay.raw.nosingles[SSG,]),
    SSS=colSums(Cay.raw.nosingles[SSS,]),
    EBG=colSums(Cay.raw.nosingles[EBG,]),
    EBS=colSums(Cay.raw.nosingles[EBS,])
    ))
  b3<-decostand(b3,"log")

  b4<-as.data.frame(rbind(#b4 level - corresponds to b3 in ADP.groups
    RP=colSums(Cay.raw.nosingles[RP,]),
    SS=colSums(Cay.raw.nosingles[SS,]),
    EB=colSums(Cay.raw.nosingles[EB,])
    ))
  b4<-decostand(b4,"log")

  #bray.b1.s<-beta.multi.abund(Cay.raw.nosingles,"bray")
  #rc.b1.s<-mean(raupcrick(Cay.raw.nosingles))#raup-crick indexbray.b2.s<-beta.multi.abund(b2)
  bray.b1.s<-beta.multi.abund(b1,"bray")
  rc.b1.s<-mean(raupcrick(b1))
  bray.b2.s<-beta.multi.abund(b2)
  rc.b2.s<-mean(raupcrick(b2))
  bray.b3.s<-beta.multi.abund(b3)
  rc.b3.s<-mean(raupcrick(b3))
  bray.b4.s<-beta.multi.abund(b4)
  rc.b4.s<-mean(raupcrick(b4))

  #grass
  #bray.b1.s.grass<-beta.multi.abund(Cay.raw.nosingles[Grass,],"bray")#only grass from Cay.raw.nosingles
  #rc.b1.s.grass<-mean(raupcrick(Cay.raw.nosingles[Grass,]))
  bray.b1.s.grass<-beta.multi.abund(b1[Grass,],"bray")#only grass from Cay.raw.nosingles
  rc.b1.s.grass<-mean(raupcrick(b1[Grass,]))
  bray.b2.s.grass<-beta.multi.abund(b2[c(1:4,9:12,17:20),])#only grass from b2
  rc.b2.s.grass<-mean(raupcrick(b2[c(1:4,9:12,17:20),]))
  bray.b3.s.grass<-beta.multi.abund(b3[c(1,3,5),])#only grass from b3
  rc.b3.s.grass<-mean(raupcrick(b3[c(1,3,5),]))
  
  #sand
  #bray.b1.s.sand<-beta.multi.abund(Cay.raw.nosingles[Sand,],"bray")#only sand from Cay.raw.nosingles
  #rc.b1.s.sand<-mean(raupcrick(Cay.raw.nosingles[Sand,]))
  bray.b1.s.sand<-beta.multi.abund(b1[Sand,],"bray")#only sand from Cay.raw.nosingles
  rc.b1.s.sand<-mean(raupcrick(b1[Sand,]))
  bray.b2.s.sand<-beta.multi.abund(b2[c(5:8,13:16,21:24),])#only sand from b2
  rc.b2.s.sand<-mean(raupcrick(Cay.raw.nosingles[c(5:8,13:16,21:24),]))
  bray.b3.s.sand<-beta.multi.abund(b3[c(2,4,6),])#only sand from b3
  rc.b3.s.sand<-mean(raupcrick(Cay.raw.nosingles[c(2,4,6),]))

  bray.df<-data.frame(rbind(B1=bray.b1.s,B2=bray.b2.s,B3=bray.b3.s,B4=bray.b4.s,
                B1.grass=bray.b1.s.grass,B2.grass=bray.b2.s.grass,B3.grass=bray.b3.s.grass,
                B1.sand=bray.b1.s.sand,B2.sand=bray.b2.s.sand,B3.sand=bray.b3.s.sand
                ))
  dissim.df<-data.frame(bray.df,raup_crick=c(rc.b1.s,rc.b2.s,rc.b3.s,rc.b4.s,rc.b1.s.grass,rc.b2.s.grass,rc.b3.s.grass,rc.b1.s.sand,rc.b2.s.sand,rc.b3.s.sand))

#Live
  #combined
  b1L<-decostand(CayL.raw.int,"log")
  b2L<-as.data.frame(rbind(#b2 level - corresponds to ADP.groups
    RPG1=colSums(CayL.raw.int[RPG[1:4],]),CayL.raw.int[RPG[5],],CayL.raw.int[RPG[6],],CayL.raw.int[RPG[7],],
    RPS1=colSums(CayL.raw.int[RPS[1:4],]),CayL.raw.int[RPS[5],],CayL.raw.int[RPS[6],],CayL.raw.int[RPS[7],],
    SSG1=colSums(CayL.raw.int[SSG[1:4],]),CayL.raw.int[SSG[5],],CayL.raw.int[SSG[6],],CayL.raw.int[SSG[7],],
    SSS1=colSums(CayL.raw.int[SSS[1:4],]),CayL.raw.int[SSS[5],],CayL.raw.int[SSS[6],],CayL.raw.int[SSS[7],],
    EBG1=colSums(CayL.raw.int[EBG[1:4],]),CayL.raw.int[EBG[5],],CayL.raw.int[EBG[6],],CayL.raw.int[EBG[7],],
    EBS1=colSums(CayL.raw.int[EBS[1:4],]),CayL.raw.int[EBS[5],],CayL.raw.int[EBS[6],],CayL.raw.int[EBS[7],]
    ))
  b2L<-decostand(b2L,"log")
  
  b3L<-as.data.frame(rbind(#b3 level - corresponds to ADP.groups
    RPG=colSums(CayL.raw.int[RPG,]),
    RPS=colSums(CayL.raw.int[RPS,]),
    SSG=colSums(CayL.raw.int[SSG,]),
    SSS=colSums(CayL.raw.int[SSS,]),
    EBG=colSums(CayL.raw.int[EBG,]),
    EBS=colSums(CayL.raw.int[EBS,])
    ))
  b3L<-decostand(b3L,"log")
  
  b4L<-as.data.frame(rbind(#b4 level - corresponds to ADP.groups
    RP=colSums(CayL.raw.int[RP,]),
    SS=colSums(CayL.raw.int[SS,]),
    EB=colSums(CayL.raw.int[EB,])
    ))
  b4L<-decostand(b3L,"log")

  #brayL.b1.s<-beta.multi.abund(CayL.raw.int,"bray")
  #rcL.b1.s<-mean(raupcrick(CayL.raw.int))#raup-crick indexbray.b2.s<-beta.multi.abund(b2)
  brayL.b1.s<-beta.multi.abund(b1,"bray")
  rcL.b1.s<-mean(raupcrick(b1))#raup-crick indexbray.b2.s<-beta.multi.abund(b2)
  brayL.b2.s<-beta.multi.abund(b2L)
  rcL.b2.s<-mean(raupcrick(b2L))
  brayL.b3.s<-beta.multi.abund(b3L)
  rcL.b3.s<-mean(raupcrick(b3L))
  brayL.b4.s<-beta.multi.abund(b4L)
  rcL.b4.s<-mean(raupcrick(b4L))
  
  #grass
  #brayL.b1.s.grass<-beta.multi.abund(CayL.raw.int[Grass,],"bray")#only grass from CayL.raw.int
  #rcL.b1.s.grass<-mean(raupcrick(CayL.raw.int[Grass,]))
  brayL.b1.s.grass<-beta.multi.abund(b1[Grass,],"bray")#only grass from CayL.raw.int
  rcL.b1.s.grass<-mean(raupcrick(b1[Grass,]))
  brayL.b2.s.grass<-beta.multi.abund(b2L[c(1:4,9:12,17:20),])#only grass from b2
  rcL.b2.s.grass<-mean(raupcrick(b2L[c(1:4,9:12,17:20),]))
  brayL.b3.s.grass<-beta.multi.abund(b3L[c(1,3,5),])#only grass from b3
  rcL.b3.s.grass<-mean(raupcrick(b3L[c(1,3,5),]))
  
  #sand
  #brayL.b1.s.sand<-beta.multi.abund(CayL.raw.int[Sand,],"bray")#only sand from CayL.raw.int
  #rcL.b1.s.sand<-mean(raupcrick(CayL.raw.int[Sand,]))
  brayL.b1.s.sand<-beta.multi.abund(b1[Sand,],"bray")#only sand from CayL.raw.int
  rcL.b1.s.sand<-mean(raupcrick(b1[Sand,]))
  brayL.b2.s.sand<-beta.multi.abund(b2L[c(5:8,13:16,21:24),])#only sand from b2
  rcL.b2.s.sand<-mean(raupcrick(b2L[c(5:8,13:16,21:24),]))
  brayL.b3.s.sand<-beta.multi.abund(b3L[c(2,4,6),])#only sand from b3
  rcL.b3.s.sand<-mean(raupcrick(CayL.raw.int[c(2,4,6),]))
  
  brayL.df<-data.frame(rbind(B1L=brayL.b1.s,B2L=brayL.b2.s,B3L=brayL.b3.s,B4L=brayL.b4.s,
                            B1L.grass=brayL.b1.s.grass,B2L.grass=brayL.b2.s.grass,B3L.grass=brayL.b3.s.grass,
                            B1L.sand=brayL.b1.s.sand,B2L.sand=brayL.b2.s.sand,B3L.sand=brayL.b3.s.sand
                            ))
  dissimL.df<-data.frame(brayL.df,raup_crick=c(rcL.b1.s,rcL.b2.s,rcL.b3.s,rcL.b4.s,rcL.b1.s.grass,rcL.b2.s.grass,rcL.b3.s.grass,rcL.b1.s.sand,rcL.b2.s.sand,rcL.b3.s.sand))
  
dist.components<-as.data.frame(rbind(dissim.df,dissimL.df))
dist.comp<-data.frame(spatial_scale=c(rep(c("B1","B2","B3","B4",rep(c("B1","B2","B3/B4"),2)),2)),subset=rep(c(rep("all",4),rep("vegetated",3),rep("unvegetated",3)),2),live_dead=c(rep("dead",10),rep("live",10)),B.BC=as.numeric(dist.components$beta.BRAY),B.bal=as.numeric(dist.components$beta.BRAY.BAL),B.gra=as.numeric(dist.components$beta.BRAY.GRA),percent_Bbal=c(as.numeric(dist.components$beta.BRAY.BAL)/as.numeric(dist.components$beta.BRAY)*100),B.RC=as.numeric(dist.components$raup_crick))
dist.components.mat<-as.matrix(dist.comp)
write.csv(dist.components.mat,"Bray_Dissimilarity_Table.csv")  #TABLE 3




#EXTRA STUFF####

#raw data
bray.b1<-beta.multi.abund(Cay.raw.nosingles,"bray")
bray.b2<-beta.multi.abund(as.data.frame(rbind(
  RPG1=colSums(Cay.raw.nosingles[RPG[1:4],]),Cay.raw.nosingles[RPG[5],],Cay.raw.nosingles[RPG[6],],Cay.raw.nosingles[RPG[7],],
  RPS1=colSums(Cay.raw.nosingles[RPS[1:4],]),Cay.raw.nosingles[RPS[5],],Cay.raw.nosingles[RPS[6],],Cay.raw.nosingles[RPS[7],],
  SSG1=colSums(Cay.raw.nosingles[SSG[1:4],]),Cay.raw.nosingles[SSG[5],],Cay.raw.nosingles[SSG[6],],Cay.raw.nosingles[SSG[7],],
  SSS1=colSums(Cay.raw.nosingles[SSS[1:4],]),Cay.raw.nosingles[SSS[5],],Cay.raw.nosingles[SSS[6],],Cay.raw.nosingles[SSS[7],],
  EBG1=colSums(Cay.raw.nosingles[EBG[1:4],]),Cay.raw.nosingles[EBG[5],],Cay.raw.nosingles[EBG[6],],Cay.raw.nosingles[EBG[7],],
  EBS1=colSums(Cay.raw.nosingles[EBS[1:4],]),Cay.raw.nosingles[EBS[5],],Cay.raw.nosingles[EBS[6],],Cay.raw.nosingles[EBS[7],])))
bray.b3<-beta.multi.abund(as.data.frame(rbind(
  RPG=colSums(Cay.raw.nosingles[RPG,]),
  RPS=colSums(Cay.raw.nosingles[RPS,]),
  SSG=colSums(Cay.raw.nosingles[SSG,]),
  SSS=colSums(Cay.raw.nosingles[SSS,]),
  EBG=colSums(Cay.raw.nosingles[EBG,]),
  EBS=colSums(Cay.raw.nosingles[EBS,]))))
bray.b4<-beta.multi.abund(as.data.frame(rbind(
  RP=colSums(Cay.raw.nosingles[RP,]),
  SS=colSums(Cay.raw.nosingles[SS,]),
  EB=colSums(Cay.raw.nosingles[EB,]))))

###### With one less Beta comparison group (B2):
ADP.groups<-read.csv(file="ADP_groups.csv",header=T,row.names=1)
Cay.ADP.S<-adipart(Cay.raw.nosingles ~ .,ADP.groups,index="richness",nsimul=999,method="r2dtable",na.rm=T)#ADP Richness
#below, Hill=1 "true diversity"
Cay.ADP.S<-adipart(Cay.raw.nosingles ~ .,ADP.groups,index="Hill1",nsimul=999,method="r2dtable",na.rm=T)#ADP Richness
Cay.ADP.S.weight<-adipart(Cay.raw.nosingles ~ .,ADP.groups,index="richness",nsimul=999,method="r2dtable",na.rm=T,weight="prop")#ADP Richness
#split by grass and sand
Cay.ADP.S.grass<-adipart(Cay.raw.nosingles[c(RPG,SSG,EBG),] ~ .,ADP.groups[c(RPG,SSG,EBG),],index="richness",nsimul=999,method="r2dtable",na.rm=T)
Cay.ADP.S.sand<-adipart(Cay.raw.nosingles[c(RPS,SSS,EBS),] ~ .,ADP.groups[c(RPS,SSS,EBS),],index="richness",nsimul=999,method="r2dtable",na.rm=T)
#weighted by sample abundance
Cay.ADP.S.grass.weight<-adipart(Cay.raw.nosingles[c(RPG,SSG,EBG),] ~ .,ADP.groups[c(RPG,SSG,EBG),],index="richness",nsimul=999,method="r2dtable",na.rm=T,weight="prop")
Cay.ADP.S.sand.weight<-adipart(Cay.raw.nosingles[c(RPS,SSS,EBS),] ~ .,ADP.groups[c(RPS,SSS,EBS),],index="richness",nsimul=999,method="r2dtable",na.rm=T,weight="prop")

ADP.Table<-data.frame(sites=c(rep("all",12),rep("grass",12),rep("sand",12)),#36
                      scale=rep(c("a1","B1","B2","B3","B4","g"),6),#36
                      weight=rep(c(rep("uniform",6),rep("weighted",6)),3),#36
                      observed=round(c(Cay.ADP.S$statistic[c(1,6:9,5)],Cay.ADP.S.weight$statistic[c(1,6:9,5)],Cay.ADP.S.grass$statistic[c(1,6:9,5)],Cay.ADP.S.grass.weight$statistic[c(1,6:9,5)],Cay.ADP.S.sand$statistic[c(1,6:9,5)],Cay.ADP.S.sand.weight$statistic[c(1,6:9,5)])),
                      percent=c(round(Cay.ADP.S$statistic[c(1,6:9,5)]/Cay.ADP.S$statistic[5]*100,2),round(Cay.ADP.S.weight$statistic[c(1,6:9,5)]/Cay.ADP.S.weight$statistic[5]*100,2),round(Cay.ADP.S.grass$statistic[c(1,6:9,5)]/Cay.ADP.S.grass$statistic[5]*100,2),round(Cay.ADP.S.grass.weight$statistic[c(1,6:9,5)]/Cay.ADP.S.grass.weight$statistic[5]*100,2),round(Cay.ADP.S.sand$statistic[c(1,6:9,5)]/Cay.ADP.S.sand$statistic[5]*100,2),round(Cay.ADP.S.sand.weight$statistic[c(1,6:9,5)]/Cay.ADP.S.sand.weight$statistic[5]*100,2)),#36
                      expected=round(c(Cay.ADP.S$oecosimu$means[c(1,6:9,5)],Cay.ADP.S.weight$oecosimu$means[c(1,6:9,5)],Cay.ADP.S.grass$oecosimu$means[c(1,6:9,5)],Cay.ADP.S.grass.weight$oecosimu$means[c(1,6:9,5)],Cay.ADP.S.sand$oecosimu$means[c(1,6:9,5)],Cay.ADP.S.sand.weight$oecosimu$means[c(1,6:9,5)])),#36
                      CI2.5=round(c(apply(Cay.ADP.S$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],apply(Cay.ADP.S.weight$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],apply(Cay.ADP.S.grass$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],apply(Cay.ADP.S.grass.weight$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],apply(Cay.ADP.S.sand$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)],apply(Cay.ADP.S.sand.weight$oecosimu$simulated, 1, quantile, probs=c(0.025))[c(1,6:9,5)])),#36
                      CI97.5=round(c(apply(Cay.ADP.S$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],apply(Cay.ADP.S.weight$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],apply(Cay.ADP.S.grass$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],apply(Cay.ADP.S.grass.weight$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],apply(Cay.ADP.S.sand$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)],apply(Cay.ADP.S.sand.weight$oecosimu$simulated, 1, quantile, probs=c(0.975))[c(1,6:9,5)])),#36
                      pval=c(Cay.ADP.S$oecosimu$pval[c(1,6:9,NA)],Cay.ADP.S.weight$oecosimu$pval[c(1,6:9,NA)],Cay.ADP.S.grass$oecosimu$pval[c(1,6:9,NA)],Cay.ADP.S.grass.weight$oecosimu$pval[c(1,6:9,NA)],Cay.ADP.S.sand$oecosimu$pval[c(1,6:9,NA)],Cay.ADP.S.sand.weight$oecosimu$pval[c(1,6:9,NA)])#36
                      )

#Master ADP table in long form for faceting in ggplot only
ADP.Table.long<-gather(ADP.Table,key="measure",value="value",c("observed","expected"))
#all in one plot, above

#master ADP dataframe
Cay.ADPdf<-data.frame(measure=c(rep("Uniform",10),rep("Weighted",10)),level=rep(c("alpha.1","beta.1","beta.2","beta.3","beta.4"),4),type=rep(c(rep("observed",5),rep("expected",5)),2),value=c(Cay.ADP.S$statistic[c(1,6:9)],Cay.ADP.S$oecosimu$means[c(1,6:9)],Cay.ADP.S.weight$statistic[c(1,6:9)],Cay.ADP.S.weight$oecosimu$means[c(1,6:9)]))
#master ADP dataframes for grass and sand separately
Cay.ADPdf.grass<-data.frame(measure=c(rep("Uniform",10),rep("Weighted",10)),level=rep(c("alpha.1","beta.1","beta.2","beta.3","beta.4"),4),type=rep(c(rep("observed",5),rep("expected",5)),2),value=round(c(Cay.ADP.S.grass$statistic[c(1,6:9)],Cay.ADP.S.grass$oecosimu$means[c(1,6:9)],Cay.ADP.S.grass.weight$statistic[c(1,6:9)],Cay.ADP.S.grass.weight$oecosimu$means[c(1,6:9)]),1))#statistic values are rounded to avoid negative numbers
Cay.ADPdf.sand<-data.frame(measure=c(rep("Uniform",10),rep("Weighted",10)),level=rep(c("alpha.1","beta.1","beta.2","beta.3","beta.4"),4),type=rep(c(rep("observed",5),rep("expected",5)),2),value=round(c(Cay.ADP.S.sand$statistic[c(1,6:9)],Cay.ADP.S.sand$oecosimu$means[c(1,6:9)],Cay.ADP.S.sand.weight$statistic[c(1,6:9)],Cay.ADP.S.sand.weight$oecosimu$means[c(1,6:9)]),1))#statistic values are rounded to avoid negative numbers

quartz()
ggplot(Cay.ADPdf, 
       aes(fill=factor(level,levels=c("beta.4","beta.3","beta.2","beta.1","alpha.1")),x=factor(type,levels=c("observed","expected")),y=value)) + 
  geom_bar(position="fill", stat="identity") +
  labs(fill= "partition level",x="",title="Cay ADP") +
  facet_grid(~ factor(measure,levels=c("Richness","Shannon")))#used factor to reorder localities in facets

quartz()
ggplot(Cay.ADPdf.weight, 
       aes(fill=factor(level,levels=c("beta.4","beta.3","beta.2","beta.1","alpha.1")),x=factor(type,levels=c("observed","expected")),y=value)) + 
  geom_bar(position="fill", stat="identity") +
  labs(fill= "partition level",x="",title="Cay ADP - Weighted by abundance") +
  facet_grid(~ factor(measure,levels=c("Richness","Shannon")))#used factor to reorder localities in facets

quartz()
ggplot(Cay.ADPdf.grass, 
       aes(fill=factor(level,levels=c("beta.4","beta.3","beta.2","beta.1","alpha.1")),x=factor(type,levels=c("observed","expected")),y=value)) + 
  geom_bar(position="fill", stat="identity") +
  labs(fill= "partition level",x="",title="Cay ADP - Grass") +
  facet_grid(~ factor(measure,levels=c("Richness","Shannon")))#used factor to reorder localities in facets

quartz()
ggplot(Cay.ADPdf.grass.weight, 
       aes(fill=factor(level,levels=c("beta.4","beta.3","beta.2","beta.1","alpha.1")),x=factor(type,levels=c("observed","expected")),y=value)) + 
  geom_bar(position="fill", stat="identity") +
  labs(fill= "partition level",x="",title="Cay ADP - Grass Weighted") +
  facet_grid(~ factor(measure,levels=c("Richness","Shannon")))#used factor to reorder localities in facets

quartz()
ggplot(Cay.ADPdf.sand, 
       aes(fill=factor(level,levels=c("beta.4","beta.3","beta.2","beta.1","alpha.1")),x=factor(type,levels=c("observed","expected")),y=value)) + 
  geom_bar(position="fill", stat="identity") +
  labs(fill= "partition level",x="",title="Cay ADP - Sand") +
  facet_grid(~ factor(measure,levels=c("Richness","Shannon")))#used factor to reorder localities in facets

quartz()
ggplot(Cay.ADPdf.sand.weight, 
       aes(fill=factor(level,levels=c("beta.4","beta.3","beta.2","beta.1","alpha.1")),x=factor(type,levels=c("observed","expected")),y=value)) + 
  geom_bar(position="fill", stat="identity") +
  labs(fill= "partition level",x="",title="Cay ADP - Sand Weighted") +
  facet_grid(~ factor(measure,levels=c("Richness","Shannon")))#used factor to reorder localities in facets

#calculating 2.5 and 97.5 confidence intervals for ADP analysis
quants.Cay.ADP<-apply(Cay.ADP.S$oecosimu$simulated, 1, quantile, probs=c(0.025, 0.975))
quants.Cay.ADP.grass<-apply(Cay.ADP.S.grass$oecosimu$simulated, 1, quantile, probs=c(0.025, 0.975))
quants.Cay.ADP.sand<-apply(Cay.ADP.S.sand$oecosimu$simulated, 1, quantile, probs=c(0.025, 0.975))

#log plot, but x labels don't work right... need to replace log t values with distances...
plot(log(Cay.decay.exp$data.x),Cay.decay.exp$data.y)

switch(Cay.decay.exp$y.type, similarities = {
  lines(fitted(Cay.decay.exp$model)[order(log(Cay.decay.exp$data.x))] ~ sort(log(Cay.decay.exp$data.x)), 
  )
}, dissimilarities = {
  lines(1 - fitted(Cay.decay.exp$model)[order(log(Cay.decay.exp$data.x))] ~ sort(log(Cay.decay.exp$data.x)), 
  )
})

#calculating distance-decay models####
#dead
Cay.decay.exp<-decay.model(Cay.sim,Cay.spat.dist,y.type="similarities",model.type="exp",perm=100)#exponential
Cay.decay.pow<-decay.model(Cay.dist,Cay.spat.dist,y.type="dissim",model.type="pow",perm=100)#power
#live
CayL.decay.exp<-decay.model(CayL.sim,Cay.spat.dist,y.type="similarities",model.type="exp",perm=100)#exponential

#plotting decay models
#par(mfrow=c(1,2))
#dead
#plot.decay(Cay.decay.exp,add=F,remove.dots=F,lwd=1,ylim=c(0.1,0.8),xlim=c(0,50))
plot.decay(Cay.decay.exp,add=F,remove.dots=F,lwd=1,ylim=c(0.1,0.8))#linear y axis
plot.decay(Cay.decay.exp,add=F,remove.dots=F,lwd=1,ylim=c(0.19,1),log="y")#log y axis
#live
plot.decay(CayL.decay.exp,add=F,remove.dots=F,lwd=1,col="blue")#linear y axis

#from: https://www.researchgate.net/post/How-can-I-estimate-the-distance-decay-of-similarity-in-communities
#dead
dis.dec<-lm(log(Cay.sim) ~ Cay.spat.dist)
#dis.dec<-lm(Cay.sim ~ Cay.spat.dist)
plot(Cay.spat.dist,log(Cay.sim))
abline(dis.dec)
dead.mantel<-mantel(Cay.sim,Cay.spat.dist,perm=999)
#dead, grass and sand separate
#dis.decG<-lm(Cay.G.sim ~ Cay.G.spat.dist)#all distances
#dis.decG.less25K<-lm(Cay.G.sim[Cay.spat.dist<25000] ~ Cay.G.spat.dist[Cay.spat.dist<25000])
dis.decS<-lm(Cay.S.sim ~ Cay.S.spat.dist)

plot(Cay.G.spat.dist,Cay.G.sim,col="red")#all distances
abline(lm(Cay.G.sim ~ Cay.G.spat.dist),col="red")
points(Cay.S.spat.dist,Cay.S.sim,col="blue")
abline(lm(Cay.S.sim ~ Cay.S.spat.dist),col="blue")
abline(lm(Cay.G.sim[Cay.spat.dist<25000] ~ Cay.G.spat.dist[Cay.spat.dist<25000]),col="red")#trying to eliminate >25,000 m distances; doesn't work
abline(lm(Cay.S.sim[Cay.spat.dist<25000] ~ Cay.S.spat.dist[Cay.spat.dist<25000]),col="blue")#trying to eliminate >25,000 m distances; doesn't work

G.mantel<-mantel(Cay.G.sim,Cay.G.spat.dist,perm=9999)#like cor test; what is significance value?
#cor.test(Cay.G.sim,Cay.G.spat.dist)
S.mantel<-mantel(Cay.G.sim,Cay.G.spat.dist,perm=999)

dis.decL<-lm(log(CayL.sim) ~ Cay.spat.dist)#log-similarity, doesn't work
dis.decL<-lm(CayL.sim ~ Cay.spat.dist)#non log-similarity, works
plot(Cay.spat.dist,log(CayL.sim))
abline(log(dis.decL))
error<-mantel(Cay.sim,Cay.spat.dist,perm=999)

#Basic Diversity Statistics (bivalves and gastropods separately)
#dead
sum(Cay.raw.preMNI)#total elements
sum(Cay.raw.preMNI[,68:ncol(Cay.raw.preMNI)])#total gastropods
sum(Cay.raw.preMNI[,1:67])#total bivalves
biv.fam<-unique(Cay.fam$Family[1:67])#total bivalve families
gast.fam<-unique(Cay.fam$Family[68:168])#total bivalve families
#bivalves
S.biv<-specnumber(Cay.raw.site[,1:67])
mean(S.biv)#mean richness among bivalves
sd(specnumber(Cay.raw.site[,1:67]))#st deviation of richness among bivalves
mean(rowSums(Cay.raw.site[,1:67]))#total individuals among bivalves
sd(rowSums(Cay.raw.site[,1:67]))#total individuals among bivalves
#gastropods
S.gast<-specnumber(Cay.raw.site[,68:168])
mean(S.gast)#mean richness among gastropods
sd(specnumber(Cay.raw.site[,68:168]))#st deviation of richness among gastropods
mean(rowSums(Cay.raw.site[,68:168]))#total individuals among gastropods
sd(rowSums(Cay.raw.site[,68:168]))#total individuals among gastropods
#bivalve evenness
diversity(Cay.raw.site[,1:67])/log(S.biv)
sd(diversity(Cay.raw.site[,1:67])/log(S.biv))#standard deviation of bivalve evenness
#gastropod evenness
diversity(Cay.raw.site[,68:168])/log(S.gast)
sd(diversity(Cay.raw.site[,68:168])/log(S.gast))
#live
sum(CayL.raw.preMNI[,68:ncol(Cay.raw.preMNI)])#total gastropods
sum(CayL.raw.preMNI[,1:67])#total bivalves
bivL.fam<-unique(Cay.fam$Family[1:67])#total bivalve families
gastL.fam<-unique(Cay.fam$Family[68:168])#total bivalve families
#bivalves
S.bivL<-specnumber(CayL.raw.site[,1:67])
mean(S.bivL)#mean richness among bivalves
sd(specnumber(CayL.raw.site[,1:67]))#st deviation of richness among bivalves
mean(rowSums(CayL.raw.site[,1:67]))#total individuals among bivalves
sd(rowSums(CayL.raw.site[,1:67]))#total individuals among bivalves
#gastropods
S.gastL<-specnumber(CayL.raw.site[,68:168])
mean(S.gastL)#mean richness among gastropods
sd(specnumber(CayL.raw.site[,68:168]))#st deviation of richness among gastropods
mean(rowSums(CayL.raw.site[,68:168]))#total individuals among gastropods
sd(rowSums(CayL.raw.site[,68:168]))#total individuals among gastropods
#bivalve evenness
diversity(CayL.raw.site[,1:67])/log(S.bivL)
sd(diversity(CayL.raw.site[,1:67])/log(S.bivL))#standard deviation of bivalve evenness
#gastropod evenness
diversity(CayL.raw.site[,68:168])/log(S.gastL)
sd(diversity(CayL.raw.site[,68:168])/log(S.gastL))

#Boxplots of Richness
Richness.D<-data.frame(site=rep(c(rep("RPG",7),rep("RPS",7),rep("SSG",7),rep("SSS",7),rep("EBG",7),rep("EBS",7)),3),subdiv=c(rep("all",42),rep("bivalve",42),rep("gastropod",42)),S=c(specnumber(Cay.raw.site),specnumber(Cay.raw.site[,c(1:67)]),specnumber(Cay.raw.site[,c(68:168)])))
Richness.D$site<-factor(Richness.D$site,levels=c("RPG","RPS","SSG","SSS","EBG","EBS"))
Richness.L<-data.frame(site=rep(c(rep("RPG",7),rep("RPS",7),rep("SSG",7),rep("SSS",7),rep("EBG",7),rep("EBS",7)),3),subdiv=c(rep("all",42),rep("bivalve",42),rep("gastropod",42)),S=c(specnumber(CayL.raw.site),specnumber(CayL.raw.site[,c(1:67)]),specnumber(CayL.raw.site[,c(68:168)])))
Richness.L$site<-factor(Richness.L$site,levels=c("RPG","RPS","SSG","SSS","EBG","EBS"))
Richness.df<-data.frame(DL=c(rep("Dead",126),rep("Live",126)),rbind(Richness.D,Richness.L))
#Richness.df$site<-factor(Richness.df$site,levels=c("RPG","RPS","SSG","SSS","EBG","EBS"))
#Boxplot, Live and Dead
ggplot(Richness.df, aes(x=site, y=S)) + 
  geom_boxplot() +
  facet_wrap(~ DL + subdiv,nrow=2,ncol=3)#facet
#Boxplot, Live-only
ggplot(Richness.L, aes(x=site, y=S)) + 
  geom_boxplot() +
  facet_wrap(~ subdiv,nrow=1,ncol=3)#facet
