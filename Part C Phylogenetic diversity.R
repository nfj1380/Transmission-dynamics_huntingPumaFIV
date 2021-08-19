#-----------------------------------------------#
#####Part C: Code to perform phylodynamic and phylogentic diversity analyses#####
#by Nick Fountain-Jones 
#-----------------------------------------------#


#------------------------------------------------------------------------
##################Reload MCC Tree##################
#------------------------------------------------------------------------

beastTree<- read.beast('CO_FIV_Gamma_LogC_ExpGrowth_Cauchy3_Pop3.tre') #load the tree

#plot the tree

phy_simple<- ggtree(beastTree) + geom_tiplab(size=2)+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 

phy_simple

#prune each dominant clade in each region. Based on FOuntain-Jones etal 2020

Gamma_treatment <- tree_subset(beastTree , node=66,levels_back = 0)
Gamma_stable <- tree_subset(beastTree , node=43,levels_back = 0)


#------------------------------------------------------------------------
##################Look at NE through time for various lineages (Karcher et al)##################
#------------------------------------------------------------------------
library(phylodyn)

#ps helps account fro preferential sampling

BSpsFR3 <- BNPR_PS(bFFV, lengthout = 100, prec_alpha = 0.01, prec_beta = 0.01,
                   beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                   derivative = FALSE, forward = TRUE)

plot_BNPR(BSpsFR3)

BSpsWS1 <- BNPR_PS(b1, lengthout = 100, prec_alpha = 0.01, prec_beta = 0.01,
                   beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                   derivative = FALSE, forward = TRUE)

plot_BNPR(BSpsWS1)

BSpsWS2 <- BNPR_PS(b2, lengthout = 100, prec_alpha = 0.01, prec_beta = 0.01,
                   beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                   derivative = FALSE, forward = TRUE)
BPFFVWS1 <- BNPR_PS(b2, lengthout = 100, prec_alpha = 0.01, prec_beta = 0.01,
                    beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                    derivative = FALSE, forward = TRUE)

plot_BNPR(BSpsWS2)

#------------------------------------------------------------------------
##################Phylogentic diversity analyses##################
#------------------------------------------------------------------------


#------------------------------------------------------------------------
##################Classify sequences into each contrast##################
#------------------------------------------------------------------------

comm <- as.data.frame(get_taxa_name(tree_view = NULL, node = NULL))

names(comm) <- c('date') 
comm <- comm %>% separate(date,c("ID","Site","Year", "month", "day"))

#for FFV
noHunt <-comm %>% filter(Year < 2011)
Hunt <-comm %>% filter(Year > 2011)

noHunt <- add_column(noHunt, noHunt = 1, Hunt=0)
Hunt <- add_column(Hunt, noHunt = 0, Hunt=1)
Combined <- bind_rows(noHunt, Hunt,id = NULL)
commDate <- Combined %>% unite(date,Year, month, day, sep='-')
comm <- commDate %>% unite(seq,ID, Site, date )
rownames(comm) <- comm$seq
comm$seq <- NULL
#needs to be transmposed to work
comm <-t(as.matrix(comm))

#forFIV

noHuntWS <-comm %>% filter(Year <= 2011, Site=='WS')
HuntWS <-comm %>% filter(Year > 2011, Site=='WS')

noHuntWS <- noHuntWS[-c(10),] #resampled individual should be removed

noHuntFR <-comm %>% filter(Year <= 2011, Site=='FR')
HuntFR <-comm %>% filter(Year > 2011, Site=='FR')

noHuntFR <- add_column(noHuntFR, noHuntFR = 1, HuntFR=0, noHuntWS=0,HuntWS=0 )
HuntFR <- add_column(HuntFR, noHuntFR = 0, HuntFR=1, noHuntWS=0,HuntWS=0 )
HuntWS <- add_column(HuntWS, noHuntFR = 0, HuntFR=0, noHuntWS=0,HuntWS=1 )
noHuntWS <- add_column(noHuntWS, noHuntFR = 0, HuntFR=0, noHuntWS=1,HuntWS=0 )

Combined <- bind_rows(noHuntWS, HuntWS,noHuntFR, HuntFR, id = NULL)
commDate <- Combined %>% unite(date,Year, month, day, sep='-')
comm <- commDate %>% unite(seq,ID, Site, date )
rownames(comm) <- comm$seq
comm$seq <- NULL

#----------------------------
#just get WS phylogeny
#----------------------------

WSnames <- Combined <- bind_rows(noHuntWS, HuntWS, id = NULL)
WSnamesDate <- WSnames  %>% unite(date,Year, month, day, sep='-')
WStips <- WSnamesDate %>% unite(seq,ID, Site, date )

WStipsT<- t(WStips)
WStipsT <- as.data.frame(WStipsT, head=T)
WStipsT <-WStipsT[-1,] 
WStipsT <-as.matrix(WStipsT) 

mergeData <- comparative.comm(b, WStipsT)
pruned.tree <- mergeData$phy
plot(pruned.tree)
str(pruned.tree)
#write the tree
write.tree(pruned.tree, file= 'WSFIVGammaTree.newick')

#comm <- select(comm,-c(seq)) #maybe a conflict here with DAMOCLES
comm <-comm[,-1]
comm <-t(as.matrix(comm))

#plot on the phylogenies  -this needs work 

par(mfrow = c(2, 2)) #phylogeney has to be a phylo object so use 'b' above
for (i in row.names(comm)) {
  plot(b, show.tip.label = FALSE, main = i)
  tiplabels(pch=22, cex=1.5,tip = which(b$tip.label %in% names(which(comm [i, ] > 0))))}

#------------------------------------------------------------------------
##################Skyline plots with phylodyn##################
#------------------------------------------------------------------------
}
#------------------------------------------------------------------------
##Extraction of the spatio-temporal information contained in 1,000 trees
#------------------------------------------------------------------------

library(diagram)
library(vioplot)
library(lubridate)
library(OutbreakTools)

#----------------------------
#subsample 1000 trees from the posterior
#----------------------------

allTrees = scan(file="FFVenvPol_relaxedC _GMRF_3.trees.txt", what="", sep="\n", quiet=TRUE)
burnIn = 10001; index1 = which(grepl("tree STATE_0 ",allTrees)) + burnIn
samplingFrequency = (length(allTrees)-index1)/1000
allTrees1 = allTrees[1:(which(grepl("tree STATE_0 ",allTrees))-1)]
allTrees2 = c(allTrees[seq(((index1-1)+samplingFrequency),length(allTrees),samplingFrequency)],"End;")
#can write this as an object which could be useful
write(c(allTrees1,allTrees2), "FFVenvPol_relaxedC _GMRF_3_1000.trees")


#----------------------------
#read Multiphy object
#----------------------------

beast1000FFV <- read.nexus('FFVenvPol_relaxedC _GMRF_3_1000.trees') #this file can be provided at request

beast1000 <- read.nexus('FIV_DG_cauchy_1000.trees')
get_taxa_name(tree_view = NULL, node = NULL)

#####change names on a multiphylo object
nn <- read.csv('newNames.csv', head=F)

beast1000renamed<-lapply(beast1000,taxa_rename, name=nn)
plot(beast1000renamed[[1]])


#------------------------------------------------------------------------
##################MultiPhyloPD##################
#------------------------------------------------------------------------
library(PhyloMeasures)

MultiPhyloPD <- function(ph, comm,...){
  out <- list(NA, ncol = nrow(comm), nrow = length(ph)) 
  
  for(i in 1:length(ph)){
    prunedphy <- prune.sample(comm,ph[[i]])
    pd <- pd.query(prunedphy, comm, standardize = TRUE, 
                   null.model="uniform", reps=1000)
    
    out[[i]] <- pd
    
  }
  out <- t(as.data.frame(out))
  colnames(out) <- row.names(comm)
  row.names(out) <- c(1:length(ph))
  as.data.frame(out)
  
  
}

pFFV <- MultiPhyloPD(beast1000FFV, comm) #has to be the transpose (done above)


#make a tidy object
df <- gather(pFFV, key, value) 
df$key <-as.factor(df$key)

#plot
plot <- ggplot(df, aes(y=value,x=key))
plot + geom_boxplot()+labs(title = 'FFV WS', y="SES.PD", x = 'Contast')+ theme_bw()
