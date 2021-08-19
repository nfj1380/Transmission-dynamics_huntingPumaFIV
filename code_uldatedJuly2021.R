
#-----------------------------------------------#
#####Transmission tree code#####

#by Nick Fountain-Jones, November 2018
#-----------------------------------------------#

#install if necesary

#library(devtools)

#devtools::install_github('xavierdidelot/TransPhylo')

#-----------------------------------------------#
#####Prelims#####
#-----------------------------------------------#

library(TransPhylo)
library(ape)
library(coda)
library(intergraph)
library(tidyverse)

library(treeio)
library(ggtree)

#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("ggtree")

#library(Rcpp)

rm(list=ls())
set.seed(0)


#-----------------------------------------------#
#####Reading in the phylogeny####
#-----------------------------------------------#


#------------------------------------------------------------------------
##################Visualize MCC Tree##################
#------------------------------------------------------------------------
#FIV - try with a .tre file

beastTree<- read.beast('CO_FIV_Gamma_LogC_ExpGrowth_Cauchy3_Pop3.tre')

#plot the tree

phy_simple<- ggtree(beastTree) + geom_tiplab(size=2)+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 

phy_simple

#prune each dominant clade in each region. Based on FOuntain-Jones etal 2020

Gamma_treatment <- tree_subset(beastTree , node=66,levels_back = 0)
Gamma_stable <- tree_subset(beastTree , node=43,levels_back = 0)

plot()
#-----------------------------------------------#
#####Treatment Region#####
#-----------------------------------------------#



#-----------------------------------------------#
#####Load in the phylogenetic tree#####
#-----------------------------------------------#

#load tree
#ptree<-ptreeFromPhylo(read.tree('WSFIVGammaTree.newick'),dateLastSample=2014.115) #read.tree is also an option.

ptree<-ptreeFromPhylo(as.phylo(Gamma_treatment),dateLastSample=2014.115)

plot(ptree) #check it worked
#IMPORTANT callibration step - need a reasonable distribution of generation time. We need to discuss this with team virus.


w.shape=2 # w.shape=1,w.scale=0.5 = 6 months on average between transmission events. Ebroch et al? #scale=1,w.shape=1 we get an Exponential distribution with mean 1 year. #scale=1,w.shape=5 we get an normal distribution with mean 5years. Gamma 2/2 shape/scale come from Courchamp et al FIV paper 1995 mean of around 4 years
w.scale=1.5 #was 1.5

#sample time (from infection to sample) The chances we sampled an individual right after infection are low, but more likely ~4 years post infection
ws.shape=2
ws.scale=1.5 #was 1.5

startPi <- 0.6 # proportion sampled

#date <- lubridate::ymd("2014-02-12")
#lubridate::decimal_date(date)

dateT=2014.117 #end of sampling date +1 day
#dateT=Inf # Epidemic over - can be good to test this assumption out.


#-----------------------------------------------#
#####Calculate transmission tree#####
#-----------------------------------------------#

#update pi.
ttree_treatment<-inferTTree(ptree,mcmcIterations=200000, thin = 2, w.shape=w.shape,ws.scale=ws.scale,ws.shape=ws.shape, w.scale=w.scale,dateT=dateT, startPi = startPi, updatePi = T, updateOff.p=F)

save(ttree_treatment, file='ttree_treatment')


load('ttree_treatment')
#set pi prior- doesn't work yet
#prior_pi_a=2 #this should focus to prior more towards 50% sampling with 0% and 100% rare
#prior_pi_b=2

#infer_multittree_share_param(ptree,mcmcIterations=10000,prior_pi_b=prior_pi_b,prior_pi_a=prior_pi_a, w.shape=w.shape,w.scale=w.scale,dateT=dateT, startPi = startPi, updatePi = T)


#check for convergence

mcmc <- convertToCoda(ttree_treatment)
effectiveSize(mcmc)

#check out mcmc diagnostics. WIll need to add some more iterationss

par(mfrow=c(2,2))
plot(sapply(ttree_treatment,function(x) x$pTTree+x$pPTree),ylab='Posterior probability',
     xlab='MCMC iterations',type='l')
plot(sapply(ttree_treatment,function(x) x$pi),ylab='Sampling proportion pi',
     xlab='MCMC iterations',type='l')
plot(sapply(ttree_treatment,function(x) x$neg),ylab='Within-host coalescent rate Ne*g',
     xlab='MCMC iterations',type='l')
plot(sapply(ttree_treatment,function(x) x$off.r),ylab='Basic reproduction number R',
     xlab='MCMC iterations',type='l')

#construct sumamry tree
cons <- consTTree(ttree_treatment)

par(mfrow=c(1,1))

plotTTree2(cons,w.shape,w.scale)
plotTTree(cons,w.shape,w.scale)


#calculate R2

R2 <- sapply(ttree_treatment,function(x) x$off.r)

summary(R2)

#within host
NeG <- sapply(ttree_treatment,function(x) x$neg)

summary(NeG)

#Other diagnotic plots

getIncidentCases(ttree_treatment,show.plot = T, burnin = 0.1)


getSamplingTimeDist(ttree_treatment,show.plot = T)

getGenerationTimeDist(ttree_treatment,show.plot = T)


cons$nam
d <- getInfectionTimeDist(ttree_treatment,k='x1952_2004-05-28',show.plot = T)
as.data.frame(d)
mean(d)

#---------------------------------------------------------------------------------

#Individual infection times
d <- getInfectionTimeDist(ttree_treatment,k='X139R1_WS_2009-05-01',show.plot = T)
e <- getInfectionTimeDist(ttree_treatment,k='X1131_WS_2010-24-02',show.plot = T)
f <- getInfectionTimeDist(ttree_treatment,k='X285_WS_2008-21-02',show.plot = T)
g <- getInfectionTimeDist(ttree_treatment,k='X1380_WS_2011-09-02',show.plot = T)
h <- getInfectionTimeDist(ttree_treatment,k='X1341_WS_2010-12-11',show.plot = T)
i <- getInfectionTimeDist(ttree_treatment,k='X1370_WS_2011-20-01',show.plot = T)
j <- getInfectionTimeDist(ttree_treatment,k='X1130_WS_2010-27-02',show.plot = T)


#-----------------------------------------------#
#####Code to drop tip simulatee###
#-----------------------------------------------#

#Treatment region

drop_tip_ptree <- function(tree){
  i <- sample(1:length(tree$tip.label), 1)
  ape::drop.tip(tree, tip=i)
}


tree<- as.phylo(Gamma_treatment)

#copy the object across 10 iterations
tree_list <- list(tree, tree, tree, tree, tree, tree, tree, tree, tree, tree) #better way to do this - but it is touchy


treelist_dropTip <- lapply(tree_list, drop_tip_ptree )

plot(treelist_dropTip[[3]]) #sanity check
#

trans_tree_extract <- function(treeList, dateLastSample = 2014.115 ){
  
  ptree<-ptreeFromPhylo(tree,dateLastSample) #read.tree is also an option.
  
  
  record<-inferTTree(ptree,mcmcIterations=100000, thin = 2, w.shape=w.shape,ws.scale=ws.scale,ws.shape=ws.shape, w.scale=w.scale,dateT=dateT, startPi = startPi, updatePi = T, updateOff.p=F)
}


trans_treeList <- lapply(treelist_dropTip, trans_tree_extract )
str(tree )
trans_treeList[[1]]

R2list <- function(trans_treeList){
  R2 <- sapply(trans_treeList,function(x) x$off.r)
  summary(R2)
}

Summary_R2_treatment <- lapply(trans_treeList, R2list)
saveRDS(Summary_R2_treatment, 'Summary_R2_treatment')


#---------------------------------------------------------------------------------
#-----------------------------------------------#
#####Turn into a graph#####
#-----------------------------------------------#

library(igraph)

tre_treatment <-computeMatWIW(ttree_treatment, burnin = 0.1)

#-----------------------------------------------#
#####Turn into a graph#####
#-----------------------------------------------#

library(network)
library(ergm)
source('igraphArrowHack.R') #need this

environment(plot.igraph2) <- asNamespace('igraph')
environment(igraph.Arrows2) <- asNamespace('igraph')

G <- as.directed(graph.adjacency(tre_treatment, weighted = T)) #make a directed graph
deg <- degree(G, mode="all") #degree

G <- igraph::delete.edges(G, which(E(G)$weight <=.05)) #removed some of the less supported links
V(G)$size <- deg*3 #makes node size more visible
plot.igraph2(G, edge.arrow.size =abs(E(G)$weight), vertex.size=(1*deg),  edge.color=c("red","green")[sign(E(G)$weight)],
             edge.width = 10 *abs(E(G)$weight), vertex.size=1+(1*deg), vertex.color='gold', edge.curved=0.6, edge.arrow.size=5 *abs(E(G)$weight), vertex.label.dist = 1, vertex.frame.color = 'gray', vertex.label.color='black')

write_graph(G, file = 'CO_treatment_network', format = c("ncol"))


#add sex as a trait
V(G)$name #need this to make sure sex is in order
obs_sex_treatment <- c('Male','Male', 'Male', 'Female', 'Female', 'Male', 'Male', 'Male', 'Male', 'Male','Male',
                'Male', 'Female', 'Male', 'Female', 'Male', 'Male'  )

Gnet <- asNetwork(G)

Gnet%v%"sex" <- obs_sex_treatment

net.size=17

#PLOT
node_colors <- rep("",net.size)
for(i in 1:net.size){
  if(get.node.attr(Gnet,"sex")[i] == "Male"){
    node_colors[i] <- "lightblue"
  }else{
    node_colors[i] <- "maroon"
  }
}
plot.network(Gnet, # our network object
             vertex.col = node_colors)


#needs an adjacency matrix

GnetIg <- asIgraph(Gnet)

GnetIgUSimple <- simplify( GnetIg, edge.attr.comb = max) #sum all edges linking individuals
plot( GnetIgUSimple )

source('subgraph_edges_homophily.R') #need this

hphly <- GnetIgUSimple  %>%
  subgraph_edges_homophily(vattr_name = "sex", heterophily = FALSE) %>%  #change heterophyly to false
  strength(mode= 'out') %>%
  as.data.frame()

malesHphly  <-cbind(strength=hphly, sex=V(GnetIgUSimple)$sex ) %>%
  group_by(sex) %>%
  summarise(average = mean(.))

malesHphly  
malesHphly[2,2] #average weighted degree for males

#compare just degree

deg <- as.numeric(deg)
Netstrength <- strength( GnetIg )

avgDeg  <-as.data.frame(cbind(degree=deg, sex=V(GnetIgUSimple)$sex ))
avgDeg$degree <- as.numeric(avgDeg$degree)

avgDeg %>%
  group_by(sex) %>%
  summarise(average = mean(degree))

str(avgDeg )

avgStrength <-as.data.frame(cbind(strength=Netstrength , sex=V(GnetIgUSimple)$sex ))
avgStrength$ strength <- as.numeric(avgStrength $ strength)

avgStrength %>%
  group_by(sex) %>%
  summarise(average = mean(strength))


#-----------------------------------------------#

###########################################################
#How does this impact network topology?
###########################################################

#weighted degree homophily

#make a dataframe of vertex names and attributes

V(G)$name #need this to make sure sex is in order
obs_sex_treatment<- c('Male','Male', 'Male', 'Female', 'Female', 'Male', 'Male', 'Male', 'Male', 'Male','Male',
                'Male', 'Female', 'Male', 'Female', 'Male', 'Male'  )

vertexNameTableObs <- data.frame(name=V(G)$name, sex = obs_sex_treatment )

#drop one tip and record weighted homophily

weightedHomophily_dropTip <- function(trans_treeList){

  tre <-computeMatWIW(trans_treeList, burnin = 0.1)
  G <- igraph::as.directed(graph.adjacency(tre, weighted = T))
  G <- igraph::delete.edges(G, which(E(G)$weight <=.05)) 
  plot(G)
  graphNames <- data.frame(names=V(G)$name)

    namesMatch <- vertexNameTableObs %>%

    filter(name %in% graphNames$names)

     sim_sex <- c(namesMatch$sex)
  
  Gnet <- asNetwork(G)
  
  #add sex as an attribute
  Gnet%v%"sex" <-   sim_sex 
  
  #extract adjacency
  GnetIg <- asIgraph(Gnet)
  
  GnetIgUSimple <- simplify( GnetIg, edge.attr.comb = max) #sum all edges linking individuals
  plot( GnetIgUSimple )
  
  hphly <- GnetIgUSimple  %>%
    subgraph_edges_homophily(vattr_name = "sex", heterophily = FALSE) %>%  #change heterophyly to false
    strength(mode= 'out') %>%
    as.data.frame()
  
  malesHphly  <-cbind(strength=hphly, sex=V(GnetIgUSimple)$sex ) %>%
    group_by(sex) %>%
    summarise(average = mean(.))
  
  malesHphly  
  #malesHphly[2,2]
  
  
}

weightedHomophilyList <- lapply(trans_treeList, weightedHomophily_dropTip )



#################################Stable region######################################

#load tree
ptree_stable<-ptreeFromPhylo(as.phylo(Gamma_stable) ,dateLastSample=2013.211) #read.tree is also an option.

plot(ptree_stable)

#IMPORTANT callibration step - need a reasonable distribution of generation time. We need to discuss this with team virus.


w.shape=2 # w.shape=1,w.scale=0.5 = 6 months on average between transmission events. Ebroch et al? #scale=1,w.shape=1 we get an Exponential distribution with mean 1 year. #scale=1,w.shape=5 we get an normal distribution with mean 5years. Gamma 2/2 shape/scale come from Courchamp et al FIV paper 1995 mean of around 4 years
w.scale=1.5

#sample time (from infection to sample) The chances we sampled an individual right after infection are low, but more likely ~4 years post infection
ws.shape=2
ws.scale=1.5

startPi <- 0.6 # proportion sampled

#date <- lubridate::ymd("2013-03-19")
#lubridate::decimal_date(date)

dateT=2014.3 #end of sampling date +1 day
#dateT=Inf # Epidemic over - can be good to test this assumption out.


#-----------------------------------------------#
#####Calculate transmission tree#####
#-----------------------------------------------#

#update pi.
ttree_stable<-inferTTree(ptree_stable,mcmcIterations=200000, thin = 2, w.shape=w.shape,ws.scale=ws.scale,ws.shape=ws.shape, w.scale=w.scale,dateT=dateT, startPi = startPi, updatePi = T, updateOff.p=F)


save(ttree_stable, file='ttree_stable')


load('ttree_stable')
#set pi prior- doesn't work yet
#prior_pi_a=2 #this should focus to prior more towards 50% sampling with 0% and 100% rare
#prior_pi_b=2

#infer_multittree_share_param(ptree,mcmcIterations=10000,prior_pi_b=prior_pi_b,prior_pi_a=prior_pi_a, w.shape=w.shape,w.scale=w.scale,dateT=dateT, startPi = startPi, updatePi = T)


#check for convergence

mcmc <- convertToCoda(ttree_stable)
effectiveSize(mcmc)

#plot tree to see if it worked
#lastIteration<-recordWSA[[length(recordWSA)]]
#plotCTree(lastIteration$ctree)

#check out mcmc diagnostics. WIll need to add some more iterationss

par(mfrow=c(2,2))
plot(sapply(ttree_treatment,function(x) x$pTTree+x$pPTree),ylab='Posterior probability',
     xlab='MCMC iterations',type='l')
plot(sapply(ttree_treatment,function(x) x$pi),ylab='Sampling proportion pi',
     xlab='MCMC iterations',type='l')
plot(sapply(ttree_treatment,function(x) x$neg),ylab='Within-host coalescent rate Ne*g',
     xlab='MCMC iterations',type='l')
plot(sapply(ttree_treatment,function(x) x$off.r),ylab='Basic reproduction number R',
     xlab='MCMC iterations',type='l')

#construct sumamry tree
cons_stable <- consTTree(ttree_stable)

par(mfrow=c(1,1))

plotTTree2(cons_stable,w.shape,w.scale)
plotTTree(cons_stable,w.shape,w.scale)

tre_stable <-computeMatWIW(ttree_stable, burnin = 0.1);str(tre)

#calculate R2

R2_stable <- sapply(ttree_stable,function(x) x$off.r)

summary(R2_stable)

#within host
NeG <- sapply(recordWSA3,function(x) x$neg)

summary(NeG)

#Other diagnotic plots

a <- getIncidentCases(recordWSA2,show.plot = T, burnin = 0.1)


b <- getSamplingTimeDist(recordWSA2,show.plot = T)

c <- getGenerationTimeDist(recordWSA2,show.plot = T)

str(cons)
cons$nam
d <- getInfectionTimeDist(recordWSA2,k='x1952_2004-05-28',show.plot = T)
as.data.frame(d)
mean(d)

a1 <- getInfectionTimeDist(recordFR,k='X1077_FR_2010-28-01',show.plot = T) #F52
b1 <- getInfectionTimeDist(recordFR,k='X1552_FR_2011-24-03',show.plot = T) #AM83
c1 <- getInfectionTimeDist(recordFR,k='X298R1_FR_2009-18-03',show.plot = T) #AF19

d1 <- getInfectionTimeDist(recordFR,k='X1647_WS_2012-23-01',show.plot = T) #M161 born WS22 April 2011
e1 <- getInfectionTimeDist(recordFR,k='X1674_FR_2013-04-01',show.plot = T) #AM98 born speing 2012

f1 <- getInfectionTimeDist(recordFR,k='X1064_FR_2010-14-01',show.plot = T) #AF54
g1 <- getInfectionTimeDist(recordFR,k='X299_FR_2009-06-03',show.plot = T) #AM20


#-----------------------------------------------#
#####Code to drop tip simulate#####
#-----------------------------------------------#
tree<- as.phylo(Gamma_stable)

plot(tree)

#copy the object across 10 iterations
tree_list <- list(tree, tree, tree, tree, tree, tree, tree, tree, tree, tree) #better way to do this - but it is touchy


treelist_dropTip <- lapply(tree_list, drop_tip_ptree )

plot(treelist_dropTip[[3]]) #sanity check
#

trans_tree_extract <- function(treeList, dateLastSample = 2013.211 ){
  
  ptree<-ptreeFromPhylo(tree,dateLastSample) #read.tree is also an option.
  
  record<-inferTTree(ptree,mcmcIterations=10000,w.shape=w.shape,w.scale=w.scale,dateT=dateT, startPi = startPi, updatePi = T)
  
}


trans_treeList_stable <- lapply(treelist_dropTip, trans_tree_extract )
str(tree )

R2list <- function(trans_treeList_stable){
  R2 <- sapply(trans_treeList,function(x) x$off.r)
  summary(R2)
}

Summary_R2_stable <- lapply(trans_treeList, R2list)

saveRDS(Summary_R2_stable, 'Summary_R2_stable')


################################ Stable region networks ######################################

G <- as.directed(graph.adjacency(tre_stable, weighted = T)) #make a directed graph
deg <- degree(G, mode="all") #degree

G <- igraph::delete.edges(G, which(E(G)$weight <=.05)) #removed some of the less supported links
V(G)$size <- deg*3 #makes node size more visible
plot.igraph2(G, edge.arrow.size =abs(E(G)$weight), vertex.size=(1*deg),  edge.color=c("red","green")[sign(E(G)$weight)],
             edge.width = 10 *abs(E(G)$weight), vertex.size=1+(1*deg), vertex.color='gold', edge.curved=0.6, edge.arrow.size=5 *abs(E(G)$weight), vertex.label.dist = 1, vertex.frame.color = 'gray', vertex.label.color='black')


#add sex as a trait
V(G)$name #need this to make sure sex is in order
obs_sex_stable <- c('Female','Male', 'Female', 'Female', 'Female', 'Female', 'Male', 'Female', 'Female', 'Male', 'Female', 'Female',
                'Female', 'Female', 'Male', 'Female', 'Male')
Gnet <- asNetwork(G)

Gnet

Gnet%v%"sex" <-obs_sex_stable 

net.size=16

#CHECK IT WORKED
node_colors <- rep("",net.size)
for(i in 1:net.size){
  if(get.node.attr(Gnet,"sex")[i] == "Male"){
    node_colors[i] <- "lightblue"
  }else{
    node_colors[i] <- "maroon"
  }
}
plot.network(Gnet, # our network object
             vertex.col = node_colors)

######################################
#calculate weighted degree for each sex
library(igraph)

GnetIg <- asIgraph(Gnet)

GnetIgUSimple <- simplify( GnetIg, edge.attr.comb = max) #sum all edges linking individuals
plot( GnetIgUSimple )

source('subgraph_edges_homophily.R')

hphly <- GnetIgUSimple  %>%
  subgraph_edges_homophily(vattr_name = "sex", heterophily = FALSE) %>%  #change heterophyly to false
  strength(mode= 'out') %>%
  as.data.frame()

malesHphly  <-cbind(strength=hphly, sex=V(GnetIgUSimple)$sex ) %>%
  group_by(sex) %>%
  summarise(average = mean(.))

malesHphly  


#compare just degree

deg <- as.numeric(deg)
Netstrength <- strength( GnetIg )

avgDeg  <-as.data.frame(cbind(degree=deg, sex=V(GnetIgUSimple)$sex ))
avgDeg$degree <- as.numeric(avgDeg$degree)

avgDeg %>%
  group_by(sex) %>%
  summarise(average = mean(degree))

str(avgDeg )

avgStrength <-as.data.frame(cbind(strength=Netstrength , sex=V(GnetIgUSimple)$sex ))
avgStrength$ strength <- as.numeric(avgStrength $ strength)

avgStrength %>%
  group_by(sex) %>%
  summarise(average = mean(strength))

###########################################################
#How does this impact network topology?
###########################################################

#weighted degree homophily

#make a dataframe of vertex names and attributes

V(G)$name #need this to make sure sex is in order
obs_sex_stable <- c('Male', 'Female', 'Female', 'Female', 'Male', 'Female',  'Male', 'Male', 'Female', 'Female','Female',
                    'Male', 'Female', 'Female', 'Female', 'Male')


vertexNameTableObs <- data.frame(name=V(G)$name, sex = obs_sex_stable )

#drop one tip and record weighted homophily


weightedHomophilyList_stable <- lapply(trans_treeList_stable , weightedHomophily_dropTip )
