#-----------------------------------------------#
#####Part B: Code to perform simulation on transmission trees#####
#-----------------------------------------------#

############################################################################################################ 

#Treatment population

############################################################################################################ 

# Load libraries
library(fitdistrplus)
library(ergm)
library(intergraph)
library(stringr) # for "str_split"
library(ggplot2)
library(reshape2)
library(tidyverse)
library(igraph)

#-----------------------------------------------------
#calculate target statistics
#-----------------------------------------------------

load( 'WStranTree')
tre <-computeMatWIW(recordWSA3, burnin = 0.1)
G_treatment <- as.directed(graph.adjacency(tre, weighted = T))

#read in transphylo graph
el <-  as.data.frame(cbind( get.edgelist(G_treatment) , prob= round( E(G)$weight, 3 )))

el$prob <- as.numeric( el$prob)

deg <- degree(G_treatment, mode="all")


### subset edgelist by probability cutoff
el2 <- subset(el, el$prob>0.01) #remove all 0 prob edges

el2$prob <- as.numeric( el2$prob)


### fit distribution to UDOI
z1 <- fitdist(el2$prob, "gamma") # gamma distribution consistently fits the best

### code used to determine which distribution consistently fits edge weight data the best
z2 <- fitdist(el2$prob, "exp")
z3 <- fitdist(el2$prob, "lnorm")

deg <- as.vector(deg)

d1 <- fitdist(deg, "pois")
d2 <- fitdist(deg, "nbinom")

hist(el2$prob, freq = F, breaks = 10)
points(x= seq(0,max(el2$prob), by = 0.001), y = dgamma(seq(0,max(el2$prob), by = 0.001), shape = z1$estimate[1], rate = z1$estimate[2]),
       col="darkgreen", lwd=2, yaxt="n")
points(x= seq(0,max(el2$prob), by = 0.001), y = dexp(seq(0,max(el2$prob), by = 0.001), rate = z2$estimate),
       col="darkblue", lwd=2, yaxt="n")
points(x= seq(0,max(el2$prob), by = 0.001), y = dlnorm(seq(0,max(el2$prob), by = 0.001), meanlog = z3$estimate[1], sdlog = z3$estimate[2]),
       col="purple", lwd=2, yaxt="n")

hist(deg, freq = F, breaks = 10)
points(x= seq(0,max(deg), by = 1), y = dpois(seq(0,max(deg), by = 1), lambda = d1$estimate[1]),
       col="darkgreen", lwd=2, yaxt="n")

points(x= seq(0,max(deg), by = 1), y = dpois(seq(0,max(deg), by = 1), lambda = d2$estimate[1]),
       col="purple", lwd=2, yaxt="n")

z1$aic #gamma fits the best for both WS and FR data
z2$aic
z3$aic

d1$aic #this is a better fit than neg binomial for FR
d2$aic#this is a better fit than npoisson for WS

### store parameters

shape <- z3$estimate[1] #for simplicity
rate <- z3$estimate[1] #best for the WS

deg <- d2$estimate[1]

temp.dens <- edge_density(G_treatment, loops = FALSE)
isos <- 4 #for larger population 25
net.size <- 17 # 17 for WS  15 for FR 


#---------------------------------------------------------------------------------------------------
# Function for simulation networks based on observed  transmission network in the treatment region


  # must unload igraph to assign edge attributes
  if( "NetPathMiner" %in% (.packages())) {
    detach("package:NetPathMiner", unload=TRUE) 
  }
  
  if( "netseg" %in% (.packages())) {
    detach("package:netseg", unload=TRUE) 
  }
    if( "igraph" %in% (.packages())) {
      detach("package:igraph", unload=TRUE) 
  }
  
  ity = 1000
  
  net_list <- list()
  
 for (i in 1:ity) {

  

  sim_norm.deg.dist <- rpois(n = net.size-isos, lambda=deg )
  # add isolates and "un-normalize" degree by multiplying by (n-1) where n is the number of vertices (nodes) in the graph
  sim_deg.dist <- c(rep(0, isos), round(sim_norm.deg.dist))
  

  # simulate network based on target degree distribution
  targets <- as.vector(table(sim_deg.dist))
  rand.degs <- sort(unique(sim_deg.dist))
  
  # for small networks, make sure simulated network is of reasonable density
# tries to maintain degree distribution AND density for large networks

      x <- network(net.size, directed = FALSE, density = temp.dens)
      g <- san(x~degree(rand.degs), target.stats = targets, constraints = ~edges) # tries to maintain mean degree

      
      # determine if network density is constrained by seasonal values or overall mean
      # if(network.density(g)>(0.75*temp.dists$net.den) & network.density(g)<(1.25*temp.dists$net.den)) break
      #if(network.density(g)>(0.75*temp.dens) & network.density(g)<(1.25*temp.dens)) break
    #}
  
  #}

  
  #sim_prob <- rgamma(n = network.edgecount(g), shape = shape, rate = rate)
  sim_prob <- rexp(n = network.edgecount(g), rate = rate)
  set.edge.attribute(g, "weight", sim_prob, e=seq_along(g$mel))
  

  #bernouli trials
  sim_sex <- rbinom(net.size, 1,.5) 
   sim_sex[sim_sex==0] <- 'Female'
   sim_sex[sim_sex==1] <- 'Male'
   
  
 g%v%"sex" <- sim_sex #add sex as an attribute
 
 name <- paste('graph:',i,sep='')
 net_list[[name]] <- g
 
 }

##############################################
  #Extract data
##############################################
  
  #double check added correctly
summary.network(net_list$`graph:1`, # the network we want to look at
  print.adj = FALSE) # if TRUE then this will print out the whole adjacency matrix.
  
  node_colors <- rep("",net.size)
   for(i in 1:net.size){
    if(get.node.attr(g,"sex")[i] == "Male"){
     node_colors[i] <- "lightblue"
  }else{
   node_colors[i] <- "maroon"
  }
  }
  plot.network(net_list$`graph:3`, # our network object
  vertex.col = node_colors)

  library(igraph)
  
  #need to convert back to igraph objects
  
  net_listIG <- lapply(net_list, asIgraph)
  
  ############################################## 
  #Function to look at subgraph homophily
  ##############################################
  
  subgraph_edges_homophily <- function(graph, vattr_name, heterophily = FALSE,
                                       drop_isolates = FALSE) {
    stopifnot( # arg checks
      igraph::is.igraph(graph) || is.character(vattr_name) || 
        length(vattr_name) == 1L || !is.na(vattr_name) || 
        vattr %in% igraph::vertex_attr_names(vattr_name)
    )
    
    vattrs <- igraph::vertex_attr(graph, name = vattr_name)
    total_el <- igraph::as_edgelist(graph, names = FALSE)
    
    # rows from total_el where the attribute of the edge source == attribute of edge target
    edges_to_keep <- vattrs[total_el[, 1L]] == vattrs[total_el[, 2L]]
    
    # for heterophilous ties, just negate the "in_group" version
    if (heterophily) edges_to_keep <- !edges_to_keep
    
    igraph::subgraph.edges(graph, 
                           eids = which(edges_to_keep), 
                           delete.vertices = drop_isolates)
  }
  
  ##############################################
  
  #needs an adjacency matrix
  
  rdata <- c()
  
  for (i in 1: length(net_listIG)){
  
  hphly <- net_listIG[[i]]  %>% 
    subgraph_edges_homophily(vattr_name = "sex", heterophily = FALSE) %>% 
    strength(mode= 'out') %>% 
    as.data.frame()

  malesHphly  <-cbind(strength=hphly, sex=V(net_listIG[[i]])$sex ) %>% 
    group_by(sex) %>% 
    summarise(average = mean(.))
  
  rdata[i] <- c(malesHphly[2,2] , rdata) 
  }
  
  rdataU <- unlist(rdata)
  
  hist(rdataU, breaks = 8)

  str(rdata)
  
  
  #calculate probabilities (bootstrap p) not quite right
  ity= 1000
  (sum(rdataU >=   0.0, na.rm=T))/ity
  #################################################################################################
 
  
  ############################################################################################################ 
  
  #Stable population
  
  ############################################################################################################ 
  
  #double check added correctly
  summary.network(net_list$`graph:1`, # the network we want to look at
                  print.adj = FALSE) # if TRUE then this will print out the whole adjacency matrix.
  
  node_colors <- rep("",net.size)
  for(i in 1:net.size){
    if(get.node.attr(g,"sex")[i] == "Male"){
      node_colors[i] <- "lightblue"
    }else{
      node_colors[i] <- "maroon"
    }
  }
  plot.network(net_list$`graph:3`, # our network object
               vertex.col = node_colors)
  
  library(igraph)
  
  #need to convert back to igraph objects
  
  net_listIG <- lapply(net_list, asIgraph)
  #needs an adjacency matrix
  
  rdata <- c()
  
  for (i in 1: length(net_listIG)){
    
    hphly <- net_listIG[[i]]  %>% 
      subgraph_edges_homophily(vattr_name = "sex", heterophily = FALSE) %>% 
      strength(mode= 'out') %>% 
      as.data.frame()
    
    malesHphly  <-cbind(strength=hphly, sex=V(net_listIG[[i]])$sex ) %>% 
      group_by(sex) %>% 
      summarise(average = mean(.))
    
    rdata[i] <- c(malesHphly[2,2] , rdata) 
  }
  
  rdataU <- unlist(rdata)
  
  hist(rdataU, breaks = 8)
  
  str(rdata)
  
  
  #calculate probabilities (bootstrap p) not quite right
  ity= 1000
  (sum(rdataU >=   0.0, na.rm=T))/ity
  
  