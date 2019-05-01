#############################################################################
######SPECIES LOSS SIMULATION#########################################
#######################################################################

##code for Frontiers of Biogeography ms

##aim: to populate islands with realistic species set and distribution of trait values, and
##to simulate the effects of extinctions.


####PACKAGES
require(FD)
require(picante) # also loads vegan and ape
require(BAT)
require(EcoSimR)
require(sars)
require(CommEcol)
require(phytools)
require(plotrix)
require(vegan)
require(ggplot2)

############################################
###INTERNAL FUNCTIONS##########################
###################################################

##work out if a species is a single-island endemic
SIE <- function(s, isl){
    ui <- unlist(isl)
    uis <- sum(ui == s)
    res <- uis == 1
    return(res)
}


##make dendogram for each island as tree object
dendo <- function(x){
    if (nrow(x) > 1) {
        r4 <- apply(x, 2, scale)
        rownames(r4) <- rownames(x)
        a4 <- vegan::vegdist(r4,method = "euclidean")
        a5 <- hclust(a4, method = "average")
        a6 <- ape::as.phylo(a5)
        a6$tip.label <- rownames(r4) 
        a6
    }
}

getk = function(areas, method = "SAR") {
    if (method == "SAR") {
        k = 20 * areas ^ 0.25  #z = 0.25, c = 20
        k = ceiling(k) #carrying capacity
    } else if (method == "habitats") {
        k = areas^(3/2) # k ~ volume of islands
    } else if (method == "mass") {
        k = areas*2 # how much mass an island can carry
    }
    k
}

getpopmass = function(bodymass) {
    density = bodymass ^ -0.22 # from Russo et al. 2003 AmNat
    density * bodymass
}

getdistances = function(species) {
    distmat = as.matrix(dist(species))
    distmat[upper.tri(distmat)] = 0
    id = expand.grid(as.numeric(rownames(distmat)), as.numeric(colnames(distmat)))
    dists = cbind(id , distance=c(distmat))
    dists = dists[dists$distance != 0,]
    dists = dists[order(dists[, 'distance']),]
    dists
}

compete = function(patches, species, areas) {
    k = getk(areas, "SAR")
    popmasses = getpopmass(species[, 'BS'])
    for (p in 1:length(patches)) {
        if (length(patches[[p]]) > 1) {
            dists = getdistances(species[patches[[p]],])
        }
        count = 1
        while (sum(popmasses[patches[[p]]]) > k[p]) {
             if (length(patches[[p]]) > 1) {
                if(count > nrow(dists)) {
                    dists = getdistances(species[patches[[p]],])
                    count = 1
                }
                while(!all(dists[count, 1:2] %in% patches[[p]]) & count < nrow(dists)) {
                    count = count + 1
                }
                if(count > nrow(dists)) {
                    dists = getdistances(species[patches[[p]],])
                    count = 1
                }
                victim = as.numeric(sample(dists[count, 1:2], 1))
                patches[[p]] = patches[[p]][patches[[p]] != victim]
                count = count + 1
            } else {
                patches[[p]] = patches[[p]][NULL]
            }
        }
    }
    patches
}

colonise = function(patches, species, areas, method = "SAR") {
    for(a in 1:length(areas)) {
        remaining = getk(areas[a], method)
        while(remaining > 0) {
            if(length(patches[[a]]) == 0) {
                s = sample(1:nrow(species), 1)
            } else {
                s = sample((1:nrow(species))[-patches[[a]]], 1)
            }                
            patches[[a]] = c(patches[[a]], s)
            if (method == "mass") {
                popsize = getpopmass(species[s, 'BS'])
            } else {
                popsize = 1
            }
            remaining = remaining - popsize
        }
    }
    lapply(patches, unique)
}

#internal function for use in speciate: adds random noise to the trait values
spec_internal <- function(traits){
  newvalues <- rnorm(2, 0, 0.8)
  traits[c(1, 3)] <- traits[c(1, 3)] + newvalues #for body size and beak add random noise
  #TM: the value inside the parentheses can be + or - so no longer always reduces dispersal
  traits[2] <- traits[2] - (traits[2] * rnorm(1, 0, 0.1))#for dispersal - reduce ##LL: too much assumption?
  return(traits)
}

speciate = function(patches, species, rate) {
    for (p in 1:length(patches)) {
        for (s in 1:length(patches[[p]])) {
            if (runif(1) <= rate) {
                traits <- species[s, ]
                dum <- spec_internal(traits)
                #check for negative trait values; if present re-run until not present
                if (any(dum <= 0)){
                  while (any(dum <= 0)) dum <- spec_internal(traits)
                  traits <- dum
                } else {
                  traits <- dum
                }
         #       newvalues <- rnorm(2, 0, 0.8)
           #     traits[c(1, 3)] <- traits[c(1, 3)] + newvalues #for body size and beak add random noise
            #    traits[2] <- traits[2] - (traits[2] * rnorm(1, 0, 0.1))#for dispersal - reduce ##LL: too much assumption?
                ##ensure no traits are negative or 0 values ##TM: occasionally runs indefinitely
         #       while (any(traits <= 0)) {
          #          newvalues <- rnorm(2, 0, 0.8)
           #         traits[c(1, 3)] <- traits[c(1, 3)] + newvalues #for body size and beak add random noise
             #       traits[2] <- traits[2] - (traits[2] * rnorm(1, 0, 0.1))#for dispersal - reduce ##LL: too much assumption?
              #  }
                species <- rbind(species, traits)
                patches[[p]] <- c(patches[[p]], nrow(species))
            }
        }
    }
    rownames(species) = 1:nrow(species)
    list(patches, species)
}
        
            
disperse = function(patches, species = NULL) {
    if (length(patches) > 1) {
        for (p in 1:length(patches)) {
            s = 1
            while (s <= length(patches[[p]])) {
                if (runif(1) <= species[patches[[p]][s], 'D']) {
                    target <- sample(1:length(patches), 1) # picking origin == failed dispersal
                    patches[[target]] = c(patches[[target]], s)
                }
                s = s + 1
            }
        }
    }
    lapply(patches, unique)
}

##############################################################
###MAIN FUNCTION###########################################
################################################################

##creates a mainland pool of 100 species with three traits
##populates the five islands
##calculates the different BG patterns using the pre-human colonisation dataset
##simulates the extinction of th % of species
##re-calculates the different BG patterns
##see Frontiers ms for more details


####Arguments
##plot_T = whether to plot the functional dendogram with extinct species highlighted
##th = the proportion of archipelago species to go extinct (e.g. th = 0.5 = 50% extinction)
##nam = name of the plotted file (for saving)
##verb = print informatino from the various functions run inside Leo

Leo <- function(plot_T = FALSE, plot_F = FALSE, th = 0.5, nam = "Fig_1.jpeg", verb = FALSE){
    
    species <- matrix(nrow = 300, ncol = 3)
    colnames(species) <- c("BS", "D", "Beak")
    rownames(species) = 1:nrow(species)
    
    species[, 1] <- rgamma(nrow(species), 1)#body size
    species[, 2] <- rbeta(nrow(species), 0.9, 1.4)#dispersal
    species[, 3] <- runif(nrow(species), 1, 8)#beak shape

    ##create islands
    isl <- vector("list", length = 5)#list to put island species in
    ar <- c(10, 20, 40, 80, 100)#island areas

    ## metacommunity dynamics:
    isl = colonise(isl, species, ar)
    radiation = speciate(isl, species, 0.1)
    isl = radiation[[1]]
    species = radiation[[2]]
    isl = disperse(isl, species)
    isl = tryCatch(compete(isl, species, ar), error = function(e) NA)
    if (length(isl) == 1) return(NA)

    ##create list with the full trait matrix for each island
    islFull <- lapply(isl, function(x) {
        if (length(species[x]) == ncol(species)) {
            species2 = t(as.matrix(species[x, ]))
        } else {
            species2 = species[x, ]
        }
        if (length(x) > 0) rownames(species2) = x
        species2})
   # print(islFull)
    dendz <- lapply(islFull, function(x) {dendo(x)})

    ##make archipelago dataset and dendogram
    allsp <- unlist(lapply(islFull, rownames))
    allsp2 <- unique(as.numeric(allsp))

    species3 <- species[allsp2,]
    rownames(species3) <- allsp2
    arcDen <- dendo(species3)

    ##calculate PA matrix for archi and each island; calculate PD
    CM <- matrix(0, nrow = nrow(species3), ncol = 1 + length(isl))
    rownames(CM) <- rownames(species3)
    colnames(CM) <- c("A", 1:length(isl))
    CM[ ,1] <- 1
    for (i in 1:length(isl)){
        dd <- isl[[i]]
        mat <- which(rownames(CM) %in% dd)
        CM[mat, (i + 1)] <- 1
    }

    if(plot_F){
      resL <- vector("list", length = 13)
    } else {
      resL <- vector("list", length = 12)
    }

    CM <- t(CM)#picante has species as columns
    resL[[1]] <- picante::pd(CM, arcDen)#get FD of each island
    ##Functional dispersion
    resL[[2]] <- FD::dbFD(x = species3, a = CM, w.abun = FALSE, stand.x = TRUE, messages = verb)

    ##c-score
    ##EcoSimR::c_score(t(CM))#EcoSimR has sites as columns
    ES <- EcoSimR::cooc_null_model(t(CM), algo = "sim9", metric = "c_score", nReps = 1000, suppressProg = !verb)#sim9 is curveball of Strona
    ##get P-values: from ecosim code
    nullmodObj <- ES
    if (nullmodObj$Obs > max(nullmodObj$Sim)){
       LP <- (length(nullmodObj$Sim) - 1)/length(nullmodObj$Sim) 
        UP <- 1/length(nullmodObj$Sim)
    } else if (nullmodObj$Obs < min(nullmodObj$Sim)){
        LP <-  1/length(nullmodObj$Sim)
       UP <- (length(nullmodObj$Sim) - 1)/length(nullmodObj$Sim) 
    } else{
        LP <- format(sum(nullmodObj$Obs >= nullmodObj$Sim)/length(nullmodObj$Sim))
        UP <- format(sum(nullmodObj$Obs <= nullmodObj$Sim)/length(nullmodObj$Sim))
    }

    ses <- format((nullmodObj$Obs - mean(nullmodObj$Sim))/sd(nullmodObj$Sim))

    resL[[3]] <- as.numeric(c(ES$Obs, mean(ES$Sim), LP, UP, ses))
    ##summary(NM)
    
    #taxonomic beta - BAT has species as columns
    resL[[4]] <- BAT::beta.multi(CM)#look at replacement?? ##LL: vegan::betadiver or even zeta diversity?

    ##beta-diversity (functional) - BAT has species as columns
    resL[[5]] <- BAT::beta.multi(CM, arcDen)#look at replacement?? ##LL: vegan::betadiver or even zeta diversity?

    ##sar z value
    sarDf <- data.frame("A" = ar, "S" = rowSums(CM)[2:6])
    resL[[6]] <- sars::sar_power(sarDf)

    
#######HUMAN IMPACT##########################
    ##EXTINCTION PROCESS

    extinct <- c()

    islFull_Ex <- islFull

    ##this iterates in turn through each island which means smaller islands
    ##lose a larger PROPORTION of species, as basically all island lose
    ##roughly same number of species, but smaller have lower starting point
    ##This is as to be expected as smaller island have > extinction rates

    j2 <- 0 
    repeat{
        for (i in 1:5){
            dum <- islFull_Ex[[i]]
            if (nrow(dum) == 1) next #always keep 1 sp on an island ##LL: why?
            ##sample random sp from ith island weighted by body size
            sdum <- sample(1:nrow(dum), 1, prob = dum[, 1])
            samSp <- dum[sdum, , drop = FALSE]
            uid <- as.numeric(rownames(samSp))
            ##check remaining sp to see if is SIE or not
            z1 <- 1:5
            z2 <- z1[-i]
            remSp <- unlist(lapply(islFull_Ex[z2], rownames))
            remSp2 <- unique(as.numeric(remSp))
            if (!uid %in% remSp2){
                extinct <- c(extinct, uid)#if species in no other island add to list (as is SIE)
                dum <- dum[-sdum, , drop = FALSE]#and remove it from islFull for ith island
            } else {
                dum <- dum[-sdum, , drop = FALSE]#just remove it from island (but still exists elsewhere)
            }#eo if
            islFull_Ex[[i]] <- dum
            if (length(extinct) == ceiling(nrow(species3) * th)) break
        }#eo for
        j2 <- j2 + 1
        if (length(extinct) == ceiling(nrow(species3) * th)) break
        if (j2 > 1000) return("NO")
    }#eo repeat

    dead2 <- rep("Present", length = nrow(species3))
    dead2[which(rownames(species3) %in% extinct)] <- "Extinct"

    ##make archipelago dataset and dendogram
    allsp_Ex <- unlist(lapply(islFull_Ex, rownames))
    allsp2_Ex <- unique(as.numeric(allsp_Ex))

    species3_Ex <- species[allsp2_Ex,]
    rownames(species3_Ex) <- allsp2_Ex
    arcDen_Ex <- dendo(species3_Ex)
    ##calculate PA matrix for archi and each island; calculate PD
    CM_Ex <- matrix(0, nrow = nrow(species3_Ex), ncol = 6)
    rownames(CM_Ex) <- rownames(species3_Ex)
    colnames(CM_Ex) <- c("A", 1:5)
    CM_Ex[ ,1] <- 1
    for (i in 1:5){
        dd <- rownames(islFull_Ex[[i]])
        mat <- which(rownames(CM_Ex) %in% dd)
        CM_Ex[mat, (i + 1)] <- 1
    }

    CM_Ex <- t(CM_Ex)#picante has species as columns
    resL[[7]] <- picante::pd(CM_Ex, arcDen_Ex)#get FD of each island
    ##Functional richness
    resL[[8]] <- FD::dbFD(x = species3_Ex , a = CM_Ex , w.abun = FALSE, stand.x = TRUE, messages = verb)#FD has species as columns

    ##c-score
    ##EcoSimR::c_score(t(CM_Ex))#EcoSimR has sites as columns
    ES_Ex <- EcoSimR::cooc_null_model(t(CM_Ex), algo = "sim9", metric = "c_score", nReps = 1000, suppressProg = !verb)#sim9 is curveball of Strona
    ##get P-values: from ecosim code
    nullmodObj <- ES_Ex
    if (nullmodObj$Obs > max(nullmodObj$Sim)){
        LP <- (length(nullmodObj$Sim) - 1)/length(nullmodObj$Sim) 
        UP <- 1/length(nullmodObj$Sim)
    } else if (nullmodObj$Obs < min(nullmodObj$Sim)){
        LP <-  1/length(nullmodObj$Sim)
        UP <- (length(nullmodObj$Sim) - 1)/length(nullmodObj$Sim) 
    } else{
        LP <- format(sum(nullmodObj$Obs >= nullmodObj$Sim)/length(nullmodObj$Sim))
        UP <- format(sum(nullmodObj$Obs <= nullmodObj$Sim)/length(nullmodObj$Sim))
    }

    ses <- format((nullmodObj$Obs - mean(nullmodObj$Sim))/sd(nullmodObj$Sim))

    resL[[9]] <- as.numeric(c(ES_Ex$Obs, mean(ES_Ex$Sim), LP, UP, ses))
    
    ##taxonomic beta - BAT has species as columns
    resL[[10]] <- BAT::beta.multi(CM_Ex)

    ##beta-diversity (functional) - BAT has species as columns
    resL[[11]] <- BAT::beta.multi(CM_Ex, arcDen_Ex)

    ##sar z value
    sarDf_Ex <- data.frame("A" = ar, "S" = rowSums(CM_Ex)[2:6])
    resL[[12]] <- sars::sar_power(sarDf_Ex)

    ##plot the dendogram with extinct species and extant species
    if (plot_T){
        #get dendogram and species list for biggest island
        bigI <- islFull[[5]]
        bigDen <- dendo(bigI)
        bigSp <- rownames(bigI)
        
        #get species names that went extinct on largest island
        bigE <- islFull_Ex[[5]]#sp remaining AFTER extinction
        ESp <- rownames(bigE)
        we <- as.numeric(which(!bigSp %in% ESp))
        dead3 <- rep("Present", length(bigSp))
        dead3[we] <- "Extinct"
        
        fmode <- as.factor(setNames(dead3, bigSp))

        jpeg(paste(nam), width = 20, height = 20, units = "cm", res = 300)
        phytools::dotTree(bigDen, fmode, colors=setNames(c("blue","red"),
                                                         c("Present", "Extinct")))
        dev.off()
    }
    
    #do PCA on trait matrix of all species, then plot functional space
    #of all species and then just extinct species
    if (plot_F){
      pca <- vegan::rda(species3)
      pcaS <- summary(pca)
      tra <- as.data.frame(pcaS$sites[ ,1:2])
      colnames(tra) <- c("A1", "A2")
      
      ct <- apply(tra, 2, function(x) cor(x, species3))
      rownames(ct) <- colnames(species3)
      
      traE <- tra[which(rownames(tra) %in% rownames(species3_Ex)),]
      
      g1 <- ggplot() + geom_point(data = tra, aes(x = A1, y = A2), colour = "blue") + xlab("PCA 1") + ylab("PCA 2") + theme_bw() +
        geom_point(data = traE, aes(x = A1, y = A2), colour = "red", shape  = 17, size = 3, alpha = 0.5) + 
        geom_polygon(data = as.data.frame(
          tra[chull(tra$A1, tra$A2),]), aes(x = A1, y = A2), size = 1,alpha = 0, colour = "blue") + 
        geom_polygon(data = as.data.frame(traE[chull(traE$A1, traE$A2),]), aes(x = A1, y = A2), size = 1.1, alpha = 0, colour = "red") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
      funcSpa <- list(g1, ct)
      resL[[13]] <- funcSpa
    }
    
    return(resL)
}


###################################################################################
#######RUN FUNCTION##############################################
#################################################################

tes <- Leo()

tes <- Leo(plot_F = TRUE)
tes[[13]][[1]]
tes[[13]][[2]]

jpeg("Figure_2.jpeg", width = 15, height = 15, res = 600, units = "cm")
tes[[13]][[1]]
dev.off()


##Run Leo N times, create a list of lists and then format it to provide average results with standard error


library(foreach)
library(doParallel)
library(cluster)

cores = 10
cl = makeCluster(cores); on.exit(stopCluster(cl))
registerDoParallel(cl)
i = 1 #Dummy line for RStudio

Leo2 = foreach(i=seq(from=1, to=100, by=1))  %dopar% { 
  require(FD)
  require(picante) # also loads vegan and ape
  require(BAT)
  require(EcoSimR)
  require(sars)
  require(CommEcol)
  require(phytools)
  require(plotrix)
  
  Leo()
}

#Leo2 <- replicate(100, Leo())

##save(Leo2, file = "Leo2.R")


anyNA(Leo2)


form_leo <- function(x = Leo2){
    
    ##all metrics for pre-colonisation data
    
    F1 <- apply(x, 2, function(y) y[[1]][ ,1])
    F1ue <- rbind(apply(F1, 1, mean), apply(F1, 1, plotrix::std.error))
    colnames(F1ue) <- c("A", 1:5)
    
    F2 <- apply(x, 2, function(y) unlist(y[[2]][c(3, 7)]))
    F2ue <- rbind(apply(F2, 1, mean, na.rm = TRUE), apply(F2, 1, plotrix::std.error))
    
    F3 <- apply(x, 2, function(y) y[[3]])
    F3ue <- rbind(apply(F3, 1, mean), apply(F3, 1, plotrix::std.error))
    colnames(F3ue) <- c("Obs", "Mean", "LP", "UP", "SES")
    
    F4 <- apply(x, 2, function(y) y[[4]][ ,1])
    F4ue <- rbind(apply(F4, 1, mean), apply(F4, 1, plotrix::std.error))
    
    F5 <- apply(x, 2, function(y) y[[5]][ ,1])
    F5ue <- rbind(apply(F5, 1, mean), apply(F5, 1, plotrix::std.error))
    
    F6 <- apply(x, 2, function(y) y[[6]]$par)
    F6ue <- rbind(apply(F6, 1, mean), apply(F6, 1, plotrix::std.error))
    
    
    l1 <- list(F1ue, F2ue, F3ue, F4ue, F5ue, F6ue)
    names(l1) <- c("FD", "FRIC_FDIS", "Chequer", "TaxoBeta", "FuncBeta", "SAR")
    
    ##all metrics for post-colonisation data
    F7 <- apply(x, 2, function(y) y[[7]][ ,1])
    F7ue <- rbind(apply(F7, 1, mean), apply(F7, 1, plotrix::std.error))
    colnames(F7ue) <- c("A", 1:5)
    
    F8 <- apply(x, 2, function(y) unlist(y[[8]][c(3, 7)]))
    F8ue <- rbind(apply(F8, 1, mean, na.rm = TRUE), apply(F8, 1, plotrix::std.error))
    
    F9 <- apply(x, 2, function(y) y[[9]])
    F9ue <- rbind(apply(F9, 1, mean), apply(F9, 1, plotrix::std.error))
    colnames(F9ue) <- c("Obs", "Mean", "LP", "UP", "SES")
    
    F10 <- apply(x, 2, function(y) y[[10]][ ,1])
    F10ue <- rbind(apply(F10, 1, mean), apply(F10, 1, plotrix::std.error))
    
    F11 <- apply(x, 2, function(y) y[[11]][ ,1])
    F11ue <- rbind(apply(F11, 1, mean), apply(F11, 1, plotrix::std.error))
    
    F12 <- apply(x, 2, function(y) y[[12]]$par)
    F12ue <- rbind(apply(F12, 1, mean), apply(F12, 1, plotrix::std.error))
    
    
    l2 <- list(F7ue, F8ue, F9ue, F10ue, F11ue, F12ue)
    names(l2) <- c("FD", "FRIC_FDIS", "Chequer", "TaxoBeta", "FuncBeta", "SAR")
    
    ##signficant tests
    t1 <- t.test(F1[1,], F7[1,]) 
    t2 <- t.test(F2[1,], F8[1,]) 
    t3 <- t.test(F3[1,], F9[1,]) 
    t4 <- t.test(F4[1,], F10[1,]) 
    t5 <- t.test(F5[1,], F11[1,])
    t6 <- t.test(F6[2,], F12[2,])
    
    ta <- list(t1, t2, t3, t4, t5, t6)
    tt <- vapply(ta, function(y){
      y1 <- y$statistic %>% round(2)
      y2 <- y$p.value %>% round(2)
      c(y1, y2)
    }, FUN.VALUE = numeric(2))
    
    l3 <- list(l1, l2, tt)
    names(l3) <- c("Pre-humans", "Post-humans", "T_tests")
    
    return(l3)  
}

TSL <- form_leo()













