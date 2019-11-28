#############################################################################
######SPECIES LOSS SIMULATION#########################################
#######################################################################

##code for Ecological Research ms

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
library(dplyr)

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
  
   #remove tolerance trait
   species <- species[,-4]
  
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
  species <- as.data.frame(species)
    for(a in 1:length(areas)) {
        remaining = getk(areas[a], method)
        while(remaining > 0) {
            if(length(patches[[a]]) == 0) {
                s = sample(1:nrow(species), 1, prob = species$D)
            } else {
                w = unique(patches[[a]])
                s2 = species$D[-w]
                s = sample((1:nrow(species))[-patches[[a]]], 1, prob = s2)
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
            #bug correction
            s2 <- patches[[p]][s]
            traits <- species[s2, ]
            t4 <- traits[4]
            rate2 <- rate * t4 #caluclate speciatino rate based on physiological tolerance tait
            if (runif(1) <= rate2) {
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

##creates a mainland pool of 300 species with three traits
##populates the five islands
##calculates the different BG patterns using the pre-human colonisation dataset
##simulates the extinction of th % of species
##re-calculates the different BG patterns
##see ms for more details


####Arguments
##plot_T = whether to plot the functional dendogram with extinct species highlighted
##plot_F = plot functional space of pre- and post- in PCA space
##th = the proportion of archipelago species to go extinct (e.g. th = 0.5 = 50% extinction)
##bs_I = logical - include body size in functional diversity and functional beta calculations
##Ext_method = method to use for extinctions, can be one of: 1) "stan" = standard way of just iterating
##across islands until threshold reached (means all islands lose similar numbers of sp, but small islands
##lose greater proportion); 2) "prob" = iterates across islands but this time a binomial distribution used
##to detemrine whether an extinction event occurs, and the probability an extinction events occurs is related 
##to island size; 3) "Lud1" = species randomly selected from pool, and then an island where this sp occurs
##randomly selected and this population goes extinct; 4) "Lud2" = same as 3 (Lud1) but the island is chosen based
##on island area.
##nam = name of the plotted file (for saving)
##verb = print information from the various functions run inside Leo


#plot_T = FALSE
#plot_F = FALSE
#th = 0.5
#bs_I = FALSE
#Ext_method = "lud2"
#nam = "Fig_1.jpeg"
#verb = FALSE


Leo <- function(plot_T = FALSE, plot_F = FALSE, th = 0.5, bs_I = FALSE,
                Ext_method = "stan", nam = "Fig_1.jpeg", verb = FALSE){
  
  if (!Ext_method %in% c("stan", "prob", "lud1", "lud2")) stop("Ext_method needs to be one of: stan, prob, lud1 or lud2")
    
    species <- matrix(nrow = 300, ncol = 4)
    colnames(species) <- c("BS", "D", "Beak", "tolerance")
    rownames(species) = 1:nrow(species)
    
    species[, 1] <- rgamma(nrow(species), 1)#body size
    species[, 2] <- rbeta(nrow(species), 0.9, 1.4)#dispersal
    species[, 3] <- runif(nrow(species), 1, 8)#beak shape
    species[, 4] <- rbeta(nrow(species), 0.8, 1.8)#tolerance
    
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
    
    # REMOVING BODY SIZE
    if (!bs_I) species3 <- species3[,-1]
    
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
    #scaling traints to mean 0, sd = 1 makes no difference to FD calculation (checked)
    resL[[1]] <- picante::pd(CM, arcDen)#get FD of each island
    ##Functional richness
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
    
    
    if (Ext_method == "lud1" || Ext_method == "lud2"){
    
    #########################################################################################
    ##ludwig method of extinction: pick a species from the pool (weighted by bs),
    #and then pick a random island popn. (Lud1; i.e. a random island on which this species is found),
    #or (Lud2) pick the island population weighted by island size;
    #and make this popn. go extinct, and so on until th threshold is met
    ###########################################################################
    
    #make an archipelago df with only unique species across the five islands
    #first add a column to each island's species' df with the species number/name
    islFull_Ex <- lapply(islFull_Ex, as.data.frame)
    iFE <- lapply(islFull_Ex, function(x){
      rn <- rownames(x)
      x$SN <- rn
      x
    })
    iFE <- lapply(iFE, as.data.frame)
    iFE2 <- dplyr::bind_rows(iFE) #then combine into one df
    iFE3 <- iFE2 %>% distinct(SN, .keep_all = TRUE) #keep only unique rows based on SN (i.e. unique species)
 
    repeat{
    
    #select a species from the archieplago based on bs
    sbs <- sample(1:nrow(iFE3), 1, prob = iFE3[, 1])
    ss <- iFE3$SN[sbs]
    #find which islands this species is found on
    wi <- vapply(islFull_Ex, function(x) ss %in% rownames(x), FUN.VALUE = logical(1))
    #if more than one island randomly select one; if just one then take that
    if (!any(wi)) stop("error in species extinction process (Ludwig)")
    if (length(which(wi)) > 1) {
      if (Ext_method == "lud2") {#select island by island area rather than area
        war <- ar[wi] 
        arP <- 1 - (war / (sum(war))) 
        si <- sample(which(wi), 1, prob = arP)
      } else { #or select it randomly
        si <- sample(which(wi), 1)
      }
    } else {
      si <- which(wi)
    }
    #if only one species left on island skip it, as need 1 species for SAr and Fd calcs etc
    nr <- nrow(islFull_Ex[[si]])
    if (nr == 1) next
    
    #make this chosen population go extinct on this particular island 
    wes <- which(rownames(islFull_Ex[[si]]) == ss)
    islFull_Ex[[si]] <- islFull_Ex[[si]][-wes,]
    #check if global extinction and if so remove from archip list and add to extinct list
    if (length(which(wi)) == 1){
      extinct <- c(extinct, ss)
      iFE3 <- iFE3[-sbs,]
    }
    
    if (length(extinct) == ceiling(nrow(species3) * th)) break
    }
    }#eo main if


    ###################################################################################
    
    ##stan method
    ##this iterates in turn through each island which means smaller islands
    ##lose a larger PROPORTION of species, as basically all island lose
    ##roughly same number of species, but smaller have lower starting point
    ##This is as to be expected as smaller island have > extinction rates

    if (Ext_method == "stan" || Ext_method == "prob"){
    
    if (Ext_method == "prob") arP <- 1 - (ar / (sum(ar))) # probability based on area
    
    
    j2 <- 0 
    repeat{
        for (i in 1:5){
            dum <- islFull_Ex[[i]]
            if (nrow(dum) == 1) next #always keep 1 sp on an island 
            ##select whether an extinction event occurs based on island size probability
            if (Ext_method == "prob") {
              rP <- rbinom(1, 1, arP[i])
              if (rP == 0) next #no extinction occurs
            }
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
    }#eo if

    dead2 <- rep("Present", length = nrow(species3))
    dead2[which(rownames(species3) %in% extinct)] <- "Extinct"

    ##make archipelago dataset and dendogram
    allsp_Ex <- unlist(lapply(islFull_Ex, rownames))
    allsp2_Ex <- unique(as.numeric(allsp_Ex))

    species3_Ex <- species[allsp2_Ex,]
    rownames(species3_Ex) <- allsp2_Ex
    
    #removing body size
    if (!bs_I)species3_Ex <- species3_Ex[,-1]
    
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
        if (bs_I == FALSE) bigI <- bigI[,-1]
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

tes <- Leo(plot_T = TRUE)


###########################################
#####functions to format the output of the main function and generate results table
####################################################

form_leo <- function(x = Leo2){
    
    ##all metrics for pre-colonisation data
    F0 <- apply(x, 2, function(y) unlist(y[[2]][1]))
    F0ue <- rbind(apply(F0, 1, mean), apply(F0, 1, plotrix::std.error))
  
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
    
    
    l1 <- list(F0ue, F1ue, F2ue, F3ue, F4ue, F5ue, F6ue)
    names(l1) <- c("SR", "FD", "FRIC_FDIS", "Chequer", "TaxoBeta", "FuncBeta", "SAR")
    
    ##all metrics for post-colonisation data
    F0b <- apply(x, 2, function(y) unlist(y[[8]][1]))
    F0bue <- rbind(apply(F0b, 1, mean), apply(F0b, 1, plotrix::std.error))
    
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
    
    
    l2 <- list(F0bue, F7ue, F8ue, F9ue, F10ue, F11ue, F12ue)
    names(l2) <- c("SR", "FD", "FRIC_FDIS", "Chequer", "TaxoBeta", "FuncBeta", "SAR")
    
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
    tt <- cbind(NA, tt)#add NA for first column (relating to species richness
    
    l3 <- list(l1, l2, tt)
    names(l3) <- c("Pre-humans", "Post-humans", "T_tests")
    
    return(l3)  
}


##save(Leo2, file = "Leo2.R")


Leo2 <- vector("list", length = 7)

d1 <- replicate(10, Leo(Ext_method = "stan"))
Leo2[[1]] <- form_leo(d1)

d2 <- replicate(10, Leo(Ext_method = "prob"))
Leo2[[2]] <- form_leo(d2)

d3 <- replicate(10, Leo(Ext_method = "lud1"))
Leo2[[3]] <- form_leo(d3)

d4 <- replicate(10, Leo(Ext_method = "lud2"))
Leo2[[4]] <- form_leo(d4)

d5 <- replicate(10, Leo(Ext_method = "stan", bs_I = TRUE))
Leo2[[5]] <- form_leo(d5)

d6 <- replicate(10, Leo(Ext_method = "stan", th = 0.3))
Leo2[[6]] <- form_leo(d6)

d7 <- replicate(10, Leo(Ext_method = "stan", th = 0.7))
Leo2[[7]] <- form_leo(d7)

anyNA(Leo2)

names(Leo2) <- c("Standard_extinction", "Probabilistic_extinction", "Ludwig_extinction1",
                 "Ludwig_extinction2", "Body_size_included", "th0.3", "th=0.7")

###function to create Table 1

ms_table <- function(z){

l1 <- z$`Pre-humans`$SR[1,] %>% round(0)
l2 <- z$`Post-humans`$SR[1,] %>% round(0)
t1 <- z$T_tests[,1]

l3 <- z$`Pre-humans`$FRIC_FDIS[1,1:6] %>% round(0)
l4 <- z$`Post-humans`$FRIC_FDIS[1,1:6] %>% round(0)
t2 <- z$T_tests[,3] %>% round(1)

l5 <- c(z$`Pre-humans`$TaxoBeta[1,], rep(NA,3)) %>% round(2)
l6 <- c(z$`Post-humans`$TaxoBeta[1,], rep(NA,3)) %>% round(2)
t3 <- z$T_tests[,5] %>% round(1)

l7 <- c(z$`Pre-humans`$FuncBeta[1,], rep(NA,3)) %>% round(2)
l8 <- c(z$`Post-humans`$FuncBeta[1,], rep(NA,3)) %>% round(2)
t4 <- z$T_tests[,6] %>% round(1)

l9 <- c(z$`Pre-humans`$SAR[3], rep(NA,5)) %>% round(2)
l10 <- c(z$`Post-humans`$SAR[3], rep(NA,5)) %>% round(2)
t5 <- z$T_tests[,7] %>% round(1)

l11 <- c(z$`Pre-humans`$Chequer[c(1,9)], rep(NA,4)) %>% round(2)
l12 <- c(z$`Post-humans`$Chequer[c(1,9)], rep(NA,4)) %>% round(2)
t6 <- z$T_tests[,4] %>% round(1)

p1 <- paste0(l1[1], " (", l1[2], ", ", l1[3], ", ", l1[4], ", ", l1[5], ", ", l1[6], ")")
p2 <- paste0(l2[1], " (", l2[2], ", ", l2[3], ", ", l2[4], ", ", l2[5], ", ", l2[6], ")")

p3 <- paste0(l3[1], " (", l3[2], ", ", l3[3], ", ", l3[4], ", ", l3[5], ", ", l3[6], ")")
p4 <- paste0(l4[1], " (", l4[2], ", ", l4[3], ", ", l4[4], ", ", l4[5], ", ", l4[6], ")")

p5 <- paste0(l5[1], " (", l5[2], ", ", l5[3], ")")
p6 <- paste0(l6[1], " (", l6[2], ", ", l6[3], ")")

p7 <- paste0(l7[1], " (", l7[2], ", ", l7[3], ")")
p8 <- paste0(l8[1], " (", l8[2], ", ", l8[3], ")")

p9 <- paste0(l9[1])
p10 <- paste0(l10[1])

p11 <- paste0(l11[1], " (SES = ", l11[2], ")")
p12 <- paste0(l12[1], " (SES = ", l12[2], ")")

rr <- rbind(c(p1,p2,t1), c(p3,p4,t2), c(p5,p6,t3), c(p7,p8,t4),
            c(p9,p10,t5),c(p11,p12,t6))

colnames(rr) <- c("Pre-human colonization", "Post-human colonization", 
                  "t-value", "P")
rownames(rr) <- c("Species richness", "Functional richness", "Taxonomic beta-diversity",
                  "Functional beta-diversity", "Slope of the SAR (z)", "C-score")

return(rr)
}


lapply(Leo2, ms_table) %>% write.csv()

