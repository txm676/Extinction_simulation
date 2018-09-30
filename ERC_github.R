#############################################################################
######ERC SPECIES LOSS SIMULATION#########################################
#######################################################################

#aim: to populate islands with realistic species set, and
#realistic impact of humans

#is neutral in regards to fitness, no fitness differences between species;
#but traits do differ

#assumes trait conservatism in speciation. Dispersal ability assumed to
#decrease with speciation

#species either disperses or speciates (0.9 - 0.1 prob); for the former
#there is 0.5 chance of inter-island colonisation, and 0.5 from pool

library(sars)
library(CommEcol)
library(vegan)
library(ape)
library(picante)
library(phytools)
library(FD)

##simulate X species with three traits

xx <- matrix(nrow = 100, ncol = 3)
colnames(xx) <- c("BS", "D", "Beak")

for (i in 1:nrow(xx)){
  xx[i, 1] <- rgamma(1, 1)#body size
  xx[i, 2] <- (rpois(1, 2.5) + 1)#dispersal
  xx[i, 3] <- runif(1, 1, 8)#beak shape
}#eo for

##create islands
isl <- vector("list", length = 5)#list to put island species in
ar <- c(0.1, 2, 4, 10, 50)#island areas
k <-  vapply(ar, function(x) 20 * x ^ 0.25, FUN.VALUE = numeric(1))  #z = 0.25, c = 20
k <- ceiling(k)#carrying capacity

##run island colonisation process
repeat{
  for (i in 1:length(isl)){
    if (length(isl[[i]]) == k[i]) next()
    dum <- isl[[i]]#vector of species on an island
    #select whether immigrant (1) or speciation event (0)
    v <- rbinom(1, 1, 0.9)
    if (length(dum) == 0 || v == 1){ #if immigration
      #select if immigrant from pool or other island (50% probability)
      i2 <- rbinom(1, 1, 0.5)
      if (length(dum) == 0 || i2 == 1){ #if immigration from pool
        d <- sample(1:100, 1, prob = xx[1:100,2])#pick a species weighted by dispersal ability (prob does not have to be 0-1)
        if (d %in% dum) next #if already on island skip
        dum <- c(dum, d)
      } else {#else dispersal from a different island
        wi <- 1:5
        wi2 <- wi[-i]#vector of island numbers excl. i
        si <- sample(wi2, 1, prob = ar[-i])#pick one island at random (weighted by area)
        ins <- isl[[si]] #species numbers on the si'th island
        xxi <- xx[ins, , drop = FALSE]#subset the xx matrix to just have species on island si
        #sample a species from this subset weighted by dispersal ability
        di <- sample(1:nrow(xxi), 1, prob = xxi[, 2])
        d <- ins[di]#work out which species this sample refers to
        if (d %in% dum) next #if already on island skip
        dum <- c(dum, d)
      }
    } else { #if speciation
      sp <- sample(dum, 1)
      sptr <- xx[sp, ]#pick a random species from whole pool of species (including speciated ones)
      #ensure no traits are negative or 0 values
      repeat{
        sptr[c(1, 3)] <- jitter(sptr[c(1, 3)])#for body size and beak add ranom noise
        if (all(sptr[c(1, 3)] > 0)) break
      }
      sptr[2] <- ifelse(sptr[2] == 1, sptr[2], (sptr[2] - 1))#for dispersal - reduce unless already 1
      if (any(sptr == 0)) next
      xx <- rbind(xx, sptr)
      dum <- c(dum, nrow(xx))#add new speciated species to vector of names
    }#eo else
    
    isl[[i]] <- dum #the species names on the ith island
  }#eo for
  #cat(length(isl[[5]]), "\n")
  if (all(vapply(isl, length, FUN.VALUE = numeric(1)) == k)) break
}#eo repeat


if(!all(vapply(isl, length, FUN.VALUE = numeric(1)) == k)) stop("richness does not equal k")


##create list with the full trait matrix for each island
islFull <- lapply(isl, function(x){
  xx2 <- xx[x, ]
  rownames(xx2) <- x
  xx2
})

##make dendogram for each island as tree object
dendo <- function(x){
  r4 <- apply(x, 2, scale)
  rownames(r4) <- rownames(x)
  a4 <- vegan::vegdist(r4,method = "euclidean")
  a5 <- hclust(a4, method = "average")
  a6 <- ape::as.phylo(a5)
  a6$tip.label <- rownames(r4) 
  a6
}

dendz <- lapply(islFull, function(x){dendo(x)})

#make archipelago dataset and dendogram
allsp <- unlist(lapply(islFull, rownames))
allsp2 <- unique(as.numeric(allsp))

xx3 <- xx[allsp2,]
rownames(xx3) <- allsp2
arcDen <- dendo(xx3)

#calculate PA matrix for archi and each island; calculate PD
CM <- matrix(0, nrow = nrow(xx3), ncol = 6)
rownames(CM) <- rownames(xx3)
colnames(CM) <- c("A", 1:5)
CM[ ,1] <- 1
for (i in 1:5){
  dd <- isl[[i]]
  mat <- which(rownames(CM) %in% dd)
  CM[mat, (i + 1)] <- 1
}

CM <- t(CM)#picante has species as columns
picante::pd(CM, arcDen)#get FD of each island
#Functional richness
FF <- FD::dbFD(x = xx3, a = CM, w.abun = FALSE, stand.x = TRUE)


#sar z value
sarDf <- data.frame("A" = ar, "S" = rowSums(CM)[2:6])
sars::sar_power(sarDf)

#tree NODF
CM2 <- CM[2:6, ]#remove first row as this is the archipelago
any(colSums(CM2) == 0)
nodf_a <- CommEcol::treeNodf(CM2, col.tree = arcDen, order.rows = TRUE, order.cols = TRUE)$rows

##REPEAT ABOVE 100 times and take average FD and treeNODF values

###########################################################
#######HUMAN IMPACT##########################
##############################################

##EXTINCTION PROCESS

extinct <- c()

islFull_Ex <- islFull

#this iterates in turn through each island which means smaller islands
#lose a larger PROPORTION of species, as basically all island lose
#roughly same number of species, but smaller have lower starting point
#This is as to be expected as smaller island have > extinction rates

repeat{
  for (i in 1:5){
    dum <- islFull_Ex[[i]]
    if (nrow(dum) == 1) next #always keep 1 sp on an island
    #sample random sp from ith island weighted by body size
    sdum <- sample(1:nrow(dum), 1, prob = dum[, 1])
    samSp <- dum[sdum, , drop = FALSE]
    uid <- as.numeric(rownames(samSp))
    #check remaining sp to see if is SIE or not
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
    if (length(extinct) == ceiling(nrow(xx3) / 2)) break
    print(i)
  }#eo for
  if (length(extinct) == ceiling(nrow(xx3) / 2)) break
}#eo repeat

dead2 <- rep("Present", length = nrow(xx3))
dead2[which(rownames(xx3) %in% extinct)] <- "Extinct"

##recreate lists excluding extinct species

#make archipelago dataset and dendogram
allsp_Ex <- unlist(lapply(islFull_Ex, rownames))
allsp2_Ex <- unique(as.numeric(allsp_Ex))

xx3_Ex <- xx[allsp2_Ex,]
rownames(xx3_Ex) <- allsp2_Ex
arcDen_Ex <- dendo(xx3_Ex)

#calculate PA matrix for archi and each island; calculate PD
CM_Ex <- matrix(0, nrow = nrow(xx3_Ex), ncol = 6)
rownames(CM_Ex) <- rownames(xx3_Ex)
colnames(CM_Ex) <- c("A", 1:5)
CM_Ex[ ,1] <- 1
for (i in 1:5){
  dd <- rownames(islFull_Ex[[i]])
  mat <- which(rownames(CM_Ex) %in% dd)
  CM_Ex[mat, (i + 1)] <- 1
}

CM_Ex <- t(CM_Ex)#picante has species as columns
picante::pd(CM_Ex, arcDen_Ex)#get FD of each island
#Functional richness
FF_Ex <- FD::dbFD(x = xx3_Ex , a = CM_Ex , w.abun = FALSE, stand.x = TRUE)

#sar z value
sarDf_Ex <- data.frame("A" = ar, "S" = rowSums(CM_Ex)[2:6])
sars::sar_power(sarDf_Ex)

#tree NODF
CM2_Ex <- CM_Ex[2:6, ]#remove first row as this is the archipelago
any(colSums(CM2_Ex) == 0)
nodf_Ex <- CommEcol::treeNodf(CM2_Ex, col.tree = arcDen_Ex, order.rows = TRUE, order.cols = TRUE)$rows

##plot the dendogram with extinct species and extant species

fmode <- as.factor(setNames(dead2,rownames(xx3)))
phytools::dotTree(arcDen, fmode, colors=setNames(c("blue","red"),
                                                 c("Present", "Extinct")))
