#Lehmann et al. 2018
#Complex responses of global insect pests to climate change
#================================================================================
#Phylogenetically corrected regression analyses
#================================================================================

rm(list=ls())
#packages required:

library(picante)
library(adephylo)
library(ade4)
library(geiger)
library(ape)
library(rotl)
library(caper)
library(MCMCglmm, quiet=TRUE)

#----------------------------------------------------------------------------------------------
# Build phylogenetic tree of the 31 pest insects
#----------------------------------------------------------------------------------------------

taxa <- tnrs_match_names(names = c("Christoneura fumiferana",
                                   "Lymantria dispar",
                                   "Operopthera brumata",
                                   "Epirrita autumnata",
                                   "Thaumetopoea pityocampa",
                                   "Leptinotarsa decemlineata",
                                   "Locusta migratoria",
                                   "Meligethes aeneus",
                                   "Plutella xylostella",
                                   "Rhopalosiphum padi",
                                   "Diuraphis noxia",
                                   "Bemisia tabaci",
                                   "Nezara viridula",
                                   "Chilo suppressalis",
                                   "Ostrinia nubialis",
                                   "Helicoverpa armigera",
                                   "Dendroctonus ponderosae",
                                   "Dendroctonus frontalis",
                                   "Ips typographus",
                                   "Bactrocera oleae",
                                   "Cydia pomonella",
                                   "Hypothenemus hampei",
                                   "Diabrotica virgifera",
                                   "Popillia japonica",
                                   "Eldana saccharina",
                                   "Leucoptera coffeella",
                                   "Marmara gulosa",
                                   "Phyllocnistis citrella",
                                   "Chilo partellus",
                                   "Myzus persicae"))

mytree <- tol_induced_subtree(ott_ids = ott_id(taxa))

plot(mytree, cex = .8, label.offset = .1, no.margin = TRUE)

write.tree(mytree, file="filename_tree.phy", digits=3, tree.names = T)

mytree <- read.tree("filename_tree.phy")
mydata <- read.table("filename_data.txt",header=T,sep="\t")

mysorteddata <- mydata[mytree$tip.label, ]

plot(mytree, show.tip.label=T)
attach(mydata)

mytree$edge.length
mytree <- compute.brlen(mytree, 1)

mytree$node.label
mytree<-makeLabel(mytree)
mytree$node.label

mytree$node.label <- NULL

set.seed(123)

pr<-list(R=list(V=1,nu=0.002),
         G=list(G1=list(V=1,nu=0.002))
)

model <- MCMCglmm(TOPT_C~1,
                  pedigree=mytree,
                  prior=pr,
                  data=mydata,
                  verbose=FALSE)

mytree$tip.label <- strip_ott_ids(mytree$tip.label, remove_underscores=TRUE)
mytree$tip.label %in% mydata$ott_name
mydata$Species[1]

TOPT.BM<-fitContinuous(mytree, mydata$TOPT_C, model="BM")

my.disp.tree=mytree
my.disp.tree$tip.label=format(mydata$TOPT_C,digits=2)
plot(my.disp.tree) 

#----------------------------------------------------------------------------------------------
# Test for associations between thermal traits and ambient temperature
#----------------------------------------------------------------------------------------------

# using the pgls function from "caper" to test for effects between 
# TOPT_C, TSUIT_past, TSUIT_future, and latitude)

data <- comparative.data(mytree, mydata, "Species", vcv=TRUE)
mod1 <- pgls(TOPT_C ~ Latitude, data, lambda = 'ML')
print(mod1)
summary(mod1)

plot(TOPT_C ~ Latitude, mydata)

abline(mod1)
profile_lambda=pgls.profile(mod1, which="lambda") # vary lambda
plot(profile_lambda)


data <- comparative.data(mytree, mydata, "Species", vcv=TRUE)
mod2 <- pgls(TSUIT_past ~ Latitude, data, lambda = 'ML')
print(mod2)
summary(mod2)

plot(TSUIT_past ~ Latitude, mydata)

abline(mod2)
profile_lambda=pgls.profile(mod2, which="lambda") # vary lambda
plot(profile_lambda)

data <- comparative.data(mytree, mydata, "Species", vcv=TRUE)
mod3 <- pgls(TSUIT_present ~ Latitude, data, lambda = 'ML')
print(mod3)
summary(mod3)

plot(TSUIT_present ~ Latitude, mydata)

abline(mod3)
profile_lambda=pgls.profile(mod3, which="lambda") # vary lambda
plot(profile_lambda)

data <- comparative.data(mytree, mydata, "Species", vcv=TRUE)
mod4 <- pgls(TSUIT_near_future ~ Latitude, data, lambda = 'ML')
print(mod4)
summary(mod4)

plot(TSUIT_near_future ~ Latitude, mydata)

abline(mod4)
profile_lambda=pgls.profile(mod4, which="lambda") # vary lambda
plot(profile_lambda)

data <- comparative.data(mytree, mydata, "Species", vcv=TRUE)
mod5 <- pgls(TSUIT_future ~ Latitude, data, lambda = 'ML')
print(mod5)
summary(mod5)

plot(TSUIT_future ~ Latitude, mydata)

abline(mod5)
profile_lambda=pgls.profile(mod5, which="lambda") # vary lambda
plot(profile_lambda)

#----------------------------------------------------------------------------------------------
#END_SCRIPT