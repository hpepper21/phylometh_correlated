getwd()
setwd("/Users/elipepper/Documents/GitHub")
##continuous
library(geiger)
library(ape)
#using examples from ape ?pic
tree.primates <- read.tree(text="((((Homo:0.21,Pongo:0.21):0.28,Macaca:0.49):0.13,Ateles:0.62):0.38,Galago:1.00);")
plot(tree.primates)
X <- c(4.09434, 3.61092, 2.37024, 2.02815, -1.46968)
Y <- c(4.74493, 3.33220, 3.36730, 2.89037, 2.30259)
names(X) <- names(Y) <- c("Homo", "Pongo", "Macaca", "Ateles", "Galago")
pic.X <- pic(X, tree.primates)
pic.Y <- pic(Y, tree.primates)
plot(pic.X[["rescaled.tree"]])

#positivize the contrasts
#create a regression through the origin
prim.mod<-lm(pic.Y~pic.X) #NULL mod
summary(prim.mod)
plot(prim.mod)
prim.mod2<-lm(pic.Y~pic.X+0) #mod w/ regression through origin
plot(pic.X,pic.Y)
abline(prim.mod2)

#the 0 term tells the lm() to fit the line through the origin
summary(prim.mod2)
anova(prim.mod,prim.mod2)

##discrete
require("corHMM")
?corHMM
data(primates)
ls()
print(primates)
require(phytools)

#review of discrete state mods
primates$trait[which(grepl("Hylobates",primates$trait[,1])),2]<-1

trait1<-primates$trait[,2]
names(trait1)<-primates$trait[,1]
plotSimmap(make.simmap(primates$tree, trait1), pts=FALSE, fsize=0.8)

rate.mat.er<-rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=2, model="ER")
print(rate.mat.er)
##what does this matrix mean?
#the func. rate.mat.maker generates and manipulates the index of the rate
#parameters to be optimized... it returns a rate matrix index;
#rate.cat specifies the number of rate categories in the hidden rates mod (HRM), which in this
#case rate.cat=1; hrm is a logical indicating whether the underlying mod is the HRM, in this case
#hrm=F (which is the default); ntraits specifies the number of traits in the data file if the 
#underlying model is not the HRM, in this case ntraits=1; nstates specifies the number of
#characters in the data file used in rayDISC, in this case nstates=2 (due to the binary-discrete dat);
#if the model is not HRM, model specifies the underlying model, in this case model="ER" - 
#"ER" stands for equal rates mod; the rate.mat.maker func. outputs the full index 
#of the rate parameters that are to be optimized,the intention 
#is that a user might want to see how the matrix is designed prior to an analysis
#and perhaps drops a few parameters beforehand due to some hypothesis (s)he might have;
#the resulting matrix can be plugged directly into corHMM, corDISC, or rayDISC

pp.er<-corHMM(primates$tree,primates$trait[,c(1,2)],rate.cat=1,rate.mat=rate.mat.er,node.states="marginal")
print(pp.er)
##what do these results mean?
#the er mod. is an equal rates mod.
#the corHMM func. takes a tree and a trait and estimates transition rates
#and ancestral states for a single binary character using the HRM

rate.mat.ard<-rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=2, model="ARD")
print(rate.mat.ard)
##what do these results mean?
#the ard mod is an all rates different mod

pp.ard<-corHMM(primates$tree,primates$trait[,c(1,2)],rate.cat=1,rate.mat=rate.mat.ard,node.states="marginal")
print(pp.ard)
##which mod is better?
#pp.er has a lower AIC, signifying that pp.er mod has a better fit

#now looking at multiple traits
#matrix w/ four states
rate.mat.er.4state<-rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=4, model="ER")
print(rate.mat.er.4state)

#converting the two binary traits into a single four character state
fourstate.trait<-rep(NA,Ntip(primates$tree))
for(i in sequence(Ntip(primates$tree))) {
  if(primates$trait[i,2]==0 && primates$trait[i,3]==0) {
    fourstate.trait[i]<-0
  }	
  if(primates$trait[i,2]==0 && primates$trait[i,3]==1) {
    fourstate.trait[i]<-1
  }	
  if(primates$trait[i,2]==1 && primates$trait[i,3]==0) {
    fourstate.trait[i]<-2
  }	
  if(primates$trait[i,2]==1 && primates$trait[i,3]==1) {
    fourstate.trait[i]<-3
  }	
}
fourstate.data<-data.frame(Genus_sp=primates$trait[,1], T1=fourstate.trait)

print(rayDISC(primates$tree, fourstate.data, ntraits=1, model="ER", node.states="marginal"))
print(rayDISC(primates$tree, fourstate.data, ntraits=1, rate.mat=rate.mat.er.4state, node.states="marginal", model="ARD"))
rate.mat.ard.4state<-rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=4, model="ARD")
print(rate.mat.ard.4state)

#making the equivalent of a GTR matrix:
rate.mat.gtr.4state<-rate.mat.ard.4state
rate.mat.gtr.4state<-rate.par.eq(rate.mat.gtr.4state, c(1,4))
rate.mat.gtr.4state<-rate.par.eq(rate.mat.gtr.4state, c(2,6))
rate.mat.gtr.4state<-rate.par.eq(rate.mat.gtr.4state, c(3,8))
rate.mat.gtr.4state<-rate.par.eq(rate.mat.gtr.4state, c(4,6))
rate.mat.gtr.4state<-rate.par.eq(rate.mat.gtr.4state, c(5,7))
rate.mat.gtr.4state<-rate.par.eq(rate.mat.gtr.4state, c(6,7))
print(rate.mat.gtr.4state)

print(rayDISC(primates$tree, fourstate.data, ntraits=1, rate.mat= rate.mat.gtr.4state, node.states="marginal", model="ARD"))

#making a mod like Pagel1994
print(rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=2, nstates=2, model="ARD"))
rate.mat.pag94<-rate.par.drop(rate.mat.ard.4state, drop.par=c(3,5,8,10))
print(rate.mat.pag94)

###Route 1
##Construct a model to test if state 1 can never be lost
#user-specified rate classes
mod <- corDISC(phy = primates$tree, data = primates$trait, ntraits = 2, rate.mat = rate.mat.pag94, model = "ARD", node.states = "marginal", diagn = FALSE)
print(mod)

##Experiment with the effects of frequencies at the root
rate.drop <- rate.par.drop(rate.mat.ard.4state, drop.par = c(1,2,3,5,6,8,9))
mod.ex <- rayDISC(primates$tree, fourstate.data, ntraits = 2, rate.mat = rate.drop, node.states = "marginal", model = "ARD", root.p = c(1,1,1,1))
mod.ex.1 <- rayDISC(primates$tree, fourstate.data, ntraits = 2, rate.mat = rate.drop, node.states = "marginal", model = "ARD", root.p = c(0,1,0,1))
mod.ex.2 <- rayDISC(primates$tree, fourstate.data, ntraits = 2, rate.mat = rate.drop, node.states = "marginal", model = "ARD", root.p = c(1,0,1,0))

##Create and use a model to see if transitions from 00 go to 11 only via 01
rate.drop <- rate.par.drop(rate.mat.ard.4state, drop.par = c(1,2,3,5,6,8,9))
mod2 <- rayDISC(primates$tree, fourstate.data, ntraits = 2, rate.mat = rate.drop, node.states = "marginal", model = "ARD", root.p = c(0.95,0.02,0.02,0.01))
print(mod2)

