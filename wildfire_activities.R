library(INLA)
library(sp) 
library(fields)
library(evd)
library(stringr)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(spdep)

DATA = "data/"

poisson.sim.occ =function(x,Dat){
  # simulation of observing at least 1 fire + selection of corresponding events
  Dat[rpois(length(x), lambda = x)>0, ]
}

load("dat/DF_Aquitaine.Rdata") 
#replace with paste0(DATA,"DF_Aquitaine_sim.Rdata")

#SAFRAN coordinates
coord = read.csv(file=paste0(DATA,"Coord_L2E_dep.csv"),sep=";")
#Keep only summer fires
DF = DF[(DF$DOY >=154) & (DF$DOY<=314),]

#remove variables not used
DF$FWI = DF$FWIfl1
DF = DF[ ,c("PIX","XL2e","YL2e","DEP","YEAR","DOY","FWI","FA","NB0.1","BA")]

# To plot the maps with departements
library(sf)
dep = st_read("map_shp/DEPARTEMENT_CARTO.shp",quiet = TRUE)
departements = st_geometry(dep[(dep$INSEE_DEP %in% c("24", "33", "40", "47")),])
dep.coord = st_transform(dep, crs = 27572) #transform coordinates
dep.coord = dep.coord[(dep.coord$INSEE_DEP %in% c("24", "33", "40", "47")),] #extract departements

#subsampling
n_ss = 20000
q_FWI = 0.9
thresh_FWI = quantile(DF$FWI,probs=q_FWI)
ids2subsample = which(DF$NB0.1 == 0)
n = length(ids2subsample)
weights_FWI = ifelse(DF$FWI[ids2subsample]>thresh_FWI,q_FWI,1-q_FWI)
set.seed(1)
ids2keep = sample(ids2subsample, size = n_ss, prob = weights_FWI)
weights_ss = c(ifelse(DF$FWI[ids2keep]>thresh_FWI,((1-q_FWI)*n)/(n_ss/2),
                      (q_FWI*n)/(n_ss/2)), rep(1, sum(DF$NB0.1 > 0)))
ids2keep = c(ids2keep, which(DF$NB0.1 > 0))
DF_sub = DF[ids2keep,]
DF_sub$weight_ss = weights_ss
# dim(DF_sub)
rm(weights_FWI, ids2subsample, ids2keep)

## INLA settings
# dimension and shape parameter
d = 1
nu = 1
alpha = nu + d/2

# FWI mesh and SPDE
knots_FWI = c(5,10,15,20, 30)
bnd_FWI = c(0,max(DF_sub$FWI))
mesh_FWI = inla.mesh.1d(knots_FWI, bnd_FWI, degree=2, boundary=c("free", "neumann"))

B_FWI = inla.spde.make.A(mesh_FWI, DF_sub$FWI)
spde_FWI = inla.spde2.pcmatern(mesh_FWI, alpha = alpha, prior.sigma = c(1, NA), prior.range = c(5, NA), constr = TRUE)
idx_FWI = inla.spde.make.index("FWI", n.spde = spde_FWI$n.spde)

# DOY mesh and SPDE
nknots = 5
knots_DOY = quantile(DF_sub$DOY, probs = (1:nknots)/(nknots+1))
bnd_DOY = c(min(DF_sub$DOY) - 7, max(DF_sub$DOY) + 7)
mesh_DOY = inla.mesh.1d(knots_DOY, bnd_DOY, degree=2, boundary=c("free","free"))

spde_DOY = inla.spde2.pcmatern(mesh_DOY, alpha = alpha, prior.sigma = c(1, NA), prior.range = c(10, NA), constr = TRUE)
idx_DOY = inla.spde.make.index("DOY", n.spde = spde_DOY$n.spde)
B_DOY = inla.spde.make.A(mesh_DOY, DF_sub$DOY)

# FA mesh and SPDE
knots_FA = c(500,1000,2000,3000,4000,5000)
bnd_FA = c(0, 6400)
mesh_FA = inla.mesh.1d(knots_FA, bnd_FA, degree=2, boundary=c("free", "free"))

spde_FA = inla.spde2.pcmatern(mesh_FA, alpha = alpha, prior.sigma = c(1, NA), prior.range = c(3000, NA), constr = TRUE)
idx_FA = inla.spde.make.index("FA", n.spde = spde_FA$n.spde)
B_FA = inla.spde.make.A(mesh_FA, DF_sub$FA)

#spatial effect
# dimension and shape parameter
d = 2
nu = 1
alpha = nu + d/2

#definition of the mesh
PIXobs = unique(DF_sub$PIX)
coordL2e = coord[match(PIXobs,coord$N_SAFRAN),]
#rescaling coordinates (otherwise, INLA may be unstable)
coordinla = cbind(coordL2e$XL2e/10^6,coordL2e$YL2e/10^6)
#define boundaries for the domain of the SPDE model
bound1 = inla.nonconvex.hull(coordinla, convex= -0.03) #central boundary
bound2 = inla.nonconvex.hull(coordinla, convex= -0.25) #external boundary
#create triangulation
mesh_spatial = inla.mesh.2d(coordinla,boundary=list(bound1,bound2),cutoff=0.001, min.angle=21, max.edge= c(0.025,0.025))

spde_PIX = inla.spde2.pcmatern(mesh = mesh_spatial, alpha = alpha, prior.range = c(.05, .05), prior.sigma = c(1,.5), constr = T)
idx_PIX = inla.spde.make.index("spatial", n.spde = spde_PIX$n.spde)
B_PIX = inla.spde.make.A(mesh_spatial, 
                         loc = cbind(DF_sub$XL2e, DF_sub$YL2e)/10^6)

#for the Besag model on sizes we need to define the neighborhood structure 
dept = as_Spatial(dep[(dep$INSEE_DEP %in% c("24", "33", "40", "47")),])
spdep::nb2INLA("SIZEadjgraphDEP.txt",spdep::poly2nb(dept,queen=F,
                                                    row.names=dept$DEP))

## Estimation of the occurrence model

fixed_effects_df = data.frame(intercept = rep(1,nrow(DF_sub))) 
stack.obs = inla.stack(data = list(y = as.numeric(DF_sub$NB0.1)), 
                       A = list(B_DOY, B_FWI, B_PIX, B_FA,1,1 ), 
                       effects = list(idx_DOY, idx_FWI,  idx_PIX, idx_FA, 
                              list(YEAR=(DF_sub$YEAR-2005)),fixed_effects_df))

hyper_prec = list(prec = list(prior = 'pcprec', param = c(1, 0.5)))

formula = y ~ -1 +  intercept + f(YEAR,model="iid", hyper = hyper_prec) + f(DOY, model = spde_DOY)  +  f(FWI, model = spde_FWI)+  f(FA, model = spde_FA) +  f(spatial, model = spde_PIX)

fit=inla(formula,
         data = inla.stack.data(stack.obs),
         family = "poisson",
         quantiles = c(0.025, 0.5, 0.975),
         control.compute = list(return.marginals=T,config=TRUE,dic=TRUE,
                                waic=TRUE),
         control.predictor = list(compute=FALSE,A=inla.stack.A(stack.obs)),
         control.inla = list(strategy="adaptive",adaptive.max=1000),
         verbose = FALSE,
         E = DF_sub$weight_ss,
         num.threads = 2)


## Estimation of the size model

DFsize = DF[DF$BA > 0.1, ] 
B_FWI = inla.spde.make.A(mesh_FWI, DFsize$FWI)
B_FA = inla.spde.make.A(mesh_FA, DFsize$FA)
B_PIX = inla.spde.make.A(mesh_spatial, loc=cbind(DFsize$XL2e,DFsize$YL2e)/10^6)

fixed_effects_df = data.frame(intercept = rep(1,nrow(DFsize))) 
stack.obs.size = inla.stack(data = list(y = DFsize$BA), 
                            A = list(B_FWI, B_FA,1,1,1 ),
                            effects = list( idx_FWI,   idx_FA,list(DEP=match(DFsize$DEP,dept[[2]])),list(YEAR=(DFsize$YEAR-2005)),
                                            fixed_effects_df))

formula.size = y ~ -1 + intercept +f(FA, model = spde_FA) +  f(FWI, model = spde_FWI) + f(YEAR,model = "iid", hyper = hyper_prec) + f(DEP,model="iid")

fit.size=inla(formula.size,
         data = inla.stack.data(stack.obs.size),
         family = "gamma",
         quantiles = c(0.025, 0.5, 0.975),
         control.compute = list(return.marginals=T,config=TRUE,dic=TRUE,
                                waic=TRUE),
         control.predictor = list(compute=FALSE,A=inla.stack.A(stack.obs.size)),
         control.inla = list(strategy="adaptive",adaptive.max=1000),
         verbose = FALSE,
         num.threads = 2)

# Wildfire simulations from the posterior distribution for the observation period (2006-2020)
## Occurrence component

nsample = 100 # number of simulations
projDF = DF[order(DF$PIX),]
projDF$FWI = projDF$FWI
# constant extrapolation
projDF$FWI[projDF$FWI > mesh_FWI$interval[2]] = mesh_FWI$interval[2]
id_SAFRAN = unique(projDF[,c('PIX','XL2e','YL2e')])
coordinla_proj = cbind(id_SAFRAN$XL2e,id_SAFRAN$YL2e)/10^6
sample = inla.posterior.sample(n = nsample, result = fit)

x.comp=list()
fixedEffectsId=rownames(fit$summary.fixed)
eff = fixedEffectsId[1]
x.comp[[eff]] = as.vector(inla.posterior.sample.eval(eff, sample))
randomEffectsId=names(fit$summary.random)
for (eff in randomEffectsId) {
  x.comp[[eff]] = unname(inla.posterior.sample.eval(eff, sample))
}
randomEffectsId = randomEffectsId[randomEffectsId != "YEAR"]
for (eff in randomEffectsId){
  meshproj = get(paste0("mesh_", eff))
  if (str_detect(eff,"spatial")){
    meshproj = inla.mesh.projector(meshproj,loc = coordinla_proj)
  }
  else{meshproj = inla.mesh.projector(meshproj, loc = projDF[,eff])}
  x.comp[[eff]]=inla.mesh.project(meshproj, x.comp[[eff]])
}

# Then simulation of fire occurrences:
niter2 = 100; nsample2 = 1
projDF$Occur_sim = 1 # put them to 1 (only lines with actual fire will be   selected below)
Liste_occur = list()
all_lambdas = matrix(NA,nrow = dim(projDF)[1],ncol = niter2)
lambdamean =  rep(0, dim(projDF)[1])
for(k in 1:niter2){
  # OFFSET
  x.proj = 0 # if offset, not the case here
  # FIXED EFFECTS
  fixedEffectsId=rownames(fit$summary.fixed)
  if (length(fixedEffectsId)>0) {
    for (eff in fixedEffectsId) {
      x.proj = x.proj + rep(x.comp[[eff]][k], dim(projDF)[1])
    }
  }
  # RANDOM EFFECTS
  randomEffectsId=names(fit$summary.random)
  if (length(randomEffectsId)>0) {
    for (eff in randomEffectsId) {
      if (str_detect(eff,"spatial")){
        idspat = match(projDF$PIX,unique(projDF$PIX))
        x.proj = x.proj + x.comp[[eff]][idspat,k]
      }
      else{
        if (eff == "YEAR"){
          idvar = projDF[[eff]] - 2005
          x.proj = x.proj + x.comp[[eff]][idvar,k]
        }
        else{
        x.proj = x.proj + x.comp[[eff]][,k]
        }
      }
    }
  }
  lambda = exp(x.proj)
  lambdamean = lambdamean + lambda/niter2
  all_lambdas[,k] = lambda
  lambdas = matrix(rep(lambda,nsample2) ,nrow=length(lambda),ncol=nsample2)
  Liste_occur=cbind(Liste_occur,apply(lambdas,MARGIN = 2,poisson.sim.occ, Dat=projDF))
}

# Results:
nsim = length(Liste_occur)
nfireDFobstot=sum(projDF$NB0.1)
nfiresim=rep(0,nsim); for (i in 1:nsim){nfiresim[i]=sum(Liste_occur[[i]]$Occur_sim)}
nfireDFsimtot=mean(nfiresim)

projDF$lambda = lambdamean
l=list(Liste_occur=Liste_occur,nsim=nsim)


## Size component 

nsim = l$nsim
Proj = l$Liste_occur
nsample.size = 1
phi = fit.size$summary.hyperpar$mean[1]
for (i in 1:length(Proj)){
  projDF.size = Proj[[i]]
  id_SAFRAN = unique(projDF.size[,c('PIX','XL2e','YL2e')])
  coordinla_proj = cbind(id_SAFRAN$XL2e,id_SAFRAN$YL2e)/10^6
  sample.size = inla.posterior.sample(n = nsample.size, result = fit.size)
  x.comp=list()
  fixedEffectsId=rownames(fit.size$summary.fixed)
  eff = fixedEffectsId[1]
  x.comp[[eff]] = as.vector(inla.posterior.sample.eval(eff, sample.size))
  randomEffectsId=names(fit.size$summary.random)
  for (eff in randomEffectsId) {
    x.comp[[eff]] = unname(inla.posterior.sample.eval(eff, sample.size))
  }
  randomEffectsId = randomEffectsId[randomEffectsId != c("YEAR","DEP")]
  for (eff in randomEffectsId){
    meshproj = get(paste0("mesh_", eff))
    meshproj = inla.mesh.projector(meshproj, loc = projDF.size[,eff])
    x.comp[[eff]]=inla.mesh.project(meshproj, x.comp[[eff]])
  }
  
  # OFFSET
  x.proj = 0 # if offset, not the case here
  # FIXED EFFECTS
  if (length(fixedEffectsId)>0) {
    for (eff in fixedEffectsId) {
      x.proj = x.proj + rep(x.comp[[eff]],dim(projDF.size)[1])
    }
  }
  # RANDOM EFFECTS
  randomEffectsId=names(fit.size$summary.random)
  if (length(randomEffectsId)>0) {
    for (eff in randomEffectsId) {
      if (eff == "YEAR"){
        idvar = projDF.size[[eff]] - 2005
        x.proj = x.proj + x.comp[[eff]][idvar]
      } else if (eff == "DEP"){
        idvar = match(projDF.size[[eff]] ,dept[[2]])
        x.proj = x.proj + x.comp[[eff]][idvar]
      } else{
        x.proj = x.proj + x.comp[[eff]]
      }
    }
  }
  surf.sim = 0 * (1:nrow(projDF.size))
  mu = exp(x.proj)
  for (j in 1:nrow(projDF.size)){
    surf.sim[j] = rgamma(1,shape=phi,rate=phi/mu[j])
  }
  idx.withrepl=rep(1:nrow(projDF.size),projDF.size$Occur_sim) 
  Proj[[i]]$Surf_sim = aggregate(surf.sim,by=list(idx.withrepl),FUN=sum)[,-1]
  Proj[[i]]$N03_sim=aggregate(surf.sim>=0.3,by=list(idx.withrepl),FUN=sum)[,-1]
  Proj[[i]]$N05_sim=aggregate(surf.sim>=0.5,by=list(idx.withrepl),FUN=sum)[,-1]
  Proj[[i]]$N1_sim=aggregate(surf.sim>=1,by=list(idx.withrepl),FUN=sum)[,-1]
  Proj[[i]]$N10_sim=aggregate(surf.sim>=10,by=list(idx.withrepl),FUN=sum)[,-1]
}
Liste_occur = Proj
surfsim=c(1:nsim);for (i in 1:nsim){surfsim[i]=sum(Liste_occur[[i]]$Surf_sim)}
meansurfsim=mean(surfsim,na.rm=T)
surfsimobs = sum(DFsize$BA)
l.size = list(Liste_occur=Liste_occur,nsim=nsim)


# Future wildfire simulations derived from climate model projections 

model_names = c("IPSL_CM5A_WRF331F","MPI_ESM_RCA4","HadGEM_RCA4","CNRM_RCA4")
nsample = niter2 = 20
sample = inla.posterior.sample(n = nsample, result = fit)
FA_proj = unique(DF[, c("PIX", "FA","DEP")])

# Size projections
nsample.size = 1
phi = fit.size$summary.hyperpar$mean[1]
sample.size = inla.posterior.sample(n=nsample.size,result = fit.size)

#spde definitions for INLA
knots_proj =  c(2030, 2050, 2070, 2090)
bnd_proj = c(2020, 2100)
mesh_proj = inla.mesh.1d(knots_proj, bnd_proj, degree=2, boundary=c("dirichlet", "dirichlet"))
year_proj = 2021:2099
B_proj = inla.spde.make.A(mesh_proj, year_proj)
d = 1
nu = 1
alpha = nu + d/2
spde_proj = inla.spde2.pcmatern(mesh_proj, alpha = alpha, prior.sigma = c(100, 0.5), prior.range = c(30, NA))
idx_proj = inla.spde.make.index("proj", n.spde = spde_proj$n.spde)
n_years = length(year_proj)
fixed_effects_df = data.frame(intercept = rep(1, n_years),time = (year_proj-2020)/10)

beta = beta.size = NULL

for (model in 1:4){
  #FWI values
  name_clim_model = model_names[model]
  load(paste0("DATA","projections_FWI_",name_clim_model,".RData"))
  
  sim_counts_45 = sim_counts_85 = matrix(NA, nrow =nrow(fwi_rc45), ncol=nsample)
  
  # add FA values
  fwi_rc45$FA = FA_proj$FA[match(fwi_rc45$PIX, FA_proj$PIX)]
  fwi_rc85$FA = FA_proj$FA[match(fwi_rc85$PIX, FA_proj$PIX)]
  
  # add DEP values
  fwi_rc45$DEP = FA_proj$DEP[match(fwi_rc45$PIX, FA_proj$PIX)]
  fwi_rc85$DEP = FA_proj$DEP[match(fwi_rc85$PIX, FA_proj$PIX)]
  
  #Rename FWIfl1 for simplicity
  fwi_rc45$FWI = fwi_rc45$FWIfl1
  fwi_rc85$FWI = fwi_rc85$FWIfl1
  fwi_rc85$FWIfl1 = fwi_rc45$FWIfl1 = NULL
  
  # constant extrapolation
  fwi_rc45$FWI[fwi_rc45$FWI > mesh_FWI$interval[2]] = mesh_FWI$interval[2]
  fwi_rc85$FWI[fwi_rc85$FWI > mesh_FWI$interval[2]] = mesh_FWI$interval[2]
  
  id_SAFRAN = unique(fwi_rc45$PIX)
  coordinla_proj = cbind(coord$XL2e, coord$YL2e)[match(id_SAFRAN, coord$N_SAFRAN),]/10^6
  
  x.comp.45 = x.comp.85 =list()
  fixedEffectsId=rownames(fit$summary.fixed)
  eff = fixedEffectsId[1]
  x.comp.45[[eff]]=x.comp.85[[eff]] = as.vector(inla.posterior.sample.eval(eff, sample))
  randomEffectsId=names(fit$summary.random)
  for (eff in randomEffectsId) {
    x.comp.45[[eff]] = x.comp.85[[eff]] = unname(inla.posterior.sample.eval(eff, sample))
  }
  randomEffectsId = randomEffectsId[randomEffectsId != "YEAR"]
  for (eff in randomEffectsId){
    meshproj = get(paste0("mesh_", eff))
    if (str_detect(eff,"spatial")){
      meshproj.45 = meshproj.85 = inla.mesh.projector(meshproj,loc = coordinla_proj)
    }
    else{
      meshproj.45 = inla.mesh.projector(meshproj, loc = fwi_rc45[,eff])
      meshproj.85 = inla.mesh.projector(meshproj, loc = fwi_rc85[,eff])
    }
    x.comp.45[[eff]]=inla.mesh.project(meshproj.45, x.comp.45[[eff]])
    x.comp.85[[eff]]=inla.mesh.project(meshproj.85, x.comp.85[[eff]])
  }
  
  # Then simulation of fire occurrences:
  niter2 = nsample; nsample2 = 1
  idspat = match(fwi_rc45$PIX,unique(fwi_rc45$PIX))
  sim_counts_45 = sim_counts_85 = matrix(NA, nrow=nrow(fwi_rc45),ncol=niter2)
  for(k in 1:niter2){
    # OFFSET
    x.proj.45 = x.proj.85 = 0 # if offset, not the case here
    # FIXED EFFECTS
    fixedEffectsId=rownames(fit$summary.fixed)
    if (length(fixedEffectsId)>0) {
      for (eff in fixedEffectsId) {
        x.proj.45 = x.proj.45 + rep(x.comp.45[[eff]][k], dim(fwi_rc45)[1])
        x.proj.85 = x.proj.85 + rep(x.comp.85[[eff]][k], dim(fwi_rc85)[1])
      }
    }
    # RANDOM EFFECTS
    randomEffectsId=names(fit$summary.random)
    if (length(randomEffectsId)>0) {
      for (eff in randomEffectsId) {
        if (str_detect(eff,"spatial")){
          x.proj.45 = x.proj.45 + x.comp.45[[eff]][idspat,k]
          x.proj.85 = x.proj.85 + x.comp.85[[eff]][idspat,k]
        }
        else{
          if (eff == "YEAR"){
            #random sampling in the posterior distribution
            x.proj.45 = x.proj.45 + sample(x.comp.45[[eff]][,k], 1)
            x.proj.85 = x.proj.85 + sample(x.comp.85[[eff]][,k], 1)
          }
          else{
            x.proj.45 = x.proj.45 + x.comp.45[[eff]][,k]
            x.proj.85 = x.proj.85 + x.comp.85[[eff]][,k]
          }
        }
      }
    }
    lambda.45 = exp(x.proj.45)
    lambda.85 = exp(x.proj.85)
    sim_counts_45[,k] = rpois(length(lambda.45),lambda.45)
    sim_counts_85[,k] = rpois(length(lambda.85),lambda.85)
  }
  rm(lamda.45,lambda.85,x.proj.45,x.proj.85,x.comp.45,x.comp.85)
  gc()
  
  #Size projections
  sim_size_45_tot=sim_size_85_tot = matrix(0,nrow=dim(fwi_rc45)[1],ncol=niter2)
  for (k in 1:niter2){
    projDF.size.45 = fwi_rc45[sim_counts_45[,k]>0,]
    projDF.size.85 = fwi_rc85[sim_counts_85[,k]>0,]
    
    x.comp.45.size = x.comp.85.size =list()
    fixedEffectsId=rownames(fit.size$summary.fixed)
    eff = fixedEffectsId[1]
    x.comp.45.size[[eff]] = x.comp.85.size[[eff]] = as.vector(inla.posterior.sample.eval(eff, sample.size))
    randomEffectsId=names(fit.size$summary.random)
    for (eff in randomEffectsId) {
      x.comp.45.size[[eff]] = x.comp.85.size[[eff]] = unname(inla.posterior.sample.eval(eff, sample.size))
    }
    randomEffectsId = randomEffectsId[randomEffectsId != c("YEAR","DEP")]
    for (eff in randomEffectsId){
      meshproj = get(paste0("mesh_", eff))
      meshproj.45 = inla.mesh.projector(meshproj, loc = projDF.size.45[,eff])
      meshproj.85 = inla.mesh.projector(meshproj, loc = projDF.size.85[,eff])
      x.comp.45.size[[eff]]=inla.mesh.project(meshproj.45, x.comp.45.size[[eff]])
      x.comp.85.size[[eff]]=inla.mesh.project(meshproj.85, x.comp.85.size[[eff]]) 
    }
    x.proj.45.size = x.proj.85.size = 0 
    # FIXED EFFECTS
    if (length(fixedEffectsId)>0) {
      for (eff in fixedEffectsId) {
        x.proj.45.size = x.proj.45.size + rep(x.comp.45.size[[eff]],dim(projDF.size.45)[1])
        x.proj.85.size = x.proj.85.size + rep(x.comp.85.size[[eff]],dim(projDF.size.85)[1])
      }
    }
    # RANDOM EFFECTS
    randomEffectsId=names(fit.size$summary.random)
    if (length(randomEffectsId)>0) {
      for (eff in randomEffectsId) {
        if (str_detect(eff,"DEP")){
          idspat.45 = match(projDF.size.45[[eff]] ,dept[[2]]) #match(projDF.size.45$PIX,unique(projDF.size.45$PIX))
          x.proj.45.size = x.proj.45.size + x.comp.45.size[[eff]][idspat.45]
          idspat.85 = match(projDF.size.85[[eff]] ,dept[[2]]) #match(projDF.size.85$PIX,unique(projDF.size.85$PIX))
          x.proj.85.size = x.proj.85.size + x.comp.85.size[[eff]][idspat.85]
        }
        else{
          if (eff == "YEAR"){
            #random sampling in the posterior distribution
            x.proj.45.size = x.proj.45.size + sample(x.comp.45.size[[eff]], 1)
            x.proj.85.size = x.proj.85.size + sample(x.comp.85.size[[eff]], 1)
          }
          else{
            x.proj.45.size = x.proj.45.size + x.comp.45.size[[eff]]
            x.proj.85.size = x.proj.85.size + x.comp.85.size[[eff]]
          }
        }
      }
    }
    mu.45 = as.vector(exp(x.proj.45.size))
    mu.85 = as.vector(exp(x.proj.85.size))
    sim_size_45 = rgamma(n=dim(projDF.size.45)[1],shape = phi,rate=as.vector(phi/mu.45))
    sim_size_85 = rgamma(n=dim(projDF.size.85)[1],shape = phi,rate=as.vector(phi/mu.85))
    
    idx.45= which(sim_counts_45[,k]>0)#rep(1:dim(fwi_rc45)[1],sim_counts_45[,k])
    idx.85=which(sim_counts_85[,k]>0)#rep(1:dim(projDF.size.85)[1],sim_counts_85[,k])
    sim_size_45_tot[idx.45,k] = sim_size_45
    sim_size_85_tot[idx.85,k] = sim_size_85
  }
  rm(idx.45,idx.85,mu.45,mu.85,x.proj.45.size,x.proj.85.size,sim_size_45,sim_size_85,x.comp.45.size,x.comp.85.size,projDF.size.85,projDF.size.45)
  gc()
  
  # Map of mean between 2070-2100
  sim_counts_45_map =as.data.frame(sim_counts_45)
  sim_counts_45_map$YEAR = fwi_rc45$YEAR
  sim_counts_45_map$PIX = fwi_rc45$PIX
  sim_counts_45_map = sim_counts_45_map[sim_counts_45_map$YEAR >= 2070,]
  sim_counts_45_map$tot = apply(sim_counts_45_map[,c(1,nsample)],1,FUN=mean)
  map.mean.45 = aggregate(.~PIX,sim_counts_45_map[c("PIX","tot")],FUN = sum)
  idx_safran =  match(map.mean.45$PIX,coord$N_SAFRAN)
  map.mean.45 = map.mean.45$tot / 30

  sim_counts_85_map =as.data.frame(sim_counts_85)
  sim_counts_85_map$YEAR = fwi_rc85$YEAR
  sim_counts_85_map$PIX = fwi_rc85$PIX
  sim_counts_85_map = sim_counts_85_map[sim_counts_85_map$YEAR >= 2070,]
  sim_counts_85_map$tot = apply(sim_counts_85_map[,c(1,nsample)],1,FUN=mean)
  map.mean.85 = aggregate(.~PIX,sim_counts_85_map,FUN = sum)
  map.mean.85 = map.mean.85$tot / 30

  data = rbind(data.frame("NB"=map.mean.45,lab = 'RCP4.5'),data.frame("NB"=map.mean.85,lab="RCP8.5"))
  data$XL2e = rep(coord$XL2e[idx_safran],2)
  data$YL2e = rep(coord$YL2e[idx_safran],2)
  file_name = paste0(DATA,"map_proj_", model, ".RData")
  save(data,file = file_name)

  #Same but for sizes
  sim_size_45_map =as.data.frame(sim_size_45_tot)
  sim_size_45_map$YEAR = fwi_rc45$YEAR
  sim_size_45_map$PIX = fwi_rc45$PIX
  sim_size_45_map = sim_size_45_map[sim_size_45_map$YEAR >= 2070,]
  sim_size_45_map$BA = apply(sim_size_45_map[,c(1,nsample)],1,FUN=mean)
  map.mean.45 = aggregate(.~PIX,sim_size_45_map[,c('PIX',"BA")],FUN = sum)
  idx_safran =  match(map.mean.45$PIX,coord$N_SAFRAN)

  sim_size_85_map =as.data.frame(sim_size_85_tot)
  sim_size_85_map$YEAR = fwi_rc85$YEAR
  sim_size_85_map$PIX = fwi_rc85$PIX
  sim_size_85_map = sim_size_85_map[sim_size_85_map$YEAR >= 2070,]
  sim_size_85_map$BA = apply(sim_size_85_map[,c(1,nsample)],1,FUN=mean)
  map.mean.85 = aggregate(.~PIX,sim_size_85_map[,c('PIX',"BA")],FUN = sum)

  data = rbind(data.frame("BA"=map.mean.45$BA/30 ,lab = 'RCP4.5'),data.frame("BA"=map.mean.85$BA/30 ,lab="RCP8.5"))
  data$XL2e = rep(coord$XL2e[idx_safran],2)
  data$YL2e = rep(coord$YL2e[idx_safran],2)
  file_name = paste0(DATA,"map_proj_size_", model, ".RData")
  save(data,file = file_name)

  
  # Yearly mean + smooth for trend
  sim_counts_85 =as.data.frame(sim_counts_85)
  sim_counts_85$YEAR = fwi_rc85$YEAR
  yearly.mean.85 = aggregate(.~YEAR,sim_counts_85,FUN = sum)
  yearly.mean.85 = apply(yearly.mean.85[,-1],1,FUN=mean)

  sim_counts_45 =as.data.frame(sim_counts_45)
  sim_counts_45$YEAR = fwi_rc45$YEAR
  yearly.mean.45 = aggregate(.~YEAR,sim_counts_45,FUN = sum)
  yearly.mean.45 = apply(yearly.mean.45[,-1],1,FUN=mean)

  sim_size_85_tot =as.data.frame(sim_size_85_tot)
  sim_size_85_tot$YEAR = fwi_rc85$YEAR
  yearly.mean.85.size = aggregate(.~YEAR,sim_size_85_tot,FUN = sum)
  yearly.mean.85.size = apply(yearly.mean.85.size[,-1],1,FUN=mean)

  sim_size_45_tot =as.data.frame(sim_size_45_tot)
  sim_size_45_tot$YEAR = fwi_rc85$YEAR
  yearly.mean.45.size = aggregate(.~YEAR,sim_size_45_tot,FUN = sum)
  yearly.mean.45.size = apply(yearly.mean.45.size[,-1],1,FUN=mean)

  gc()
  
  #smoothing curves with INLA
  stack = inla.stack(data = list(y_45 = yearly.mean.45, y_85 = yearly.mean.85, y_45.size =  yearly.mean.45.size,y_85.size =  yearly.mean.85.size), A = list(B_proj, 1), effects = list(idx_proj, fixed_effects_df)) 
  hyper_prec=list(theta=list(initial = log(1), prior="pc.prec", 
                             param = c(1, 0.5)))
  my_formula = y_45 ~ -1 + intercept + time  + f(proj, model = spde_proj)
  fit_45 = inla(my_formula,
                data = inla.stack.data(stack),
                family = "lognormal",
                control.family = list(hyper = hyper_prec),
                control.predictor = list(compute=FALSE, A=inla.stack.A(stack)), verbose = FALSE, num.threads = 2)
  my_formula = y_85 ~ -1 +  intercept + time + f(proj, model = spde_proj)
  fit_85 = inla(my_formula,
                data = inla.stack.data(stack),
                family = "lognormal",
                control.family = list(hyper = hyper_prec),
                control.predictor = list(compute = FALSE, A = inla.stack.A(stack)),verbose = FALSE, num.threads = 2)
  my_formula = y_45.size ~ -1 + intercept + time  + f(proj, model = spde_proj)
  fit_45.size = inla(my_formula,
                     data = inla.stack.data(stack),
                     family = "lognormal",
                     control.family = list(hyper = hyper_prec),
                     control.predictor = list(compute = FALSE, A = inla.stack.A(stack)),verbose = FALSE, num.threads = 2)
  my_formula = y_85.size ~ -1 + intercept + time  + f(proj, model = spde_proj)
  fit_85.size = inla(my_formula,
                     data = inla.stack.data(stack),
                     family = "lognormal",
                     control.family = list(hyper = hyper_prec),
                     control.predictor = list(compute = FALSE, A = inla.stack.A(stack)),verbose = FALSE, num.threads = 2)
  
  #save slope parameters
  beta = rbind(beta,cbind(fit_45$summary.fixed[2,c(1,3,5)],
                          "lab"=paste(name_clim_model,"4.5")))
  beta = rbind(beta,cbind(fit_85$summary.fixed[2,c(1,3,5)],
                          "lab"=paste(name_clim_model,"8.5")))
  beta.size = rbind(beta.size,cbind(fit_45.size$summary.fixed[2,c(1,3,5)],
                                      "lab"=paste(name_clim_model,"4.5")))
  beta.size = rbind(beta.size,cbind(fit_85.size$summary.fixed[2,c(1,3,5)],
                                    "lab"=paste(name_clim_model,"8.5")))
  
  #plot the results
  eff = "proj"
  grid = inla.mesh.projector(mesh_proj, loc = year_proj)
  ngrid = length(grid$x)
  file_name = paste0("figures/plot_smooth_", model, ".png")
  png(file_name)
  ylim = range(yearly.mean.45, yearly.mean.85)
  curve_mean_45 = 100 * exp(fit_45$summary.linear.predictor$mean[1:ngrid])
  curve_upper_45 = 100*exp(fit_45$summary.linear.predictor$`0.975quant`[1:ngrid])
  curve_lower_45 = 100 * exp(fit_45$summary.linear.predictor$`0.025quant`[1:ngrid])
  curve_mean_85 = 100 * exp(fit_85$summary.linear.predictor$mean[1:ngrid])
  curve_upper_85 = 100 * exp(fit_85$summary.linear.predictor$`0.975quant`[1:ngrid])
  curve_lower_85 = 100 * exp(fit_85$summary.linear.predictor$`0.025quant`[1:ngrid])

  ggplot()+
    geom_point(aes_string(x=year_proj,y=yearly.mean.45),color = "gray50",size=2)+
    ggtitle(name_clim_model)+
    geom_point(aes_string(x=year_proj,y=yearly.mean.85),color='pink',size=2)+
    geom_line(aes_string(x=grid$x,y=curve_mean_45),color='gray20',size=1.5)+
    ylim(70,470)+
     geom_line(aes_string(x=grid$x,y=curve_lower_45),color='gray20',size=1.5,
               linetype = 2)+
    geom_line(aes_string(x=grid$x,y=curve_upper_45),color='gray20',size=1.5,
              linetype = 2)+
    geom_line(aes_string(x=grid$x,y=curve_mean_85),color='red',size=1.5)+
     geom_line(aes_string(x=grid$x,y=curve_lower_85),color='red',size=1.5,
               linetype = 2)+
    geom_line(aes_string(x=grid$x,y=curve_upper_85),color='red',size=1.5,
              linetype = 2)+
    theme_bw()+
    labs( y = 'Annual count', x = "Year")+
    theme(axis.text=element_text(size=19.3),
          axis.title=element_text(size=24),
          title=element_text(size=18))
  dev.off()

  file_name = paste0("figures/plot_smooth_size_", model, ".png")
  png(file_name)
    curve_mean_45 = 100 * exp(fit_45.size$summary.linear.predictor$mean[1:ngrid])
  curve_upper_45 = 100*exp(fit_45.size$summary.linear.predictor$`0.975quant`[1:ngrid])
  curve_lower_45 = 100 * exp(fit_45.size$summary.linear.predictor$`0.025quant`[1:ngrid])
  curve_mean_85 = 100 * exp(fit_85.size$summary.linear.predictor$mean[1:ngrid])
  curve_upper_85 = 100 * exp(fit_85.size$summary.linear.predictor$`0.975quant`[1:ngrid])
  curve_lower_85 = 100 * exp(fit_85.size$summary.linear.predictor$`0.025quant`[1:ngrid])
   
  ggplot()+
    geom_point(aes_string(x=year_proj,y=yearly.mean.45.size),color = "gray50",size=2)+
    ylim(0,2500)+
    ggtitle(name_clim_model)+
    geom_point(aes_string(x=year_proj,y=yearly.mean.85.size),color='pink',size=2)+
    geom_line(aes_string(x=grid$x,y=curve_mean_45),color='gray20',size=1.5)+
     geom_line(aes_string(x=grid$x,y=curve_lower_45),color='gray20',size=1.5,linetype = 2)+
    geom_line(aes_string(x=grid$x,y=curve_upper_45),color='gray20',size=1.5,linetype = 2)+
    geom_line(aes_string(x=grid$x,y=curve_mean_85),color='red',size=1.5)+
     geom_line(aes_string(x=grid$x,y=curve_lower_85),color='red',size=1.5,linetype = 2)+
    geom_line(aes_string(x=grid$x,y=curve_upper_85),color='red',size=1.5,linetype = 2)+
    theme_bw()+
    labs( y = 'Annual burnt area (ha)', x = "Year")+
    theme(axis.text=element_text(size=19.3),
          axis.title=element_text(size=24),
          title=element_text(size=18))
  dev.off()

  rm(fwi_rc45,fwi_rc85,id_SAFRAN,coordinla_proj,yearly.mean.85,yearly.mean.45,sim_counts_45,sim_counts_85,fit_85,fit_45,stack,yearly.mean.85.size,yearly.mean.45.size,sim_size_45_tot,sim_size_85_tot,fit_85.size,fit_45.size,idspat, meshproj.45,meshproj.85,sim_counts_45_map,sim_counts_85_map,sim_size_45_map,sim_size_85_map)
  gc()
}





