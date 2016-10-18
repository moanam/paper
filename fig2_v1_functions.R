####################################################################################################
# FUNCTION FOR FITTING AMPLIFICATION TRAJECTORIES TO LOGISTIC EQUATION

fit2logistic <- function(traj, times){
  # Fit trajectory to logistic function
  half_max = mean(c(min(traj), max(traj)))
  t_dat = tail(times[which(traj < half_max)], 2)
  d_dat = tail(traj[which(traj < half_max)], 2)
  r0 = log((d_dat[2] - min(traj))/(d_dat[1] - min(traj)))/(t_dat[2] - t_dat[1])
  p0 = exp(-r0 * t_dat[1])
  
  tryCatch({
    weights = c(rep(1, sum(traj < half_max)), rep(1, sum(traj >= half_max)))
    fit = nls(traj ~ (K*P0*exp(r*times))/(1 + P0*(exp(r*times) - 1)),
              start=c(K=max(traj), P0=p0, r=r0),
              control=c(maxiter=500), 
              weights=weights)
    
    coeffs = coef(fit)
    return(fit)
  }, error = function(err){print(err)}
  )
}

####################################################################################################
# FUNCTION FOR CALCULATING Ct FROM A LOGISTIC FIT

calc_Ct <- function(log_fit, threshold){
  coeffs = coef(log_fit)
  
  K = coeffs["K"]
  P0 = coeffs["P0"]
  r = coeffs["r"]
  
  ct = log((threshold*(1-P0))/(P0*(K-threshold)))/r
  
  return(ct)
}
####################################################################################################
# FUNCTION FOR CALCULATING STANDARD ERROR OF THE MEAN
stdErr <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x))
####################################################################################################
# FIT LOGISTIC FUNCTION TO POOLED PARTICLE COLONIZATION TRAJECTORIES (GLOBAL FIT)
# Note that this function fits the log-transformed data.  If you fit the data on a linear scale, the
# larger values are weighted much more heavily than the lower values.

globalFit2Logistic <- function(data.df){
  fit = nls(log10(Y) ~ log10((K*P0*exp(r*X))/(1 + P0*(exp(r*X) - 1))),
            data=data.df,
            start=c(K=1e5, P0=1.637e-04, r=1.80e-01),
            control=c(maxiter=500))
  
  return(fit)
}

####################################################################################################
# FUNCTION TO PLOT ERROR BARS ON PLOTS
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
####################################################################################################
# FUNCTION TO SMOOTH ABSOLUTE ABUNDANCE TRAJECTORIES

smooth_traj <- function(traj, times){
  library(caTools)
  
  traj.smooth <- NULL
  
  for (i in 1:length(times)){
    
    t = times[(i-1):(i+1)]
    y = traj[(i-1):(i+1)]
    traj.smooth = c(traj.smooth, predict(lm(y ~ t))[2])
  }
  
  return(traj.smooth)
}
####################################################################################################
# FUNCTION TO CALCULATE COM
calc_COM <- function(OTU_traj){
  times <- c(0,8,12,16,20,24,28,36,44,52,60,68,76,92,116,140)
  trajSum = sum(OTU_traj)
  
  COM = sum(OTU_traj * times)/trajSum
  return(COM)
}

####################################################################################################
# Convert read count data to relative abundances

count2rel <- function(dat){
    nRows = dim(dat)[1]
    nCols = dim(dat)[2]
    
    dat.rel = as.data.frame(matrix(nrow=nRows, ncol=nCols))
    row.names(dat.rel) = row.names(dat)
    colnames(dat.rel) = colnames(dat)
    
    for (i in 1:ncol(dat)){
        dat.rel[,i] = dat[,i]/sum(dat[,i])
    }
    
    return(dat.rel)
}
####################################################################################################
# Convert relative abundances to absolute abundances
rel2abs <- function(tab, abs.vals){
  nRows = dim(tab)[1]
  nCols = dim(tab)[2]
  tab.abs = as.data.frame(matrix(nrow=nRows, ncol=nCols))
  row.names(tab.abs) = row.names(tab)
  colnames(tab.abs) = colnames(tab)
  
  for (j in 1:nrow(tab)){
    otu.rel.vals = unlist(tab[j, ])
    otu.abs.vals = otu.rel.vals * abs.vals
    tab.abs[j, ] = otu.abs.vals
  }
  
  return(tab.abs)
}

####################################################################################################
# Normalize absolute abundances
normAbsAbund <- function(tab){
    nRows = dim(tab)[1]
    nCols = dim(tab)[2]
    tab.norm = as.data.frame(matrix(nrow=nRows, ncol=nCols))
    row.names(tab.norm) = row.names(tab)
    colnames(tab.norm) = colnames(tab)
    
    for (j in 1:nrow(tab)){
        otu.abs.vals = unlist(tab[j, ])
        otu.abs.norm.vals = otu.abs.vals/max(otu.abs.vals)
        tab.norm[j, ] = otu.abs.norm.vals
    }
    
    return(tab.norm)
}
####################################################################################################
# Calculate median, smoothed trajectories for all OTUs
medSmoothOTUs <- function(OTU.names, dat.M1.abs, dat.M2.abs, dat.M3.abs){
  dat.OTUs.abs.med.smooth.norm <- NULL
  dat.OTUs.abs.med.smooth <- NULL
  
  abs.vars = ls()[grepl(ls(), pattern="dat.M.*.abs")]
  
  for (i in 1:length(OTU.names)){
    OTU = OTU.names[i]
    
    trajs <- NULL
    for (j in 1:length(abs.vars)){
      temp.abs = get(abs.vars[j])
      temp.OTU.abs = unlist(temp.abs[OTU, ])
      trajs = rbind(trajs, temp.OTU.abs)
    }
    
    trajs.med = apply(trajs, 2, median)
    trajs.med.smooth = runmed(trajs.med, k=3, endrule="keep")
    trajs.med.smooth.norm = trajs.med.smooth / max(trajs.med.smooth)
    
    dat.OTUs.abs.med.smooth = rbind(dat.OTUs.abs.med.smooth, trajs.med.smooth)
    dat.OTUs.abs.med.smooth.norm = rbind(dat.OTUs.abs.med.smooth.norm, trajs.med.smooth.norm)
  }
  
  row.names(dat.OTUs.abs.med.smooth.norm) = OTU.names
  colnames(dat.OTUs.abs.med.smooth.norm) = times
  dat.OTUs.abs.med.smooth.norm = data.frame(dat.OTUs.abs.med.smooth.norm, check.names=FALSE)
 
  row.names(dat.OTUs.abs.med.smooth) = OTU.names
  colnames(dat.OTUs.abs.med.smooth) = times
  dat.OTUs.abs.med.smooth = data.frame(dat.OTUs.abs.med.smooth, check.names=FALSE)
  
  output <- list(dat.OTUs.abs.med.smooth, dat.OTUs.abs.med.smooth.norm)
  
  return(output)
}
####################################################################################################
diversityRarefaction = function(x, sampPerSize=10, sampleSizes=floor(seq(2,sum(x), length.out=100)), method='diversity', plot=FALSE){
  require(vegan)
  mat = matrix(nrow=length(sampleSizes),ncol=sampPerSize)
  i=1
  for (size in sampleSizes){
    samp = sapply(1:sampPerSize, FUN = function(a) sample(rep(1:length(x), times=x),size=size) )
    if (method == 'speciesRichness')
      mat[i,] = apply(samp, MARGIN=2, FUN=function(x) { length(unique(x)) } )
    if (method == 'diversity')
      mat[i,] = apply(samp, MARGIN=2, FUN=function(x) diversity(table(x)))
    
    i=i+1
  }
  if (plot){
    plot(rep(sampleSizes,each=sampPerSize), as.vector(t(mat)),pch=16,cex=.5)
    abline(h=diversity(x))
  }
  return( list(x=sampleSizes, y=rowMeans(mat)) )
}
####################################################################################################
smoothOTUs <- function(OTU.names, dat){
  dat.OTUs.abs.smooth.norm <- NULL
  dat.OTUs.abs.smooth <- NULL
    
  for (i in 1:length(OTU.names)){
    OTU = OTU.names[i]
    
    traj = dat[OTU, ]
    traj.smooth = runmed(traj, k=3)
    traj.smooth.norm = traj.smooth / max(traj.smooth)
    
    dat.OTUs.abs.smooth = rbind(dat.OTUs.abs.smooth, traj.smooth)
    dat.OTUs.abs.smooth.norm = rbind(dat.OTUs.abs.smooth.norm, traj.smooth.norm)
  }
  
  row.names(dat.OTUs.abs.smooth.norm) = OTU.names
  colnames(dat.OTUs.abs.smooth.norm) = times
  dat.OTUs.abs.smooth.norm = data.frame(dat.OTUs.abs.smooth.norm, check.names=FALSE)
 
  row.names(dat.OTUs.abs.smooth) = OTU.names
  colnames(dat.OTUs.abs.smooth) = times
  dat.OTUs.abs.smooth = data.frame(dat.OTUs.abs.smooth, check.names=FALSE)
  
  output <- list(dat.OTUs.abs.smooth, dat.OTUs.abs.smooth.norm)
  
  return(output)
}