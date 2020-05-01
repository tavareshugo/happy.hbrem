.packageName <- "happy.hbrem_2.2"


hbrem <- function( RX, HaploidInd, Ndip, Nstrain, Nind, Npost=2000, Nbin, Ry ) {

  brem <- .Call( "hbrem", RX, HaploidInd, Ndip, Nstrain, Nind, Npost, Nbin, Ry,  PACKAGE="happy.hbrem" )

  return(brem)
}


hbrem.true <- function( RX, HaploidInd, Ndip, Nstrain, Nind, Npost=2000, Nbin, Ry ) {
  
  brem.true <- .Call( "hbrem_true", RX, HaploidInd, Ndip, Nstrain, Nind, Npost, Nbin, Ry, PACKAGE="happy.hbrem" )

  return(brem.true)
}


hbrem.locus <- function(m, g, model, Ry, cc, HaploidInd, Ndip, Nstrain, Nind, Npost, Nbin) {

  if ( class(g) == "condensed.happy" ) {

    cum.mark <- 0
    chr <- 0
    while ( m > cum.mark) {
      chr = ( chr + 1 )
      cum.mark = ( cum.mark + length(g[[model]]$chr[[chr]]) )
    }

    cum.gen <- ( cum.mark - length(g[[model]]$chr[[chr]]) )
    chr.mark <- ( m - cum.gen )

    d <- g[[model]]$chr[[chr]][chr.mark][[1]][[1]]
    d.cc <- d[cc,]

  } else {
    d <- hdesign(g, m, model=model)
    d.cc <- d[cc,]
  }

  hb <- hbrem(RX=d.cc, HaploidInd=HaploidInd, Ndip=Ndip, Nstrain=Nstrain, Nind=Nind, Npost=Npost, Nbin=Nbin, Ry=Ry)

  reg.full.lm <- lm(Ry ~ d.cc)
  reg.null.lm <- lm(Ry ~ 1)
  Ftest <- anova(reg.null.lm, reg.full.lm)
  pval <- Ftest$"Pr(>F)"[2]
    
                                        #    cat( m, -log10(pval), "\n" )
  return( c(  -log10(pval), hb[[1]], hb[[2]], hb[[3]], hb[[4]], hb[[5]] ))
}


hbrem.perm.locus <- function(m, g, model, Ry, cc, HaploidInd, Ndip, Nstrain, Nind, Npost, Nbin) {

  if ( class(g) == "condensed.happy" ) {

    cum.mark <- 0
    chr <- 0
    while ( m > cum.mark) {
      chr = ( chr + 1 )
      cum.mark = ( cum.mark + length(g[[model]]$chr[[chr]]) )
    }

    cum.gen <- ( cum.mark - length(g[[model]]$chr[[chr]]) )
    chr.mark <- ( m - cum.gen )

    d <- g[[model]]$chr[[chr]][chr.mark][[1]][[1]]
    d.cc <- d[cc,]

  } else {
    d <- hdesign(g, m, model=model)
    d.cc <- d[cc,]
  }

  nrow.d <- length(d.cc[,1])
  ncol.d <- length(d.cc[1,])
  d.best <- matrix(rep(0, nrow.d*ncol.d), nrow.d, ncol.d)

  best.cc <- c()
  for ( i in 1:nrow.d ) {
    
    best <- which( d.cc[i,] == max(d.cc[i,]) )
    if ( length(best) == 1 ) {
      best.cc[i] = best
    } else {
      best.cc[i] = sample(best, 1)
    }

    d.best[ i, best.cc[i] ] = 1
  }

  hb <- hbrem.true(RX=d.best, HaploidInd=HaploidInd, Ndip=Ndip, Nstrain=Nstrain, Nind=Nind, Npost=Npost, Nbin=Nbin, Ry=Ry)

  reg.full.lm <- lm(Ry ~ d.cc)
  reg.null.lm <- lm(Ry ~ 1)
  Ftest <- anova(reg.null.lm, reg.full.lm)
  pval <- Ftest$"Pr(>F)"[2]

                                        #    cat( m, -log10(pval), "\n" )
  return( c(  -log10(pval), hb[[1]], hb[[2]], hb[[3]], hb[[4]], hb[[5]] ))
}


hbrem.merge.locus <- function(m, sdp, g, model, Ry, cc, HaploidInd, Ndip, Nstrain, Nind, Npost, Nbin) {

  d <- hdesign(g, m, model=model)
  d.cc <- d[cc,]

  if ( model == "additive" ) {

    all.0 <- numeric(length(d[,1]))
    all.1 <- numeric(length(d[,1]))
    for ( i in 1:Nstrain ) {
      if ( sdp[i] == 0 ) {
        all.0 = ( all.0 + d[,i] )
      } else if ( sdp[i] == 1 ) {
        all.1 = ( all.1 + d[,i] )
      } else {
        cat("sdp not 0 or 1\n")
      }
    }

    d.merge <- cbind(all.0, all.1)

  } else if ( model == "full" ) {

    sdp.matrix <- matrix(kronecker(sdp, sdp, paste, sep=""), nrow=Nstrain)
    sdp.vector <- c( diag(sdp.matrix), sdp.matrix[upper.tri(sdp.matrix, diag=FALSE)])
    sdp.full <- rep(1,36)
    sdp.full[sdp.vector == "00"] = 0
    sdp.full[sdp.vector == "11"] = 2

    dip.0 <- numeric(d[,1])
    dip.1 <- numeric(d[,1])
    dip.2 <- numeric(d[,1])
    for ( i in 1:Ndip ) {
      if ( sdp.full[i] == 0 ) {
        dip.0 = ( dip.0 + d[,i] )
      } else if ( sdp.full[i] == 1 ) {
        dip.1 = ( dip.1 + d[,i] )
      } else if ( sdp.full[i] == 2 ) {
        dip.2 = ( dip.2 + d[,i] )
      } else {
        cat("sdp not 0,1 or 2\n")
      }
    }

    d.merge <- cbind(dip.0,dip.1,dip.2)

  } else {
    cat("model not specified: must be one of additive, full\n")
  }


  hb <- hbrem(RX=d.merge, HaploidInd=HaploidInd, Ndip=3, Nstrain=2, Nind=Nind, Npost=Npost, Nbin=Nbin, Ry=Ry)

      reg.full.lm <- lm(Ry ~ d.merge)
      reg.null.lm <- lm(Ry ~ 1)
      Ftest <- anova(reg.null.lm, reg.full.lm)
      pval <- Ftest$"Pr(>F)"[2]

  #    cat( m, -log10(pval), "\n" )
      return( c(  -log10(pval), hb[[1]], hb[[2]], hb[[3]], hb[[4]], hb[[5]] ))

}


hbrem.region <- function(g, markers, Ry, cc, HaploidInd,  Npost, Nbin, Nperm=1000, thres.quick=c(0.5, 0.1, 0.05), thres.precise=c( 0.05, (1/length(markers)) ), thres.method="none", mc.cores=1) {
  
  if ( HaploidInd == 1 ) {
    model = "additive"
    Ndip = length(g$strains)
    Nstrain = Ndip
    Nind = length(Ry)
    marker.names <- g$additive$genome$marker[markers]
  }
  else if ( HaploidInd == 0 ) {
    model = "full"
    ns = length(g$strains)
    Ndip = ns*(ns+1)/2
    Nstrain = ns
    Nind = length(Ry)
    marker.names <- g$full$genome$marker[markers]
  }

  nmark <- length(markers)

  
  if ( mc.cores == 1 ) 
    res =  t(sapply ( markers, hbrem.locus, g, model, Ry, cc, HaploidInd, Ndip, Nstrain, Nind, Npost, Nbin) )
  else {
    res.list=mclapply ( markers, hbrem.locus, g, model, Ry, cc, HaploidInd, Ndip, Nstrain, Nind, Npost, Nbin, mc.cores=mc.cores)
    res = t(do.call( "cbind", res.list ) )
  }

  mark.pars.df <- data.frame(res[,1:46])
  names(mark.pars.df) = c( "F.logPval", "Hbar", "sd.Ni", "BIC.qtl", "BIC.null", "BF", "logBF", "DIC.qtl", "DIC.null", "DIC.diff", "pd.qtl", "pd.null", "mode.k", "ga", "gb", "mode.var", "med.k", "med.mu", "med.var", "mean.k", "mean.mu", "mean.var", "hpd.k.50.lower", "hpd.k.50.upper", "hpd.mu.50.lower", "hpd.mu.50.upper", "hpd.var.50.lower", "hpd.var.50.upper", "hpd.k.75.lower", "hpd.k.75.upper", "hpd.mu.75.lower", "hpd.mu.75.upper", "hpd.var.75.lower", "hpd.var.75.upper", "hpd.k.95.lower", "hpd.k.95.upper", "hpd.mu.95.lower", "hpd.mu.95.upper", "hpd.var.95.lower", "hpd.var.95.upper", "hpd.k.99.lower", "hpd.k.99.upper", "hpd.mu.99.lower", "hpd.mu.99.upper", "hpd.var.99.lower", "hpd.var.99.upper")
  offset = ncol(mark.pars.df)
  mark.pars.df$Name=as.character(marker.names)
  if ( class(g) == "happy.genome" ) {
    idx = match( marker.names, g[[model]]$genome$marker )
    mark.pars.df$Chr = g[[model]]$genome$chromosome[idx]
    mark.pars.df$Bp = g[[model]]$genome$bp[idx]

    bp2 = mark.pars.df$Bp[2:length(mark.pars.df$Bp)]
    bidx = which(bp2 < mark.pars.df$Bp[1:(length(mark.pars.df$Bp)-1)])-1

    mark.pars.df$CumBp = rep(0,nrow(mark.pars.df))
    mark.pars.df$CumBp[bidx+2] = mark.pars.df$Bp[bidx+1]
    mark.pars.df$CumBp = cumsum(mark.pars.df$CumBp) + mark.pars.df$Bp
         
  }

  cat("max logP = ", max(mark.pars.df$F.logPval), "\n")
  cat("max mode(k) = ", max(mark.pars.df$mode.k), "\n")


  hap.means.df <- data.frame( res[,(offset+1):(offset+Ndip)])
  hap.sdevs.df <- data.frame( res[,(offset+Ndip+1):(offset+2*Ndip)])
  hap.avNis.df <- data.frame( res[,(offset+2*Ndip+1):(offset+3*Ndip)])
  strain.means.df <- data.frame( res[,(offset+3*Ndip+1):(offset+3*Ndip+Nstrain)])

  if ( HaploidInd == 1 ) {
   names(hap.means.df) <- g$strains
   names(hap.sdevs.df) <- g$strains
   names(hap.avNis.df) <- g$strains
   names(strain.means.df) <- g$strains
  }
  else if ( HaploidInd == 0 ) {
   strain.names = g$strains
   num.strains = length(g$strains)
   diplotype.names <- matrix(kronecker(strain.names, strain.names, paste, sep="."), nrow=num.strains)
   names.full.symmetric <- c( diag(diplotype.names), diplotype.names[upper.tri(diplotype.names, diag=FALSE)])
   names(hap.means.df) <- names.full.symmetric
   names(hap.sdevs.df) <- names.full.symmetric
   names(hap.avNis.df) <- names.full.symmetric
   names(strain.means.df) <- g$strains
 }


  permuted.y=NULL
  if ( is.numeric(Nperm) & Nperm>0 ) {
    permuted.y = replicate( Nperm, sample(Ry) )
  }
  FlogP.thres=NULL
  modek.thres=NULL

  if ( thres.method == "precise" ) {

    FlogP.region = matrix( rep(0, nmark*Nperm), Nperm, nmark )
    modek.region = matrix( rep(0, nmark*Nperm), Nperm, nmark )

    for ( i in 1:Nperm ) {

      cat("perm ", i, "\n")
      
      if ( mc.cores == 1 )
        res =  t(sapply ( markers, hbrem.locus, g, model, permuted.y[,i], cc, HaploidInd, Ndip, Nstrain, Nind, Npost, Nbin) )
      else {
        res.list=mclapply ( markers, hbrem.locus, g, model, permuted.y[,i], cc, HaploidInd, Ndip, Nstrain, Nind, Npost, Nbin, mc.cores=mc.cores)
        res = t(do.call( "cbind", res.list ) )
      }

      FlogP.vector <- res[,1]
      modek.vector <- res[,13]
      for ( j in 1:nmark ) {
        FlogP.region[i,j] = FlogP.vector[j]
        modek.region[i,j] = modek.vector[j]
      }
      
    }

    for ( j in 1:nmark ) {
      FlogP.region[,j] <- sort.list(-FlogP.region[,j])
      modek.region[,j] <- sort.list(-modek.region[,j])
    }

    FlogP.thres <- FlogP.region[Nperm*thres.precise,]
    modek.thres <- modek.region[Nperm*thres.precise,]
    
  } else if ( thres.method == "quick" ) {

    FlogP.regionwide.distn <- numeric(Nperm)
    modek.regionwide.distn <- numeric(Nperm)
  
    for ( i in 1:Nperm ) {

      cat("perm ", i, "\n")
      
      if ( mc.cores == 1 )
        res =  t(sapply ( markers, hbrem.perm.locus, g, model, permuted.y[,i], cc, HaploidInd, Ndip, Nstrain, Nind, Npost, Nbin) )
      else {
        res.list=mclapply ( markers, hbrem.perm.locus, g, model, permuted.y[,i], cc, HaploidInd, Ndip, Nstrain, Nind, Npost, Nbin, mc.cores=mc.cores)
        res = t(do.call( "cbind", res.list ) )
      }

      FlogP.vector <- res[,1]
      modek.vector <- res[,12]
      max.FlogP <- max(FlogP.vector)
      max.modek <- max(modek.vector)

      FlogP.regionwide.distn[i] = max.FlogP
      modek.regionwide.distn[i] = max.modek
    }

    FlogP.sort <- FlogP.regionwide.distn[ sort.list(-FlogP.regionwide.distn) ]
    modek.sort <- modek.regionwide.distn[ sort.list(-modek.regionwide.distn) ]
    
    FlogP.thres <- FlogP.sort[Nperm*thres.quick]
    modek.thres <- modek.sort[Nperm*thres.quick]
  }
  
 hbrem.region.list <- list(Summary.Parameters=mark.pars.df, Diplo.Means=hap.means.df, Diplo.StDevs=hap.sdevs.df, Diplo.ExpCounts=hap.avNis.df, Strain.Means=strain.means.df, F.logP.thres=FlogP.thres, Mode.k.thres=modek.thres)

  return(hbrem.region.list)
}
