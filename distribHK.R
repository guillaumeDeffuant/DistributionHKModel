#Author: Guillaume Deffuant
# May 2024.
# Hegselmann and Krause Model in disrtribution


#------------------------------------------
# Function running the model on distribution
# epsilon: confidence bound (in (0,1))
# ne : discretization of epsilon
# nstepM : maximum number of steps for running the model
# odd: if odd = 0, discretization of grid can be even or odd, if odd = 1, discretization  of grid is odd, if odd = 2, discretizartion is even
# rd = 0 no noise, if rd > 0, factor of noise (small receomended e.g. 0.01)
runDistModelHK = function(epsilon, ne, nstepsM, odd = 0, visu = TRUE, save = TRUE, rd = 0){
  
  #print(paste("save", save))
  if (visu) save = TRUE
  nGrid = round(ne / epsilon)
  if (odd == 1){
    if ((nGrid %% 2) == 0) nGrid = nGrid+1
  }
  if (odd == 2) {
    if ((nGrid %% 2) == 1) nGrid = nGrid+1
  }
  ops = 0.5/nGrid + ((0:(nGrid-1))/nGrid)
  #states[(1:N0)*ng, 1] = rep(1,N0)
  states = rep(1,nGrid)
  nh = trunc(nGrid / 2)
  #nrd = round(ne/2)
  nrd = ne
  if (rd>0) {
    if ((nGrid %%2) == 1) states[(nh +1- nrd):(nh+1+nrd)] = states[(nh+1 - nrd):(nh+1+nrd)] + runif(2*nrd+1, -rd, rd)
    else states[(nh+1 - nrd):(nh+nrd)] = states[(nh +1- nrd):(nh+nrd)] + runif(2*nrd, -rd, rd)
  }
  traj = list(states)
  nh = trunc(nGrid/2)
  par = nGrid %% 2
  changes = 1
  tt = 0
  # print(nGrid)
  while((tt < nstepsM) & (changes > 0)){
    tt = tt+1
    delt = rep(0, nGrid)
    changes = 0
    if (rd>0) nnt = nGrid
    else nnt = nh+par
    for(j in (1:nnt)){
      # for(j in (1:nGrid)){
      opj = states[j]
      if (opj >  0){
        km = max(1, j - ne)
        kM =  min(nGrid, j + ne)
        vv = sum(states[km:kM]*(km:kM)) / sum(states[km:kM])
        if (opj > 0) {
          sup0 = which(states[km:kM] > 0)
          if ((length(sup0) > 3) | (((sup0[length(sup0)]-sup0[1]) > (length(sup0)-1))) ) changes = changes+1
        }
        if (rd == 0) {if (j > nh) opj = opj / 2}
        tvv =trunc(vv)
        w = vv - tvv
        delt[tvv] = delt[tvv] + opj*(1-w)
        delt[tvv+1] = delt[tvv+1] + opj*w
      }
    }
    # if (par == 1) delt[nh+1] = delt[nh+1]/2
    if (rd == 0) delt = delt+ rev(delt)
    #  diff = states - delt
    #  changes = sum(abs(diff))
    states = delt
    if (save) traj[[tt+1]] = states
    cat("\r",tt,"      ", changes,"       ", nGrid, "    ")
  }
  nsteps = tt
  if (visu) plotTrajHK(traj, ne)
  return(list('steps' = tt, 'statesLast' = states, 'traj' = traj))
}

#example:
# ll = runDistModelHK(epsilon = 0.15, ne = 100, nstepsM = 1000)

#------------------------------------
# visualises trajectories of clusters as soon as clusters are formed.
plotTrajHK = function(traj, ne, maxt = 200, elines = FALSE, leg = FALSE){
  
  xs = NULL
  ys = NULL
  cols = NULL
  pchs = NULL
  cexs = NULL
  ng = length(traj[[1]])
  ieps = ng/ne
  if (traj[[1]][trunc(ng/2)] == 1) rd = 0
  else rd = 1
  ttraj = length(traj)
  div = max(1,trunc((ttraj-1) / maxt))
  nn = trunc((ttraj -1) / div)
  # print(nn)
  for(p in (0:nn)){
    tt = 1 + div*p
    dist = traj[[tt]]
    yps = which(dist > 0)
    if (length(yps) > nn / 10){
      xs = c(xs, rep(tt, length(yps)))
      #ys = c(ys, (yps-0.5)/ng)
      ys = c(ys, ieps*(((yps-0.5)/ng)-0.5))
      #for(yp in yps) cols = c(cols, colDens(dist[yp]/(2*ne)))
      for(yp in yps) cols = c(cols, colPw(dist[yp]/(2*ne)))
      for(yp in yps) pchs = c(pchs, pchPw(dist[yp]/(2*ne)))
      for(yp in yps) cexs = c(cexs, cexPw(dist[yp]/(2*ne)))
    }
    else{
      cls = clustFromAttHK(dist)
      xs = c(xs, rep(tt, length(cls$x)))
      ys = c(ys, ieps*(cls$x-0.5))
      for(j in 1:length(cls$x)) cols = c(cols, colPw(cls$val[j]/(2*ne)))
      for(j in 1:length(cls$x)) pchs = c(pchs, pchPw(cls$val[j]/(2*ne)))
      for(j in 1:length(cls$x)) cexs = c(cexs, cexPw(cls$val[j]/(2*ne)))
    }
  }
  dist = traj[[ttraj]]
  cls = clustFromAttHK(dist)
  for (k in(1:round(maxt*0.2))){
    xs = c(xs, rep(tt+k*div, length(cls$x)))
    ys = c(ys, ieps*(cls$x-0.5))
    for(j in 1:length(cls$x)) cols = c(cols, colPw(cls$val[j]/(2*ne)))
    for(j in 1:length(cls$x)) pchs = c(pchs, pchPw(cls$val[j]/(2*ne)))
    for(j in 1:length(cls$x)) cexs = c(cexs, cexPw(cls$val[j]/(2*ne)))
  }
  #print(length(xs))
  par(mar=c(2.1, 2.1, 0.5, 0.5), cex = 1.5)
  plot(xs, ys, xlab = 'iterations (N(t) interactions)', ylab = 'opinion', main = '', col = cols, pch = pchs, cex = cexs)
  #y1 = 0.5 + ne / ng
  if (elines){
    y1 = 1
    lines(c(-1, 2*maxt*div), c(y1, y1), lwd = 2)
    #y1 = 0.5 - ne / ng
    y1 = -1
    lines(c(-1, 2*maxt*div), c(y1, y1), lwd = 2)
  }
  texts = c('1', '0.75', '0.5', '0.25', '<0.01')
  if (leg)  legend(tt, 1, texts, pch = c(19, 19, 19, 19, 3), col = c(colPw(1), colPw(0.75), colPw(0.5), colPw(0.25), colPw(0.01)), cex = 0.8)
  
}

#------------------
colPw = function(pw){
  
  pw = min(1, pw)
  if (pw > 0.5){
    vv = (pw - 0.5) * 2
    return(rgb(vv, 1 - vv, 0))
  } 
  else{
    vv = (0.5 - pw)*2
    return(rgb(0, (1- vv), vv))
  }
}

#------------------
pchPw = function(pw){
  
  if (pw > 0.01) return (19)
  return (3)
}

#------------------
cexPw = function(pw){
  
  if (pw > 0.01) return (0.4)
  return (0.4)
}

#---------------------------------
# Explores the final positions of clusters (with max number of steps grid size^2).
# ne: discretization of epsilon
# epsInv1 : min value of 1/epsilon
# epsInv2 : max value of 1/epsilon
# nEps number of values of epsilon tested
# rd = 0 no noise, if rd > 0, factor of noise (small receomended e.g. 0.01)
#------------------------------
explEpsilonDHK = function(ne, epsInv1, epsInv2, nEps, odd = 0, rd = 0){
  
  evs = epsInv1 + ((epsInv2 - epsInv1)*(0:(nEps-1))/(nEps-1))
  nte = NULL
  npe = NULL
  cvl = NULL
  states = list()
  #nGrid = nGrid0
  res = list()
  for (i in (1:nEps)){
    epsilon = 0.5/evs[i]
    nGrid = round(ne / mod$epsilon)
    if ((nGrid %% 2) == 0) nGrid = nGrid +1
    steps = nGrid^2
    cat("\r                                             ", evs[i], "      ", i, "     ")
    sd =runDistModelHK(epsilon, ne, steps, odd = odd, visu = FALSE, save = FALSE, rd)

    states[[length(states)+1]] = sd$statesLast
    cls = clustFromAttHK(sd$statesLast)
    res[[i]] = cls
    names(res)[i]= epsilon
    nte = c(nte, length(cls$x))
    npe = c(npe, length(cls$pw[cls$pw > 0.6]))
    cvl = c(cvl, sd$steps)
  }
  # plot(evs, nte, ylim = c(0, max(nte)), xlab = '1/(2 epsilon)', ylab = 'n clusters')
  # points(evs, npe, col = 'red')
  plotAttPositionsCont2(res)
  return(list('nt' = nte, 'np' = npe, 'all' = res, 'cvl' = cvl, 'states' = states))
}

#Example:
#ll = explEpsilonDHK(100,  1.5, 4.5, 150, rd = 0.01)

#-----------------
#plot positions of final clusters when 1/epsilon varies.
#-------------------
plotAttPositionsCont2 = function(exEp, xl, yl, cx = 0.3, leg = TRUE){
  
  ne = length(exEp)
  xs = NULL
  ys = NULL
  cols = NULL
  pchs = NULL
  #cexs = NULL
  epss = 1/ (2* (as.numeric(names(exEp))))
  if (missing(yl)) yl = epss[ne]
  if (missing(xl)) xl = c(epss[1],epss[ne])
  for (i in (1:ne)){
    cls = exEp[[i]]
    xx = epss[i]
    #print(xx)
    #pw = cls$pw
    # pw = cls$val * xx
    pw = cls$val / (2*ne)
    for (j in (1:length(cls$x))){
      xs = c(xs, xx)
      yy = (cls$x[j] - 0.5)* 2 * xx
      ys = c(ys, yy)
      cols = c(cols, colPw(cls$pw[j]))
      pchs = c(pchs, pchPw(cls$pw[j]))
      # cexs = c(cexs, cexPw(cls$pw[j]))
    }
  }
  par(mar=c(2.1, 2.1, 0.5, 0.5), cex = 1.5)
  ld = 1
  plot(xs, ys, xlim = xl, ylim = c(-yl,yl), col = cols, pch = pchs, lwd = ld, cex = cx)
  if(leg) legend(xl[1], 0.95*yl, c(">=1", "0.5","0.1", "<0.01"), pch = c(19, 19, 19, 3), lwd = c(1, 1, 1, 1), col = c("red", "green", "blue","blue"), lty = c(NA, NA, NA, NA), cex = 0.9)
  
}

