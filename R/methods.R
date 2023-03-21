plot.dayup=function(x,...){
  if(class(x)!='dayup'){stop('Object must be of type dayup')}
  plotdayup(x,...)
}

plot.dayupmulti=function(x,...){
  if(class(x)!='dayupmulti'){stop('Object must be of type dayupmulti')}
  plotdayupmulti(x,...)
}

plot.yearup=function(x,...){
  if(class(x)!='yearup'){stop('Object must be of type yearup')}
  plotyearup(x,...)
}

plotdayup=function(obsdata, ytype='Alt',moonphase=TRUE,Name=''){
  layout(rbind(1,2))
  par(mar=c(0,0,0,0))
  par(oma=c(3.1,3.1,3.1,2.1))
  if(moonphase){mooncol=hsv(h=0, s=0, v=1-obsdata$moonphase)}else{mooncol='darkgrey'}
  if(ytype=='Alt'){
    magplot(obsdata$obs$LTPOSIX, obsdata$obs$Alt, xaxt='n', type='l', ylim=c(-10,90), xlab='', ylab='Alt / deg', tcl=0.5, mgp=c(2,0.5,0), col='blue', prettybase=30)
    lines(obsdata$obs$LTPOSIX, obsdata$obs$AltMoon, col=mooncol)
    lines(obsdata$obs$LTPOSIX, obsdata$obs$AltSun, col='orange')
    axis.POSIXct(1, obsdata$at, obsdata$at, format = '%H', tcl=0.5, mgp=c(2,0.5,0))
    abline(v=obsdata$dayline+c(-24,-12,0,12,24)*3600, lty=2, col='grey')
    abline(v=obsdata$dayline+c(-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,20,21,22,23)*3600, lty=3, col='grey')
    abline(v=obsdata$rise, lty=2, col='orange')
    abline(v=obsdata$set, lty=2, col='orange')
    abline(v=obsdata$obs$LTPOSIX[which.max(obsdata$obs$Alt)], lty=3, col='blue')
    abline(h=c(-10,0,10,20,30,40,50,60,70,80,90), lty=c(3,1,3,3,2,3,3,2,3,3,1), col='grey')
    text(obsdata$obs$LTPOSIX[which.max(obsdata$obs$Alt)], obsdata$obs$Alt[which.max(obsdata$obs$Alt)],Name)
  }
  if(ytype=='AM'){
    magplot(obsdata$obs$LTPOSIX, obsdata$obs$AirMass, xaxt='n', type='l', ylim=c(3,1), xlab='', ylab='Air Mass', tcl=0.5, mgp=c(2,0.5,0),col='blue')
    lines(obsdata$obs$LTPOSIX, obsdata$obs$AirMassMoon, col=mooncol)
    lines(obsdata$obs$LTPOSIX, obsdata$obs$AirMassSun, col='orange')
    axis.POSIXct(1, obsdata$at, obsdata$at, format = '%H', tcl=0.5, mgp=c(2,0.5,0))
    abline(v=obsdata$dayline+c(-24,-12,0,12,24)*3600, lty=2, col='grey')
    abline(v=obsdata$dayline+c(-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,20,21,22,23)*3600, lty=3, col='grey')
    abline(v=obsdata$rise, lty=2, col='orange')
    abline(v=obsdata$set, lty=2, col='orange')
    abline(v=obsdata$obs$LTPOSIX[which.max(obsdata$obs$Alt)], lty=3, col='blue')
    abline(h=airmass(c(-10,0,10,20,30,40,50,60,70,80,90)), lty=c(3,1,3,3,2,3,3,2,3,3,1), col='grey')
    text(obsdata$obs$LTPOSIX[which.max(obsdata$obs$Alt)], obsdata$obs$AirMass[which.max(obsdata$obs$Alt)],Name)
  }
  
  title(main=paste('Y', obsdata$obs$LT.year[1], 'M', obsdata$obs$LT.mon[1], 'D', obsdata$obs$LT.mday[1], ', RA ', round(obsdata$RAdeg),', Dec ', round(obsdata$Decdeg), ', Lon ', round(obsdata$Lon), ', Lat ', round(obsdata$Lat),', Mphase ',round(obsdata$moonphase,2),', Msep ', round(obsdata$moonsep), ', Ssep ', round(obsdata$sunsep),sep=''), outer=TRUE)
  
  select=which(abs(diff(obsdata$obs$Az))>200)
  moonselect=which(abs(diff(obsdata$obs$AzMoon))>200)
  sunselect=which(abs(diff(obsdata$obs$AzSun))>200)
  xtarget=obsdata$obs$LTPOSIX
  ytarget=obsdata$obs$Az
  xmoon=obsdata$obs$LTPOSIX
  ymoon=obsdata$obs$AzMoon
  xsun=obsdata$obs$LTPOSIX
  ysun=obsdata$obs$AzSun
  xtarget[select]=NA
  ytarget[select]=NA
  xmoon[moonselect]=NA
  ymoon[moonselect]=NA
  xsun[sunselect]=NA
  ysun[sunselect]=NA
  magplot(xtarget, ytarget, xaxt='n', yaxt='n', type='l', xlab='Local Time / Hours', ylab='Az / deg', ylim=c(0,380), col='blue')
  lines(xmoon, ymoon, col=mooncol)
  lines(xsun, ysun, col='orange')
  axis.POSIXct(1, obsdata$at, obsdata$at, format = '%H', tcl=0.5, mgp=c(2,0.5,0))
  magaxis(2,prettybase=90,labels=FALSE)
  axis(2, at=c(0,90,180,270,360), labels=c('0N','90E','180S','270W','360N'), tick=FALSE, tcl=0.5, mgp=c(2,0.5,0))
  abline(v=obsdata$dayline+c(-24,-12,0,12,24)*3600, lty=2, col='grey')
  abline(v=obsdata$dayline+c(-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,20,21,22,23)*3600, lty=3, col='grey')
  abline(v=obsdata$rise, lty=2, col='orange')
  abline(v=obsdata$set, lty=2, col='orange')
  abline(v=obsdata$obs$LTPOSIX[which.max(obsdata$obs$Alt)], lty=3, col='blue')
  abline(h=c(0,30,60,90,120,150,180,210,240,270,300,330,360), lty=c(1,3,3,2,3,3,2,3,3,2,3,3,1), col='grey')
  legend('bottomright', legend=c('Target',' Sun', 'Moon'), col=c('blue','orange',mooncol),lty=1)
  
}

plotdayupmulti=function(obsdata, ytype='Alt',moonphase=TRUE,Names){
  layout(rbind(1,2))
  par(mar=c(0,0,0,0))
  par(oma=c(3.1,3.1,3.1,2.1))
  if(moonphase){mooncol=hsv(h=0, s=0, v=1-obsdata[[1]]$moonphase)}else{mooncol='darkgrey'}
  if(ytype=='Alt'){
    magplot(obsdata[[1]]$obs$LTPOSIX, obsdata[[1]]$obs$Alt, xaxt='n', type='l', ylim=c(-10,90), xlab='', ylab='Alt / deg', tcl=0.5, mgp=c(2,0.5,0), col='blue', prettybase=30)
    for(i in 2:(length(obsdata)-1)){lines(obsdata[[i]]$obs$LTPOSIX, obsdata[[i]]$obs$Alt, col='blue')}
    lines(obsdata[[1]]$obs$LTPOSIX, obsdata[[1]]$obs$AltMoon, col=mooncol)
    lines(obsdata[[1]]$obs$LTPOSIX, obsdata[[1]]$obs$AltSun, col='orange')
    axis.POSIXct(1, obsdata[[1]]$at, obsdata[[1]]$at, format = '%H', tcl=0.5, mgp=c(2,0.5,0))
    abline(v=obsdata[[1]]$dayline+c(-24,-12,0,12,24)*3600, lty=2, col='grey')
    abline(v=obsdata[[1]]$dayline+c(-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,20,21,22,23)*3600, lty=3, col='grey')
    abline(v=obsdata[[1]]$rise, lty=2, col='orange')
    abline(v=obsdata[[1]]$set, lty=2, col='orange')
    abline(h=c(-10,0,10,20,30,40,50,60,70,80,90), lty=c(3,1,3,3,2,3,3,2,3,3,1), col='grey')
    for(i in 1:(length(obsdata)-1)){
      abline(v=obsdata[[i]]$obs$LTPOSIX[which.max(obsdata[[i]]$obs$Alt)], lty=3, col='blue')
      text(obsdata[[i]]$obs$LTPOSIX[which.max(obsdata[[i]]$obs$Alt)], obsdata[[i]]$obs$Alt[which.max(obsdata[[i]]$obs$Alt)],obsdata$names[i])
    }
  }
  if(ytype=='AM'){
    magplot(obsdata[[1]]$obs$LTPOSIX, obsdata[[1]]$obs$AirMass, xaxt='n', type='l', ylim=c(3,1), xlab='', ylab='Air Mass', tcl=0.5, mgp=c(2,0.5,0),col='blue')
    for(i in 2:(length(obsdata)-1)){lines(obsdata[[i]]$obs$LTPOSIX, obsdata[[i]]$obs$AirMass, col='blue')}
    lines(obsdata[[1]]$obs$LTPOSIX, obsdata[[1]]$obs$AirMassMoon, col=mooncol)
    lines(obsdata[[1]]$obs$LTPOSIX, obsdata[[1]]$obs$AirMassSun, col='orange')
    axis.POSIXct(1, obsdata[[1]]$at, obsdata[[1]]$at, format = '%H', tcl=0.5, mgp=c(2,0.5,0))
    abline(v=obsdata[[1]]$dayline+c(-24,-12,0,12,24)*3600, lty=2, col='grey')
    abline(v=obsdata[[1]]$dayline+c(-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,20,21,22,23)*3600, lty=3, col='grey')
    abline(v=obsdata[[1]]$rise, lty=2, col='orange')
    abline(v=obsdata[[1]]$set, lty=2, col='orange')
    abline(v=obsdata[[1]]$obs$LTPOSIX[which.max(obsdata[[1]]$obs$Alt)], lty=3, col='blue')
    abline(h=airmass(c(-10,0,10,20,30,40,50,60,70,80,90)), lty=c(3,1,3,3,2,3,3,2,3,3,1), col='grey')
    for(i in 1:(length(obsdata)-1)){
      abline(v=obsdata[[i]]$obs$LTPOSIX[which.max(obsdata[[i]]$obs$Alt)], lty=3, col='blue')
      text(obsdata[[i]]$obs$LTPOSIX[which.max(obsdata[[i]]$obs$Alt)], obsdata[[i]]$obs$AirMass[which.max(obsdata[[i]]$obs$Alt)],obsdata$names[i])
    }
  }
  
  title(main=paste('Y', obsdata[[1]]$obs$LT.year[1], 'M', obsdata[[1]]$obs$LT.mon[1], 'D', obsdata[[1]]$obs$LT.mday[1], ', Lon ', round(obsdata[[1]]$Lon), ', Lat ', round(obsdata[[1]]$Lat),', Mphase ',round(obsdata[[1]]$moonphase,2),sep=''), outer=TRUE)
  
  select=which(abs(diff(obsdata[[1]]$obs$Az))>200)
  moonselect=which(abs(diff(obsdata[[1]]$obs$AzMoon))>200)
  sunselect=which(abs(diff(obsdata[[1]]$obs$AzSun))>200)
  xtarget=obsdata[[1]]$obs$LTPOSIX
  ytarget=obsdata[[1]]$obs$Az
  xmoon=obsdata[[1]]$obs$LTPOSIX
  ymoon=obsdata[[1]]$obs$AzMoon
  xsun=obsdata[[1]]$obs$LTPOSIX
  ysun=obsdata[[1]]$obs$AzSun
  xtarget[select]=NA
  ytarget[select]=NA
  xmoon[moonselect]=NA
  ymoon[moonselect]=NA
  xsun[sunselect]=NA
  ysun[sunselect]=NA
  magplot(xtarget, ytarget, xaxt='n', yaxt='n', type='l', xlab='Local Time / Hours', ylab='Az / deg', ylim=c(0,380), col='blue')
  for(i in 2:(length(obsdata)-1)){
    select=which(abs(diff(obsdata[[i]]$obs$Az))>200)
    xtarget=obsdata[[i]]$obs$LTPOSIX
    ytarget=obsdata[[i]]$obs$Az
    xtarget[select]=NA
    ytarget[select]=NA
    lines(xtarget, ytarget, col='blue')
  }
  lines(xmoon, ymoon, col=mooncol)
  lines(xsun, ysun, col='orange')
  axis.POSIXct(1, obsdata[[1]]$at, obsdata[[1]]$at, format = '%H', tcl=0.5, mgp=c(2,0.5,0))
  magaxis(2,prettybase=90,labels=FALSE)
  axis(2, at=c(0,90,180,270,360), labels=c('0N','90E','180S','270W','360N'), tick=FALSE, tcl=0.5, mgp=c(2,0.5,0))
  abline(v=obsdata[[1]]$dayline+c(-24,-12,0,12,24)*3600, lty=2, col='grey')
  abline(v=obsdata[[1]]$dayline+c(-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,20,21,22,23)*3600, lty=3, col='grey')
  abline(v=obsdata[[1]]$rise, lty=2, col='orange')
  abline(v=obsdata[[1]]$set, lty=2, col='orange')
  abline(h=c(0,30,60,90,120,150,180,210,240,270,300,330,360), lty=c(1,3,3,2,3,3,2,3,3,2,3,3,1), col='grey')
  for(i in 1:(length(obsdata)-1)){
    abline(v=obsdata[[i]]$obs$LTPOSIX[which.max(obsdata[[i]]$obs$Alt)], lty=3, col='blue')
  }
  legend('bottomright', legend=c('Target',' Sun', 'Moon'), col=c('blue','orange',mooncol),lty=1)
  
}

plotyearup=function(obsdata){
  layout(rbind(1,2))
  par(mar=c(0,0,0,0))
  par(oma=c(3.1,3.1,3.1,3.1))
  magplot(obsdata$obs$LTPOSIX, obsdata$obs$up30/6, xaxt='n', type='l', ylim=c(0,24), xlab='', ylab='Time Up / Hrs', tcl=0.5, mgp=c(2,0.5,0), col='blue', prettybase=30)
  abline(v=obsdata$yearlines, lty=3, col='grey')
  abline(h=seq(0,24,by=3), lty=c(2,3), col='grey')
  lines(obsdata$obs$LTPOSIX, obsdata$obs$up60/6, col='red')
  legend('topright',legend=c('Alt > 30 deg','Alt > 60 deg'),lty=1,col=c('blue','red'))
  axis.POSIXct(1, obsdata$at, obsdata$at, tcl=0.5, mgp=c(2,0.5,0))
  lims=par()$usr
  par(usr=c(lims[1:2],0,200))
  moonphase=spline(obsdata$obs$LTPOSIX,obsdata$obs$moonphase)$y
  moonphase[moonphase<0]=0
  moonphase[moonphase>1]=1
  magplot(spline(obsdata$obs$LTPOSIX, obsdata$obs$moonsep), xaxt='n', type='l', ylim=c(0,180), xlab=paste(obsdata$startyear,'-',obsdata$startyear+1,sep=''), ylab='Moon Sep / Deg', tcl=0.5, mgp=c(2,0.5,0), col='grey', prettybase=30)
  abline(h=seq(0,180,by=30), lty=c(2,3), col='grey')
  abline(v=obsdata$yearlines, lty=3, col='grey')
  points(spline(obsdata$obs$LTPOSIX, obsdata$obs$moonsep), col=hsv(h=0, s=0, v=1-moonphase), pch=16)
  axis.POSIXct(1, obsdata$at, obsdata$at, tcl=0.5, mgp=c(2,0.5,0))
}
