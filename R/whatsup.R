gettelescope=function(name){
  telescopes = NULL
  data('telescopes',envir = environment())
  allownames=tolower(as.character(telescopes[,'Name']))
  if(tolower(name) %in% allownames==FALSE){stop(paste('Provided telescope name is not allowed, must be one of',paste(as.character(telescopes[,'Ref']),sep='',collapse=', '),' (case insensitive). See ?telescopes for details.'))}
  out=telescopes[allownames==tolower(name),]
  names(out)=colnames(telescopes)
  out=as.vector(out)
  return(out)
}

gettarget=function(name){
  targets = NULL
  data('targets',envir = environment())
  allownames=tolower(as.character(targets[,'Name']))
  if(tolower(name) %in% allownames==FALSE){stop('Provided target name is not allowed.')}
  out=targets[allownames==tolower(name),]
  names(out)=colnames(targets)
  out=as.vector(out)
  return(out)
}

airmass=function(Alt=30){
  out=1/(sin(Alt*pi/180)+0.50572*(6.07995+Alt)^(-1.6364))
  out[Alt<0]=NA
  return(out)
}

date2jd=function(year=1968, mon=5, mday=23, hour=12){
  year[mon %in% 1:2]=year[mon %in% 1:2]-1
  mon[mon %in% 1:2]=mon[mon %in% 1:2]+12
  A = floor(year/100)
  B = floor(A/4)
  C = floor(2-A+B)
  E = floor(365.25*(year+4716))
  F = floor(30.6001*(mon+1))
  return(C+mday+E+F+hour/24-1524.5)
}
jd2date=function(JD=2440000){
    Q = JD+0.5
    Z = floor(Q)
    W = floor((Z - 1867216.25)/36524.25)
    X = floor(W/4)
    A = floor(Z+1+W-X)
    B = floor(A+1524)
    C = floor((B-122.1)/365.25)
    D = floor(365.25*C)
    E = floor((B-D)/30.6001)
    F = floor(30.6001*E)
    fracmday=B-D-F+(Q-Z)
    hour= (fracmday-floor(fracmday))*24
    mday= floor(fracmday)
    mon = floor(E-1)
    mon = mon%%12
    year = floor(C-4715)
    year[mon>=3]=year[mon>=3]-1
    return(list(year=year, mon=mon, mday=mday, hour=hour))
}

whatsup=function(RA="12:30:16", Dec="-30:13:15", Target='user', Date='get', Time=c(12,0,0), Lon=115.8178, Lat=-31.97333, Loc='user', UTCdiff=8, Altitude=10, Pressure=1000, Temp=20, step=0.1){
  #no change
  if(Target != 'user'){
    RA=as.character(gettarget(Target)[1,2])
    Dec=as.character(gettarget(Target)[1,3])
  }
  if(Loc != 'user'){
    obs=gettelescope(Loc)
    Lon=as.numeric(obs[1,'Lon'])
    Lat=as.numeric(obs[1,'Lat'])
    Altitude=as.numeric(obs[1,'Height'])
  }
  if(UTCdiff=='guess'){UTCdiff=round(Lon/15)}
  RAdeg=hms2deg(RA)
  Decdeg=dms2deg(Dec)
  options(longitude=Lon, latitude=Lat)
  localEoT=lt()[1]-lst(lambda=Lon)[1]
  if(Date[1]=='get'){
    datetime=as.POSIXlt(Sys.time())
    year=datetime$year+1900
    mon=datetime$mon+1
    mday=datetime$mday
  }else{
    year=Date[1]
    mon=Date[2]
    mday=Date[3]
  }
  if(Time[1]=='get'){
    datetime=as.POSIXlt(Sys.time())
    hour=datetime$hour
    min=datetime$min
    sec=datetime$sec
  }else{
    hour=Time[1]
    min=Time[2]
    sec=Time[3]
  }
  startjd=date2jd(year, mon, mday, hour=hour+min/60+sec/3600)-UTCdiff/24
  moon=moonpos(startjd)
  sun=sunpos(startjd)
  moonsep=acos(sum(sph2car(RAdeg,Decdeg)*sph2car(moon$ra, moon$dec)))*180/pi
  sunsep=acos(sum(sph2car(RAdeg,Decdeg)*sph2car(sun$ra, sun$dec)))*180/pi
  moonphase=mphase(startjd)
  out=data.frame()
  for(i in seq(0,1,by=step/24)){
    sink("aux")
    jd=startjd+i
    Tempout=eq2hor(RAdeg, Decdeg, jd , lat=Lat, lon=Lon, altitude=Altitude, pres=Pressure, temp=Temp+273.15)
    Tempoutmoon=eq2hor(moon$ra, moon$dec, jd , lat=Lat, lon=Lon, altitude=Altitude, pres=Pressure, temp=Temp+273.15)
    Tempoutsun=eq2hor(sun$ra, sun$dec, jd , lat=Lat, lon=Lon, altitude=Altitude, pres=Pressure, temp=Temp+273.15)
    sink(NULL)
    outLT=jd2date(jd+UTCdiff/24)
    outLST=lst(jd, lambda=Lon)
    outLTPOSIX=ISOdatetime(outLT$year, outLT$mon, outLT$mday, as.numeric(deg2hms(outLT$hour*15)[1]), as.numeric(deg2hms(outLT$hour*15)[2]), as.numeric(deg2hms(outLT$hour*15)[3]))
    outTemp=data.frame(JD=startjd+i, LST=outLST[1], LT=outLT, LTPOSIX=outLTPOSIX, Alt=Tempout$alt, Az=Tempout$az, HA=Tempout$ha, AirMass=airmass(Tempout$alt), AltMoon=Tempoutmoon$alt, AzMoon=Tempoutmoon$az, HAMoon=Tempoutmoon$ha, AirMassMoon=airmass(Tempoutmoon$alt), AltSun=Tempoutsun$alt, AzSun=Tempoutsun$az, HASun=Tempoutsun$ha, AirMassSun=airmass(Tempoutsun$alt))
    r=as.POSIXct(round(range(outLTPOSIX,na.rm = TRUE), "hours"))
    out=rbind(out,outTemp)
  }
  N=dim(out)[1]
  r=as.POSIXct(round(range(out$LTPOSIX), "hours"))
  ataxis=seq(r[1], r[2], by = "hour")
  dayline=ISOdatetime(out$LT.year[N], out$LT.mon[N], out$LT.mday[N], 0, 0, 0)
  obsGST=gst(date2jd(out$LT.year[N], out$LT.mon[N], out$LT.mday[N]))[1]
  LSTdiff=(gst(out$JD[floor(N/2)],0)[1]+Lon/15-UTCdiff)
  Tempsun1=sun.rst(jday = floor(out$JD[1]), phi = Lat)
  Tempsun2=sun.rst(jday = floor(out$JD[N]), phi = Lat)
  rise1=(as.numeric(deg2hms((Tempsun1[[1]][1]*15-LSTdiff*15) %% 360)))
  rise2=(as.numeric(deg2hms((Tempsun2[[1]][1]*15-LSTdiff*15) %% 360)))
  set1=(as.numeric(deg2hms((Tempsun1[[3]][1]*15-LSTdiff*15) %% 360)))
  set2=(as.numeric(deg2hms((Tempsun2[[3]][1]*15-LSTdiff*15) %% 360)))
  rise=c(ISOdatetime(out$LT.year[1], out$LT.mon[1], out$LT.mday[1], rise1[1], rise1[2], rise1[3]), ISOdatetime(out$LT.year[N], out$LT.mon[N], out$LT.mday[N], rise2[1], rise2[2], rise2[3]))
  set=c(ISOdatetime(out$LT.year[1], out$LT.mon[1], out$LT.mday[1], set1[1], set1[2], set1[3]), ISOdatetime(out$LT.year[N], out$LT.mon[N], out$LT.mday[N], set2[1], set2[2], set2[3]))
  out=list(obs=out, UTCdiff=UTCdiff, RA=RA, Dec=Dec, RAdeg=RAdeg, Decdeg=Decdeg, Lon=Lon, Lat=Lat, Altitude=Altitude, Pressure=Pressure, Temp=Temp, at=ataxis, dayline=dayline, rise=rise, set=set, moonsep=moonsep, moonphase=moonphase, sunsep=sunsep)
  class(out)='whatsup'
  return=out
}

plot.whatsup=function(x,...){
  if(class(x)!='whatsup'){stop('Object must be of type whatsup')}
  plotwhatsup(x,...)
}

plotwhatsup=function(obsdata, ytype='Alt',moonphase=TRUE){
  layout(rbind(1,2))
  par(mar=c(0,0,0,0))
  par(oma=c(3.1,3.1,3.1,2.1))
  if(moonphase){mooncol=hsv(v=0, s=0, alpha=obsdata$moonphase)}else{mooncol='darkgrey'}
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
    abline(h=c(0,10,20,30,40,50,60,70,80,90), lty=c(1,3,3,2,3,3,2,3,3,1), col='grey')
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
  }

  title(main=paste('Y', obsdata$obs$LT.year[1], 'M', obsdata$obs$LT.mon[1], 'D', obsdata$obs$LT.mday[1], ', RA ', round(obsdata$RAdeg),', Dec ', round(obsdata$Decdeg), ', Lon ', round(obsdata$Lon), ', Lat ', round(obsdata$Lat),', Mphase ',round(obsdata$moonphase,2),', Msep ', round(obsdata$moonsep), ', Ssep ', round(obsdata$sunsep),sep=''), outer=TRUE)

  magplot(obsdata$obs$LTPOSIX, obsdata$obs$Az, xaxt='n', yaxt='n', type='l', xlab='Local Time from Now', ylab='Az / deg', ylim=c(0,380), col='blue')
  lines(obsdata$obs$LTPOSIX, obsdata$obs$AzMoon, col=mooncol)
  lines(obsdata$obs$LTPOSIX, obsdata$obs$AzSun, col='orange')
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
