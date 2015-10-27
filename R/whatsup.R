gettelescopes=function(name){
  telescopes = NULL
  data('telescopes',envir = environment())
  allownames=tolower(as.character(telescopes[,'Name']))
  if(tolower(name) %in% allownames==FALSE){stop(paste('Provided name name is not allowed, must be one of',paste(as.character(telescopes[,'Ref']),sep='',collapse=', '),' (case insensitive). See ?telescopes for details.'))}
  out=telescopes[allownames==tolower(name),]
  names(out)=colnames(telescopes)
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

whatsup=function(RA="12:30:16", Dec="-30:13:15", Date='get', Time=c(12,0,0), lon=115.8178, lat=-31.97333, UTCdiff=8, altitude=10, pressure=1000, temp=20, step=0.5){
  if(is.character(lon)){
    obs=gettelescopes(lon)
    lon=obs['lon']
    lat=obs['lat']
    altitude=obs['Height']
  }
  RAdeg=hms2deg(RA)
  Decdeg=dms2deg(Dec)
  options(longitude=lon, latitude=lat)
  localEoT=lt()[1]-lst(lambda=lon)[1]
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
    tempout=eq2hor(RAdeg, Decdeg, jd , lat=lat, lon=lon, altitude=altitude, pres=pressure, temp=temp+273.15)
    tempoutmoon=eq2hor(moon$ra, moon$dec, jd , lat=lat, lon=lon, altitude=altitude, pres=pressure, temp=temp+273.15)
    tempoutsun=eq2hor(sun$ra, sun$dec, jd , lat=lat, lon=lon, altitude=altitude, pres=pressure, temp=temp+273.15)
    sink(NULL)
    outLT=jd2date(jd+UTCdiff/24)
    outLST=lst(jd, lambda=lon)
    outLTPOSIX=ISOdatetime(outLT$year, outLT$mon, outLT$mday, as.numeric(deg2hms(outLT$hour*15)[1]), as.numeric(deg2hms(outLT$hour*15)[2]), as.numeric(deg2hms(outLT$hour*15)[3]))
    outtemp=data.frame(JD=startjd+i, LST=outLST[1], LT=outLT, LTPOSIX=outLTPOSIX, Alt=tempout$alt, Az=tempout$az, HA=tempout$ha, AirMass=airmass(tempout$alt), AltMoon=tempoutmoon$alt, AzMoon=tempoutmoon$az, HAMoon=tempoutmoon$ha, AirMassMoon=airmass(tempoutmoon$alt), AltSun=tempoutsun$alt, AzSun=tempoutsun$az, HASun=tempoutsun$ha, AirMassSun=airmass(tempoutsun$alt))
    r=as.POSIXct(round(range(outLTPOSIX,na.rm = TRUE), "hours"))
    out=rbind(out,outtemp)
  }
  N=dim(out)[1]
  r=as.POSIXct(round(range(out$LTPOSIX), "hours"))
  ataxis=seq(r[1], r[2], by = "hour")
  dayline=ISOdatetime(out$LT.year[N], out$LT.mon[N], out$LT.mday[N], 0, 0, 0)
  obsGST=gst(date2jd(out$LT.year[N], out$LT.mon[N], out$LT.mday[N]))[1]
  LSTdiff=(gst(out$JD[floor(N/2)],0)[1]+lon/15-UTCdiff)
  tempsun1=sun.rst(jday = floor(out$JD[1]), phi = lat)
  tempsun2=sun.rst(jday = floor(out$JD[N]), phi = lat)
  rise1=(as.numeric(deg2hms((tempsun1[[1]][1]*15-LSTdiff*15) %% 360)))
  rise2=(as.numeric(deg2hms((tempsun2[[1]][1]*15-LSTdiff*15) %% 360)))
  set1=(as.numeric(deg2hms((tempsun1[[3]][1]*15-LSTdiff*15) %% 360)))
  set2=(as.numeric(deg2hms((tempsun2[[3]][1]*15-LSTdiff*15) %% 360)))
  rise=c(ISOdatetime(out$LT.year[1], out$LT.mon[1], out$LT.mday[1], rise1[1], rise1[2], rise1[3]), ISOdatetime(out$LT.year[N], out$LT.mon[N], out$LT.mday[N], rise2[1], rise2[2], rise2[3]))
  set=c(ISOdatetime(out$LT.year[1], out$LT.mon[1], out$LT.mday[1], set1[1], set1[2], set1[3]), ISOdatetime(out$LT.year[N], out$LT.mon[N], out$LT.mday[N], set2[1], set2[2], set2[3]))
  out=list(obs=out, UTCdiff=UTCdiff, RA=RA, Dec=Dec, RAdeg=RAdeg, Decdeg=Decdeg, lon=lon, lat=lat, altitude=altitude, pressure=pressure, temp=temp, at=ataxis, dayline=dayline, rise=rise, set=set, moonsep=moonsep, moonphase=moonphase, sunsep=sunsep)
  class(out)='whatsup'
  return=out
}

plot.whatsup=function(x,...){
  if(class(x)!='whatsup'){stop('Object must be of type whatsup')}
  plotwhatsup(x,...)
}

plotwhatsup=function(obsdata, ytype='Alt'){
  layout(rbind(1,2))
  par(mar=c(0,0,0,0))
  par(oma=c(3.1,3.1,3.1,2.1))
  if(ytype=='Alt'){
    magplot(obsdata$obs$LTPOSIX, obsdata$obs$Alt, xaxt='n', type='l', ylim=c(-30,90), xlab='', ylab='Alt / deg', tcl=0.5, mgp=c(2,0.5,0))
    lines(obsdata$obs$LTPOSIX, obsdata$obs$AltMoon, col='grey')
    lines(obsdata$obs$LTPOSIX, obsdata$obs$AltSun, col='orange')
    axis.POSIXct(1, obsdata$at, obsdata$at, format = '%H', tcl=0.5, mgp=c(2,0.5,0))
    abline(v=obsdata$dayline, lty=1)
    abline(v=obsdata$rise, lty=2, col='orange')
    abline(v=obsdata$set, lty=2, col='orange')
    abline(v=obsdata$obs$LTPOSIX[which.max(obsdata$obs$Alt)], lty=3, col='blue')
    abline(h=c(-15,0,30,60,90), lty=c(3,1,3,3,1))
  }
  if(ytype=='AM'){
    magplot(obsdata$obs$LTPOSIX, obsdata$obs$AirMass, xaxt='n', type='l', ylim=c(3,1), xlab='', ylab='Air Mass', tcl=0.5, mgp=c(2,0.5,0))
    lines(obsdata$obs$LTPOSIX, obsdata$obs$AirMassMoon, col='grey')
    lines(obsdata$obs$LTPOSIX, obsdata$obs$AirMassSun, col='orange')
    axis.POSIXct(1, obsdata$at, obsdata$at, format = '%H', tcl=0.5, mgp=c(2,0.5,0))
    abline(v=obsdata$dayline, lty=1)
    abline(v=obsdata$rise, lty=2, col='orange')
    abline(v=obsdata$set, lty=2, col='orange')
    abline(v=obsdata$obs$LTPOSIX[which.max(obsdata$obs$Alt)], lty=3, col='blue')
    abline(h=c(2,1.15,1), lty=c(3,3,1))
  }

  title(main=paste('Y', obsdata$obs$LT.year[1], 'M', obsdata$obs$LT.mon[1], 'D', obsdata$obs$LT.mday[1], ', RA ', round(obsdata$RAdeg),', Dec ', round(obsdata$Decdeg), ', Lon ', round(obsdata$lon), ', Lat ', round(obsdata$lat),', Mphase ',round(obsdata$moonphase,2),', Msep ', round(obsdata$moonsep), ', Ssep ', round(obsdata$sunsep),sep=''), outer=TRUE)

  magplot(obsdata$obs$LTPOSIX, obsdata$obs$Az, xaxt='n', type='l', xlab='Local Time from Now', ylab='Az (0N90E180S270W) / deg',ylim=c(min(obsdata$obs$Az),max(obsdata$obs$Az)+diff(range(obsdata$obs$Az))*0.2))
  axis.POSIXct(1, obsdata$at, obsdata$at, format = '%H', tcl=0.5, mgp=c(2,0.5,0))
  abline(v=obsdata$dayline, lty=1)
  abline(v=obsdata$rise, lty=2, col='orange')
  abline(v=obsdata$set, lty=2, col='orange')
  abline(v=obsdata$obs$LTPOSIX[which.max(obsdata$obs$Alt)], lty=3, col='blue')
  abline(h=c(0,90,180,270), lty=2)
}
