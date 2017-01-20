whatsupeuler=function (ai, bi, select, fk4 = F, radian = F) 
{
    twopi = 2 * pi
    fourpi = 4 * pi
    rad_to_deg = 180/pi
    if (fk4) {
        equinox = "(b1950)"
        psi = c(0.57595865315, 4.9261918136, 0, 0, 0.11129056012, 
            4.7005372834)
        stheta = c(0.88781538514, -0.88781538514, 0.39788119938, 
            -0.39788119938, 0.86766174755, -0.86766174755)
        ctheta = c(0.46019978478, 0.46019978478, 0.9174369467, 
            0.9174369467, 0.49715499774, 0.49715499774)
        phi = c(4.9261918136, 0.57595865315, 0, 0, 4.7005372834, 
            0.11129056012)
    }
    else {
        equinox = "(j2000)"
        psi = c(0.574770433, 4.9368292465, 0, 0, 0.11142137093, 
            4.71279419371)
        stheta = c(0.88998808748, -0.88998808748, 0.39777715593, 
            -0.39777715593, 0.86766622025, -0.86766622025)
        ctheta = c(0.45598377618, 0.45598377618, 0.91748206207, 
            0.91748206207, 0.49714719172, 0.49714719172)
        phi = c(4.9368292465, 0.574770433, 0, 0, 4.71279419371, 
            0.11142137093)
    }
    i = select
    if (radian) {
        ao = ai - phi[i]
        bo = bi
    }
    else {
        ao = ai/rad_to_deg - phi[i]
        bo = bi/rad_to_deg
    }
    sb = sin(bo)
    cb = cos(bo)
    cbsa = cb * sin(ao)
    bo = -stheta[i] * cbsa + ctheta[i] * sb
    tmp = bo
    tmp[tmp > 1] = 1
    bo = asin(tmp)
    if (!radian) 
        bo = bo * rad_to_deg
    ao = atan2(ctheta[i] * cbsa + stheta[i] * sb, cb * cos(ao))
    ao = ((ao + psi[i] + fourpi)%%twopi)
    if (!radian) 
        ao = ao * rad_to_deg
    return(list(ao = ao, bo = bo))
}

getplanet=function(name='mars', JD=2440000){
  name=tolower(name)
  planID=which(c('mercury','venus','earth','mars','jupiter','saturn','uranus','neptune','pluto')==name)
  tmp = helio(JD, planID, radian = TRUE)
  rad = tmp$hrad
  lon = tmp$hlong
  lat = tmp$hlat
  tmp = helio(JD, 3, radian = TRUE)
  rade = as.numeric(tmp$hrad)
  lone = as.numeric(tmp$hlong)
  late = as.numeric(tmp$hlat)
  x = rad * cos(lat) * cos(lon) - rade * cos(late) * cos(lone)
  y = rad * cos(lat) * sin(lon) - rade * cos(late) * sin(lone)
  z = rad * sin(lat) - rade * sin(late)
  lambda = atan2(y, x) * 180/pi
  beta = atan2(z, sqrt(x * x + y * y)) * 180/pi
  tmp = whatsupeuler(lambda, beta, 4)
  RAdeg=as.numeric(tmp$ao)
  Decdeg=as.numeric(tmp$bo)
  return(data.frame(Name=name, RA=deg2hms(RAdeg,type='cat'), Dec=deg2dms(Decdeg,type='cat'), RAdeg=RAdeg, Decdeg=Decdeg, Type='P'))
}

gettelescope=function(name){
  telescopes = NULL
  data('telescopes',envir = environment())
  telescopes$Name=as.character(telescopes$Name)
  allownames=tolower(telescopes$Name)
  if(tolower(name) %in% allownames==FALSE){stop(paste('Provided telescope name is not allowed, must be one of',paste(telescopes$Name,sep='',collapse=', '),' (case insensitive). See ?telescopes for details.'))}
  out=telescopes[allownames==tolower(name),]
  names(out)=colnames(telescopes)
  out=as.vector(out)
  return(out)
}

gettarget=function(name, JD=2440000){
  targets = NULL
  data('targets',envir = environment())
  allownames=c(tolower(targets$Name),c('mercury','venus','earth','mars','jupiter','saturn','uranus','neptune','pluto'))
  if(tolower(name) %in% allownames==FALSE){
    name=gsub(' ','_',name)
    out=as.vector(nameresolve(name))
    if(is.na(out['RA'])){
      out=as.vector(data.frame(Name=name, RA="0:0:0", Dec="0:0:0", RAdeg=0, Decdeg=0, Type='R'))
      cat('Provided target name is not allowed!\n\n')
    }
  }else{
    if(tolower(name) %in% c('mercury','venus','earth','mars','jupiter','saturn','uranus','neptune','pluto')){
      out=getplanet(name,JD)
      out=as.vector(out)
    }else{
      out=targets[allownames==tolower(name),]
      names(out)=colnames(targets)
      out=as.vector(out)
    }
  }
  return(out)
}

nameresolve=function(name="M31"){
  RAdeg=rep(NA,length(name))
  Decdeg=rep(NA,length(name))
  RAhms=rep(NA,length(name))
  Decdms=rep(NA,length(name))
  for(i in 1:length(name)){
    temp=GET("http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/-oxp/~SNV",query=name[i])
    temp2=xmlTreeParse(as.character(temp))
    check=length(grep('Nothing',unlist(temp2)))
    check2=length(grep('refused',unlist(temp2)))
    check3=length(grep('Multiple',unlist(temp2)))
    if(check3>0){cat('Ambiguous name, using first returned from Seaseme!\n\n')}
    if(check>0 | check2>0){
      RAdeg[i]=NA
      Decdeg[i]=NA
    }else{
      RAdeg[i]=as.numeric(xmlValue(temp2$doc$children$Sesame[['Target']][['Resolver']][['jradeg']]))
      Decdeg[i]=as.numeric(xmlValue(temp2$doc$children$Sesame[['Target']][['Resolver']][['jdedeg']]))
      RAhms[i]=deg2hms(RAdeg[i], type='cat')
      Decdms[i]=deg2dms(Decdeg[i], type='cat')
    }
  }
  return(data.frame(Name=name, RA=RAhms, Dec=Decdms, RAdeg=RAdeg, Decdeg=Decdeg, Type='R'))
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
    mon[mon>12]=mon[mon>12]-12
    year = floor(C-4715)
    year[mon>=3]=year[mon>=3]-1
    return(list(year=year, mon=mon, mday=mday, hour=hour))
}

dayup=function(RA="12:30:16", Dec="-30:13:15", Target='user', Date='get', Time=c(12,0,0), Lon=115.8178, Lat=-31.97333, Loc='user', UTCdiff='get', Altitude=10, Pressure=1000, Temp=20, step=0.1){
  #no change
  planets=c('mercury','venus','earth','mars','jupiter','saturn','uranus','neptune','pluto')
  if(Loc != 'user'){
    obs=gettelescope(Loc)
    Lon=as.numeric(obs[1,'Lon'])
    Lat=as.numeric(obs[1,'Lat'])
    Altitude=as.numeric(obs[1,'Height'])
  }
  if(UTCdiff=='get'){UTCdiff=(lt()[1]-gmt()[1])}
  if(UTCdiff=='guess'){UTCdiff=round(Lon/15)}
  options(longitude=Lon, latitude=Lat)
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
  if(Target != 'user'){
    temp=gettarget(Target,startjd)
    RA=as.character(temp[1,'RA'])
    Dec=as.character(temp[1,'Dec'])
    RAdeg=temp[1,'RAdeg']
    Decdeg=temp[1,'Decdeg']
  }else{
    RAdeg=hms2deg(RA)
    Decdeg=dms2deg(Dec)
  }
  moon=moonpos(startjd)
  sun=sunpos(startjd)
  moonsep=acos(sum(sph2car(RAdeg,Decdeg)*sph2car(moon$ra, moon$dec)))*180/pi
  sunsep=acos(sum(sph2car(RAdeg,Decdeg)*sph2car(sun$ra, sun$dec)))*180/pi
  moonphase=mphase(startjd)
  jdadd=seq(0,1,by=step/24)
  Njd=length(jdadd)
  jd=startjd+jdadd
  sink("aux")
  Tempout=suppressWarnings(eq2hor(rep(RAdeg,Njd), rep(Decdeg,Njd), jd, lat=Lat, lon=Lon, altitude=Altitude, pres=Pressure, temp=Temp+273.15))
  Tempoutmoon=suppressWarnings(eq2hor(rep(moon$ra,Njd), rep(moon$dec,Njd), jd , lat=Lat, lon=Lon, altitude=Altitude, pres=Pressure, temp=Temp+273.15))
  Tempoutsun=suppressWarnings(eq2hor(rep(sun$ra,Njd), rep(sun$dec,Njd), jd , lat=Lat, lon=Lon, altitude=Altitude, pres=Pressure, temp=Temp+273.15))
  sink(NULL)
  outLT=jd2date(jd+UTCdiff/24)
  newdatetime=jd2date(jd)
  outLST=rep(NA,Njd)
  for(i in 1:Njd){outLST[i]=lst(date2jd(newdatetime$year[i],newdatetime$mon[i],newdatetime$mday[i],0), newdatetime$hour[i])}
  outLTPOSIX=ISOdatetime(outLT$year, outLT$mon, outLT$mday, as.numeric(deg2hms(outLT$hour*15)[,1]), as.numeric(deg2hms(outLT$hour*15)[,2]), as.numeric(deg2hms(outLT$hour*15)[,3]))
  out=data.frame(JD=jd, LST=as.numeric(outLST), LT=outLT, LTPOSIX=outLTPOSIX, Alt=Tempout$alt, Az=Tempout$az, HA=Tempout$ha, AirMass=airmass(Tempout$alt), AltMoon=Tempoutmoon$alt, AzMoon=Tempoutmoon$az, HAMoon=Tempoutmoon$ha, AirMassMoon=airmass(Tempoutmoon$alt), AltSun=Tempoutsun$alt, AzSun=Tempoutsun$az, HASun=Tempoutsun$ha, AirMassSun=airmass(Tempoutsun$alt))
  
  N=dim(out)[1]
  r=as.POSIXct(round(range(out$LTPOSIX, na.rm=TRUE), "hours"))
  ataxis=seq(r[1], r[2], by = "hour")
  dayline=ISOdatetime(out$LT.year[N], out$LT.mon[N], out$LT.mday[N], 0, 0, 0)
  obsGST=gst(date2jd(out$LT.year[N], out$LT.mon[N], out$LT.mday[N]))[1]
  LSTdiff=(gst(out$JD[floor(N/2)],0)[1]+Lon/15-UTCdiff) %% 24
  Tempsun1=sun.rst(jday = floor(out$JD[1]), phi = Lat)
  Tempsun2=sun.rst(jday = floor(out$JD[N]), phi = Lat)
  if(any(Tempoutsun$alt<0)){
  rise1=(as.numeric(deg2hms((Tempsun1[[1]][1]*15-LSTdiff*15) %% 360)))
  rise2=(as.numeric(deg2hms((Tempsun2[[1]][1]*15-LSTdiff*15) %% 360)))
  set1=(as.numeric(deg2hms((Tempsun1[[3]][1]*15-LSTdiff*15) %% 360)))
  set2=(as.numeric(deg2hms((Tempsun2[[3]][1]*15-LSTdiff*15) %% 360)))
  rise=c(ISOdatetime(out$LT.year[1], out$LT.mon[1], out$LT.mday[1], rise1[1], rise1[2], rise1[3]), ISOdatetime(out$LT.year[N], out$LT.mon[N], out$LT.mday[N], rise2[1], rise2[2], rise2[3]))
  set=c(ISOdatetime(out$LT.year[1], out$LT.mon[1], out$LT.mday[1], set1[1], set1[2], set1[3]), ISOdatetime(out$LT.year[N], out$LT.mon[N], out$LT.mday[N], set2[1], set2[2], set2[3]))
  }else{
    rise=c(NA,NA)
    set=c(NA,NA)
  }
  out=list(obs=out, UTCdiff=UTCdiff, RA=RA, Dec=Dec, RAdeg=RAdeg, Decdeg=Decdeg, Lon=Lon, Lat=Lat, Altitude=Altitude, Pressure=Pressure, Temp=Temp, at=ataxis, dayline=dayline, rise=rise, set=set, moonsep=moonsep, moonphase=moonphase, sunsep=sunsep)
  class(out)='dayup'
  return=out
}

dayupmulti=function(Targets=cbind(Names=c('G02', 'G09', 'G12', 'G15', 'G23'), RA=c("02:00:00", "09:00:00", "12:00:00", "14:30:00", "23:00:00"), Dec=c("-05:00:00", "01:00:00", "-01:00:00", "01:00:00", "-32:30:00")),...){
  out=list()
  for(i in 1:length(Targets[,1])){
    out=c(out,list(dayup(RA=Targets[i,tolower(colnames(Targets))=='ra'], Dec=Targets[i,tolower(colnames(Targets))=='dec'],...)))
  }
  out=c(out,list(names=Targets[,tolower(colnames(Targets))=='names']))
  class(out)="dayupmulti"
  return=out
}

plot.dayup=function(x,...){
  if(class(x)!='dayup'){stop('Object must be of type dayup')}
  plotdayup(x,...)
}

plot.dayupmulti=function(x,...){
  if(class(x)!='dayupmulti'){stop('Object must be of type dayupmulti')}
  plotdayupmulti(x,...)
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

nowup=function(Ncut=20, Azlim=c(0,360), Altlim=c(60,90), Date='get', Time='get', Lon=115.8178, Lat=-31.97333, Loc='user', UTCdiff='get', Altitude=10, Pressure=1000, Temp=20, select=c('G','S','C')){
  if(Loc != 'user'){
    obs=gettelescope(Loc)
    Lon=as.numeric(obs[1,'Lon'])
    Lat=as.numeric(obs[1,'Lat'])
    Altitude=as.numeric(obs[1,'Height'])
  }
  if(UTCdiff=='get'){UTCdiff=(lt()[1]-gmt()[1])}
  if(UTCdiff=='guess'){UTCdiff=round(Lon/15)}
  data('targets',envir = environment())
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
  exactjd=date2jd(year, mon, mday, hour=hour+min/60+sec/3600)-UTCdiff/24
  exactlst=lst(date2jd(year,mon,mday,0), hour+min/60+sec/3600-UTCdiff, lambda=Lon)[1]
  cutcat=targets[targets$Type %in% select,]
  HA=exactlst*15-cutcat$RAdeg
  HA=HA%%360
  altaz=hadec2altaz(HA, cutcat$Decdeg, lat=Lat, ws=FALSE)
  if(Azlim[1]<Azlim[2]){Azsel=altaz$az>Azlim[1] & altaz$az<Azlim[2]}
  if(Azlim[1]>Azlim[2]){Azsel=altaz$az>Azlim[1] | altaz$az<Azlim[2]}
  Altsel=altaz$alt>Altlim[1] & altaz$alt<Altlim[2]
  cutcat=cbind(cutcat[Azsel & Altsel,], HA=HA[Azsel & Altsel]/15, Alt=altaz$alt[Azsel & Altsel], Az=altaz$az[Azsel & Altsel])
  cutcat=cutcat[order(cutcat$Alt, decreasing=TRUE),]
  N=dim(cutcat)[1]
  if(N>Ncut){cutcat=cutcat[1:Ncut,]}
  return(cutcat)
}

nowwhere=function(RA="12:30:16", Dec="-30:13:15", Target='user', Date='get', Time='get', Lon=115.8178, Lat=-31.97333, Loc='user', UTCdiff='get', Altitude=10, Pressure=1000, Temp=20, Name='Auto'){
  if(Loc != 'user'){
    obs=gettelescope(Loc)
    Lon=as.numeric(obs[1,'Lon'])
    Lat=as.numeric(obs[1,'Lat'])
    Altitude=as.numeric(obs[1,'Height'])
  }
  if(UTCdiff=='get'){UTCdiff=(lt()[1]-gmt()[1])}
  if(UTCdiff=='guess'){UTCdiff=round(Lon/15)}

  options(longitude=Lon, latitude=Lat)
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
  
  exactjd=date2jd(year, mon, mday, hour=hour+min/60+sec/3600)-UTCdiff/24
  exactlst=lst(date2jd(year,mon,mday,0), hour+min/60+sec/3600-UTCdiff, lambda=Lon)[1]
  
  if(Target != 'user'){
    temp=gettarget(Target,exactjd)
    RA=as.character(temp[1,'RA'])
    Dec=as.character(temp[1,'Dec'])
  }
  RAdeg=hms2deg(RA)
  Decdeg=dms2deg(Dec)
  
  moon=moonpos(exactlst)
  sun=sunpos(exactlst)
  moonsep={}
  sunsep={}
  for(i in 1:length(RAdeg)){
    moonsep=c(moonsep,acos(sum(sph2car(RAdeg[i],Decdeg[i])*sph2car(moon$ra, moon$dec)))*180/pi)
    sunsep=c(sunsep,acos(sum(sph2car(RAdeg[i],Decdeg[i])*sph2car(sun$ra, sun$dec)))*180/pi)
  }
  moonphase=mphase(exactlst)
  sink("aux")
  targetobs=suppressWarnings(eq2hor(RAdeg, Decdeg, rep(exactjd,length(RAdeg)), lat=Lat, lon=Lon, altitude=Altitude, pres=Pressure, temp=Temp+273.15))
  moonobs=suppressWarnings(eq2hor(moon$ra, moon$dec, exactjd, lat=Lat, lon=Lon, altitude=Altitude, pres=Pressure, temp=Temp+273.15))
  sunobs=suppressWarnings(eq2hor(sun$ra, sun$dec, exactjd, lat=Lat, lon=Lon, altitude=Altitude, pres=Pressure, temp=Temp+273.15))
  sink(NULL)
  if(Name[1]=='Auto'){Name=paste('Target',1:length(RAdeg),sep='-')}
  obs=data.frame(
    Name=Name,
    RA=deg2hms(RAdeg, type='cat'),
    Dec=deg2dms(Decdeg, type='cat'),
    Alt=targetobs$alt,
    AirMass=airmass(targetobs$alt),
    Az=targetobs$az,
    HA=targetobs$ha/15,
    MoonSep=moonsep,
    Sunsep=sunsep
  )
  moonsun=data.frame(
    Name=c('Sun', 'Moon'),
    RA=deg2hms(c(sun$ra, moon$ra), type='cat'),
    Dec=deg2dms(c(sun$dec, moon$dec), type='cat'),
    Alt=c(sunobs$alt, moonobs$alt),
    AirMass=airmass(c(sunobs$alt, moonobs$alt)),
    Az=c(sunobs$az, moonobs$az),
    HA=c(sunobs$ha, moonobs$ha)/15
  )
  return(list(obs=obs, moonsun=moonsun, JD=exactjd, RJD=exactjd-2400000, LST=exactlst, LT.year=year, LT.mon=mon, LT.mday=mday, LT.hour=hour+min/60+sec/3600, UTCdiff=UTCdiff, Lon=Lon, Lat=Lat, Altitude=Altitude, Pressure=Pressure, Temp=Temp))
}

yearup=function(RA="12:30:16", Dec="-30:13:15", Target='user', Date='get', Lon=115.8178, Lat=-31.97333, Loc='user', UTCdiff='get', Altitude=10, Pressure=1000, Temp=20, select=c('G','S')){
  if(Loc != 'user'){
    obs=gettelescope(Loc)
    Lon=as.numeric(obs[1,'Lon'])
    Lat=as.numeric(obs[1,'Lat'])
    Altitude=as.numeric(obs[1,'Height'])
  }
  if(UTCdiff=='get'){UTCdiff=(lt()[1]-gmt()[1])}
  if(UTCdiff=='guess'){UTCdiff=round(Lon/15)}
  options(longitude=Lon, latitude=Lat)
  data('targets',envir = environment())
  if(Date[1]=='get'){
    datetime=as.POSIXlt(Sys.time())
    year=datetime$year+1900
    mon=datetime$mon+1
    mday=1
  }else{
    year=Date[1]
    mon=Date[2]
    mday=1
  }
  startjd=date2jd(year, mon, mday, hour=12)
  if(Target != 'user'){
    temp=gettarget(Target,startjd)
    RA=as.character(temp[1,'RA'])
    Dec=as.character(temp[1,'Dec'])
  }
  RAdeg=hms2deg(RA)
  Decdeg=dms2deg(Dec)
  jd=floor(startjd)+seq(0,365,by=5)
  uptime30={}
  uptime60={}
  moonphase={}
  moonsep={}
  for(i in jd){
  tempsun=sun.rst(jday = i, phi = Lat)
  gstlist=gst(i,hour=12,length=24*6,by=1/6) %% 24
  HA=gstlist*15-RAdeg
  HA=HA%%360
  altaz=hadec2altaz(HA, Decdeg, lat=Lat, ws=FALSE)
  if(tempsun$set>tempsun$rise){
    tempsun$set=tempsun$set-24
    gstlist[gstlist>tempsun$rise]=gstlist[gstlist>tempsun$rise]-24
  }
  afterset=gstlist>tempsun$set
  beforerise=gstlist<tempsun$rise
  uptime30=c(uptime30, length(which(afterset & beforerise & altaz$alt>30)))
  uptime60=c(uptime60, length(which(afterset & beforerise & altaz$alt>60)))
  moon=moonpos(i)
  moonphase=c(moonphase, mphase(i))
  moonsep=c(moonsep, acos(sum(sph2car(RAdeg,Decdeg)*sph2car(moon$ra, moon$dec)))*180/pi)
  }
  outdate=jd2date(jd)
  outLTPOSIX=ISOdatetime(outdate$year, outdate$mon, outdate$mday, 12, 0, 0)
  r=as.POSIXct(round(range(outLTPOSIX, na.rm=TRUE), "days"))
  ataxis=seq(r[1], r[2], by = "month")
  yearseq=expand.grid(c(year,year+1),1:12)
  yearlines=ISOdatetime(yearseq[,1], yearseq[,2], 1, 12, 0, 0)
  out=list(obs=data.frame(jd=jd, LTPOSIX=outLTPOSIX, up30=uptime30, up60=uptime60, moonsep=moonsep, moonphase=moonphase), at=ataxis, startyear=year, yearlines=yearlines)
  class(out)='yearup'
  return(out)
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

plot.yearup=function(x,...){
  if(class(x)!='yearup'){stop('Object must be of type yearup')}
  plotyearup(x,...)
}
