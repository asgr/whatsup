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
  if(!(all(Tempoutsun$alt<0) | all(Tempoutsun$alt>0))){
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
  return(out)
}

dayupmulti=function(Targets=cbind(Names=c('G02', 'G09', 'G12', 'G15', 'G23'), RA=c("02:00:00", "09:00:00", "12:00:00", "14:30:00", "23:00:00"), Dec=c("-05:00:00", "01:00:00", "-01:00:00", "01:00:00", "-32:30:00")),...){
  out=list()
  for(i in 1:length(Targets[,1])){
    out=c(out,list(dayup(RA=Targets[i,tolower(colnames(Targets))=='ra'], Dec=Targets[i,tolower(colnames(Targets))=='dec'],...)))
  }
  out=c(out,list(names=Targets[,tolower(colnames(Targets))=='names']))
  class(out)="dayupmulti"
  return(out)
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
  
  moon=moonpos(exactjd)
  sun=sunpos(exactjd)
  moonsep={}
  sunsep={}
  for(i in 1:length(RAdeg)){
    moonsep=c(moonsep,acos(sum(sph2car(RAdeg[i],Decdeg[i])*sph2car(moon$ra, moon$dec)))*180/pi)
    sunsep=c(sunsep,acos(sum(sph2car(RAdeg[i],Decdeg[i])*sph2car(sun$ra, sun$dec)))*180/pi)
  }
  moonphase=mphase(exactjd)
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
