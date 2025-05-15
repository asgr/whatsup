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
  data('telescopes', envir = environment())
  telescopes$Name = as.character(telescopes$Name)
  allownames = tolower(telescopes$Name)
  if(all(tolower(name) %in% allownames) == FALSE){stop(paste('Provided telescope name is not allowed, must be one of',paste(telescopes$Name,sep='',collapse=', '),' (case insensitive). See ?telescopes for details.'))}
  #out=telescopes[allownames==tolower(name),]
  #names(out)=colnames(telescopes)
  #out=as.vector(out)
  return(telescopes[allownames %in% tolower(name),])
}

gettarget=function(name, JD=2440000){
  targets = NULL
  data('targets',envir = environment())
  allownames=c(tolower(targets$Name),
               c('mercury','venus','earth','mars','jupiter','saturn','uranus','neptune','pluto'),
               'moon',
               'sun'
               )
  if(tolower(name) %in% allownames==FALSE){
    name=gsub(' ','_',name)
    out = nameresolve(name)
    if(is.na(out['RA'])){
      out = data.frame(Name=name, RA="0:0:0", Dec="0:0:0", RAdeg=0, Decdeg=0, Type='R')
      cat('Provided target name is not allowed!\n\n')
    }
  }else{
    if(tolower(name) %in% c('mercury','venus','earth','mars','jupiter','saturn','uranus','neptune','pluto')){
      out = getplanet(name,JD)
    }else if(tolower(name) == 'moon'){
      out = moonpos(JD)
      out = data.frame(Name=name, RA=deg2hms(out$ra, type='cat'), Dec=deg2dms(out$dec, type='cat'), RAdeg=out$ra, Decdeg=out$dec, Type='M')
    }else if(tolower(name) == 'sun'){
      out = sunpos(JD)
      out = data.frame(Name=name, RA=deg2hms(out$ra, type='cat'), Dec=deg2dms(out$dec, type='cat'), RAdeg=out$ra, Decdeg=out$dec, Type='S')
    }else{
      out=targets[allownames==tolower(name),]
      names(out)=colnames(targets)
      out=as.vector(out)
    }
  }
  return(out)
}

nameresolve=function(name="M31"){
  if(!requireNamespace("XML", quietly = TRUE)){
    stop('The XML package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  if(!requireNamespace("httr", quietly = TRUE)){
    stop('The httr package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  RAdeg=rep(NA,length(name))
  Decdeg=rep(NA,length(name))
  RAhms=rep(NA,length(name))
  Decdms=rep(NA,length(name))
  for(i in 1:length(name)){
    temp=httr::GET("http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/-oxp/~SNV",query=name[i])
    temp2=XML::xmlTreeParse(as.character(temp))
    check=length(grep('Nothing',unlist(temp2)))
    check2=length(grep('refused',unlist(temp2)))
    check3=length(grep('Multiple',unlist(temp2)))
    if(check3>0){cat('Ambiguous name, using first returned from Sesame!\n\n')}
    if(check>0 | check2>0){
      RAdeg[i]=NA
      Decdeg[i]=NA
    }else{
      RAdeg[i]=as.numeric(XML::xmlValue(temp2$doc$children$Sesame[['Target']][['Resolver']][['jradeg']]))
      Decdeg[i]=as.numeric(XML::xmlValue(temp2$doc$children$Sesame[['Target']][['Resolver']][['jdedeg']]))
      RAhms[i]=deg2hms(RAdeg[i], type='cat')
      Decdms[i]=deg2dms(Decdeg[i], type='cat')
    }
  }
  return(data.frame(Name=name, RA=RAhms, Dec=Decdms, RAdeg=RAdeg, Decdeg=Decdeg, Type='R'))
}