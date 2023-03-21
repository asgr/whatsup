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