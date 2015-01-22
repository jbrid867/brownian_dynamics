from math import *

def vp(x,y,z,xx,yy,zz,vx,vy,vz):
     drx=xx-x
     dry=yy-y
     drz=zz-z
     dr=(drx**2+dry**2+drz**2)**0.5
     drx=drx/dr
     dry=dry/dr
     drz=drz/dr
     print dr
     v=vx*drx+vy*dry+vz*drz
     return v


vperp=vp(-2.29703e-10,1.95291e-10,-2.87041e-11,-2e-10,2e-10,6.66667e-11,-458.526,-257.297,507.839)



