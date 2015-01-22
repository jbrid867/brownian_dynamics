from math import *

def vp(x,y,z,xx,yy,zz,vx,vy,vz):
	vi=[vx,vy,vz]
	drx=xx-x
	dry=yy-y
	drz=zz-z
	dr=(drx**2+dry**2+drz**2)**0.5
	print dr

	drx=drx/dr
	dry=dry/dr
	drz=drz/dr

	rhat=[drx,dry,drz]
	print rhat
	vpar=[0,0,0]
	


	v=vx*drx+vy*dry+vz*drz
	for i in range(3):
		vpar[i]=vi[i]-v*rhat[i]

	print vpar[i]

	print v
	print vx
	print v*rhat[0]+vpar[0]

	v1f=[vpar[0],vpar[1],vpar[2]]
	v2f=[v*rhat[0],v*rhat[1],v*rhat[2]]
	print v1f, v2f

	Ei=0
	E1f=0
	E2f=0;

	for i in range(3):
		Ei+=0.5*vi[i]*vi[i]
		E1f+=0.5*v1f[i]*v1f[i] 
		E2f+= 0.5*v2f[i]*v2f[i]

	print Ei
	print E1f
	print E2f

	return v



vperp=vp(8.15937e-11,-2.3652e-10,2.20108e-10,6.66667e-11,-3.33333e-10,2e-10,136.736,-731.749,324.214)




