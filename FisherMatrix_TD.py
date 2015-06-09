import numpy as np
import scipy.integrate
import scipy.misc

## fiducial cosmology

Om=0.3
OL=0.7
Ok=1-Om-OL
h=0.7
w0=-1.0
wa=0.0

par_f=np.array((h,OL,Ok,w0,wa))

## experiment design

N_lens=4000; #LSST
fTc=0.0064 #d_Tc=0.64%, uncertainty
#fTc=0.00037 # Suyu 2014

## mock lenses

zl_m=0.5
zl_s=0.15
zs_m=2.0
zs_s=0.75

##### setting ends here #####


## bootstrap
lens=np.zeros((N_lens,2))

i=0
while i<(N_lens):
	lens_g=zl_s*np.random.randn()+zl_m
	sour_g=zs_s*np.random.randn()+zs_m
	

	if lens_g<sour_g:
		lens[i,0]=lens_g
		lens[i,1]=sour_g
		i=i+1


##Time delay measurable
def Tc(zl,zs,par):
#par=5x1 array of comsological parameters (h,OL,Ok,w0,wa)
 	h,OL,Ok,w0,wa=par[0],par[1],par[2],par[3],par[4]
	Om=1-Ok-OL
	
	#inverse E(z)
	rEz=lambda x: (Om*(1.+x)**3+Ok*(1+x)**2+OL*(1+x)**(3*(1+w0+wa))*np.exp(-3.0*wa*x/(1+x)))**(-0.5)

	#dimentionless ratio E_ratio
	El,er=scipy.integrate.quad(rEz,0,zl)
	Es,er=scipy.integrate.quad(rEz,0,zs)
	Els,er=scipy.integrate.quad(rEz,zl,zs)
	
	#curvature
	if Ok<0.0:
		cur=(np.absolute(Ok))**0.5
		El=np.sin(cur*El)/cur
		Es=np.sin(cur*Es)/cur
		Els=np.sin(cur*Els)/cur
	elif Ok>0.0:
		cur=(np.absolute(Ok))**0.5
		El=np.sinh(cur*El)/cur
		Es=np.sinh(cur*Es)/cur
		Els=np.sinh(cur*Els)/cur

	E_ratio=El*Es/Els

	#Tc
	Tc=E_ratio/(100*h)
	return Tc

######

# Tc of fiducial cosmology
Tc_f=np.zeros(N_lens)

for i in range(N_lens):
	Tc_f[i]=Tc(lens[i,0],lens[i,1],par_f)
#	print Tc_f[i]

######
 
## 2D Fisher matrix

# chi-square
def chisq(par):
#par=5x1 array of comsological parameters (h,OL,Ok,w0,wa)
 	#h,OL,Ok,w0,wa=par[0],par[1],par[2],par[3],par[4]

	dv=par-par_f #difference to fiducial cosmology
	sig2=np.dot(dv,dv)

	chi=0.
    #fTc=0.0064 #d_Tc=0.64%
	for i in range(N_lens):
	
		Tc_i=Tc(lens[i,0],lens[i,1],par)
		#print Tc_i, Tc_f[i]
#		chi=chi+(Tc_i-Tc_f[i])**2/((Tc_i*fTc)**2+sig2)
		chi=chi+(Tc_i-Tc_f[i])**2/((Tc_i*fTc)**2)


		
	return chi


# change parameter values

def cpar(df,i):
# d=interval fraction, i=index of vary parameter, vary from fiducial
# return 2 arrays: x-dx & x+dx and dx
	#h,OL,Ok,w0,wa=par[0],par[1],par[2],par[3],par[4]

	par_m=np.array((h,OL,Ok,w0,wa)) #minus array
	par_p=np.array((h,OL,Ok,w0,wa)) #plus array

	#if par_m[i] !=0.0:
	#	dx=par_m[i]*df
	#else:
	dx=df
	
	par_m[i]=par_m[i]-dx
	par_p[i]=par_p[i]+dx

	return par_m,par_p

def cpar2(d,i,j):
# return 4 arrays: par(dx,dy),(-dx,dy),(dx,-dy),(-dx,-dy)
	par_1=np.array((h,OL,Ok,w0,wa))
	par_2=np.array((h,OL,Ok,w0,wa))
	par_3=np.array((h,OL,Ok,w0,wa))
	par_4=np.array((h,OL,Ok,w0,wa))

	par_1[i],par_1[j]=par_1[i]+d,par_1[j]+d
	par_2[i],par_2[j]=par_2[i]-d,par_2[j]+d
	par_3[i],par_3[j]=par_3[i]+d,par_3[j]-d
	par_4[i],par_4[j]=par_4[i]-d,par_4[j]-d

	return par_1,par_2,par_3,par_4

## chi-square differential

#d^2/dx^2

def f_ii(i):
# i: index of dx
	dx=0.001
	#v=cpar(0.011,i)
	v=cpar(dx,i)
	par_m,par_p=v[0],v[1]
	
	a,b= chisq(par_p),chisq(par_m)
#	print chisq(par_p),chisq(par_m)

	der=(b+a)/(1000*(2*dx)**2)

	return der

def f_ij(i,j):
# i: index of dx, j: index of dy
	dx=0.001
	v=cpar2(dx,i,j)
	par_1,par_2,par_3,par_4=v[0],v[1],v[2],v[3]

	a,b,c,d=chisq(par_1),chisq(par_2),chisq(par_3),chisq(par_4)
#	print a-b-c+d

	der=(a-b-c+d)/(1000*(4*dx)**2)
	
	return der



############

#print f_ii(2)/2.
#print f_ij(2,1)/2.

# Fisher Matrix
F=np.zeros((5,5))

for i in range(5):
	for j in range(5):
		if i == j: #diagonal
			F[i,j]=f_ii(i)/2.
        elif i<j:
			F[i,j]=f_ij(i,j)/2.
			F[j,i]=F[i,j]

print F


 
	
