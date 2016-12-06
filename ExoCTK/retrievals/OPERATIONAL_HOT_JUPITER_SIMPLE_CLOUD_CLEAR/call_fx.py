from fm import *
from matplotlib.pyplot import *
wlgrid=np.arange(0.3035,20,0.035)

# limits of parameters
# [T: 0 inf, Pref: 1 -4,  gases: -6  5, scatt amp: -5 5, scat slope(power): 0 5
# [T   Pref Pc  log(alphaH2O) log(alphaCH4) log(alphaCO) log(alphaCO2) log(alphaNH3) log(alphaNaK)	log(scatt. amp.)	scat. slope(power)
x=[1250, -1,  1,         2.,         -1.,        -1.,             -1,        -1,		-2.  , 		        0.,		   4.       ]


#calling fx
#for i in range(5):
y_mod,y_hires,wnocrop=fx(x,wlgrid) 

plot(1E4/wnocrop, y_hires)
plot(wlgrid, y_mod)
axis([.3, 2.5, .99*np.min(y_hires), 1.01*np.max(y_hires)])
show()





