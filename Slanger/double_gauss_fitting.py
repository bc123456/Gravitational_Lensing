#from sklearn import mixture
import matplotlib.pyplot as plt
import matplotlib.mlab
import numpy as np
from pylab import *
from scipy.optimize import leastsq

#data = np.genfromtxt('/home/lkit/Programming/Python/gaussian_fit.dat', skiprows = 1)
#x = data[:, 0]
#y = data[:, 1]

def z_to_v(Max_z):

 PI = 3.1415927
 Omega  = 0.3
 Lambda = 0.7
 DH     = 3000.0/0.7  #units of h^-1   DH = 3000/h, h=0.7
 N_z    = 10000
 incr_z = Max_z/float(N_z)
 z_int  = 0.0
 Dc     = 0.0
 for i_z in range (1,N_z):
         z_int = z_int + incr_z
         E_z=math.sqrt((((1.0+z_int)**3)*Omega) + Lambda)
         Dc = Dc + DH*incr_z/E_z

 v=Dc*70

 return v;
##############################################################

data=[0]*179
f=open('/home/lkit/Data/2015Dec_ref/nasa_z.dat')
i=0
for line in f:
 data[i]=float(line)
 i+=1
f.close()

diff_z= [0]*179
for i in range (0,178): 
 diff_z[i]=z_to_v(Max_z=data[i])-z_to_v(Max_z=0.373189359903)

sel=diff_z

n,bins,patches=plt.hist(sel, bins=40, range=(-5000,5000), histtype='stepfilled',align='mid') #define range
#plt.show()

n_bins=[0]*(len(bins)-1)
for i in range(0,len(bins)-1):
 n_bins[i]=(bins[i]+bins[i+1])/2.

#print n,n_bins,patches
#print "n= ",len(n)
#print "n_bins= ",len(n_bins)
#print "patches=",len(patches)
                  
mean = sum(n)/len(n)      
su=0
for i in range(0,len(n)):
 su=su+ n[i]*((n_bins[i]-mean)**2)/len(n)		#sigma=sqrt(sum[(x_i-mean)^2]/N)
sigma = np.sqrt(su)				 #must use np. :P

def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2));

print "sigma= ",sigma

x=n_bins
y=n
def gaussian(x,c1, mu1, sigma1):
    res =   c1 * np.exp( - (x - mu1)**2.0 / (2.0 * sigma1**2.0) ) 
    return res


def double_gaussian( x, params ):
    c1, mu1, sigma1, c2, mu2, sigma2 = params
    res =   c1 * np.exp( - (x - mu1)**2.0 / (2.0 * sigma1**2.0) ) \
          + c2 * np.exp( - (x - mu2)**2.0 / (2.0 * sigma2**2.0) )
    return res

def double_gaussian_fit( params ):
    fit = double_gaussian( x, params )
    return (fit - y_proc)

print 'x,y= ',x,y

# Remove background.
y_proc = np.copy(y)
#y_proc[y_proc < 5] = 0.0

# Least squares fit. Starting values found by inspection.
fit = leastsq( double_gaussian_fit, [2.0,500.0,1500.0,2.0,0.0,1500.0] )	
#plot( x, y, c='b' )
a,b,c,d,e,f=fit[0]
print a,b,c,d,e,f
#plot( x, double_gaussian( x, fit[0] ), c='r' )
plot( x, gaussian( x, a,b,c ), c='r' )
plot( x, gaussian( x, d,e,f ), c='g' )
plt.show()


#amplitude,mean,sigma (80 bins)
# [4.0,-1000.0,1500.0,7.0,750.0,1500.0] 
#1.21567206803, 520.316785619, 2231.18908677 
#2.00215007678, 102.698232931, 993.850855972


