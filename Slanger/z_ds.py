from math import sqrt as sqrt
import numpy as np

def z_ds(z_src):

    z_lens =0.375        # redshift of the lens 
    z_fid= 5.0   #fiducial redshift of delensed image
    #z_src= 0.93  #<===== change source redshift here
    Omega = 0.3          # Omega_m
    Lambda  = 0.7          # Omega_Lambda 
    Omega_K = 0.0
    DH      = 3000.0 
    if z_src==0: z_src=z_src+1e-6
    def integra(z):
        z_int  = 0.0
        Dc     = 0.0
        N_z    = 10000
        incr_z = z*1./N_z
        for i_z in range(N_z):
            z_int = z_int + incr_z 
            E_z=sqrt(Omega*(1+z_int)**3 + Lambda)         
            Dc = Dc + DH*incr_z/E_z 
        return Dc;
    Dm1 = integra(z_lens) 
    Dm2_true = integra(z_src)
    Ds_true = Dm2_true/(1.0 + z_src)   #<===============
    Dm2_fid = integra(z_fid)
    Ds_fid = Dm2_fid/(1.0 + z_fid)   #<=========

    DH      = 3000.0              
    Dds_true     = (Dm2_true-Dm1)/(1.0 + z_src)	#<=============
    Dds_fid     = (Dm2_fid-Dm1)/(1.0 + z_fid)	#<=============

    weight_true=Dds_true/Ds_true
    weight_fid=Dds_fid/Ds_fid
    return weight_true/weight_fid;

print "Redhisft= ",
red=float(raw_input())
print "Ratio= ",z_ds(red)

'''
data=[0.81,0.72,1.42,1.27,1.2215,1.0035,2.533,3.4294,1.5856]
for i in range(9):
	print z_ds(data[i])
'''

