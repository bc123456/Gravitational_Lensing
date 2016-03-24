PRO MHW_filter



                                ;Sent by Chema
                                ;2015

                                ;This is for filtering an image with Mexian Hat Wavelet
                                ;Need to be used with MHW.pro
  
  

N=512 ; I assume the map you are using is 512x512. Change if necessary

scale=3.0 ; pixels. You can play with this scale

dat = '/Users/Mandy/desktop/Fortran/Case10/Case7/OptSN_GalFit_a/20.dat' 

map = fltarr(512,512)
mapS = fltarr(512,512)
openr, lun, dat, /GET_LUN
readf, lun, map
free_lun, lun

MHW, N, scale, MHW, FFT_MHW, power_MHW

FFTm = FFT(map,-1)*FFT_MHW  ; map is your NxN 2D map with the cloud of clusters
mapS = FLOAT(FFT(FFTm,1)) ; Smoothed map in real space
mapS = -mapS
writefits,'/Users/Mandy/desktop/Fortran/Case10/Case8/OptSN_GalFit_1.3.3/20_MHW.fits',mapS

END
