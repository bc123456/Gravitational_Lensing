PRO MakeAlpha_ACSresolution_Fast


                                ;Mandy
                                ;2016

                                ;I wrote this for calculating alpha
                                ;map with high resolution (ACS resol)
                                ;in a fast way (FFT method)
  

  N = 1600
  FOV_arcmin = 0.8
  pix2rad    = ((FOV_arcmin/60.0)/FLOAT(N))*(3.1415927/180.0)
  
  
  lightfits = '~/FiducialMass/814CleanBCG_v2.fits'

  light = readfits(lightfits)
  light = GAUSS_SMOOTH(light,3)

  Xo = 30
  WINDOW,0,XS=N,YS=N,XP=0*N + 1*Xo,YP=500, RETAIN=2
  WINDOW,1,XS=N,YS=N,XP=1*N + 2*Xo,YP=500, RETAIN=2
;  
;  
;  
;  
;  
                                ; --------------------------------------------------
                                ;Make kernel
                               
  Kernel_X = FLTARR(N,N)
  Kernel_y = FLTARR(N,N)
  
; Create Kernel in radians units.
; ===============================

FOR i=0,N-1 DO BEGIN
    FOR j=0,N-1 DO BEGIN
        dX = (FLOAT(i-N/2) + 0.01)*pix2rad
        dY = (FLOAT(j-N/2) + 0.01)*pix2rad
        d2 = dX^2 + dY^2 
        Kernel_X(i,j) = dX/d2 
        Kernel_Y(i,j) = dY/d2 
    ENDFOR
ENDFOR    
; Center to border
; ================
AUX = Kernel_X
Kernel_X(0:N/2-1,0:N/2-1) = AUX(N/2:N-1,N/2:N-1)
Kernel_X(0:N/2-1,N/2:N-1) = AUX(N/2:N-1,0:N/2-1)
Kernel_X(N/2:N-1,0:N/2-1) = AUX(0:N/2-1,N/2:N-1)
Kernel_X(N/2:N-1,N/2:N-1) = AUX(0:N/2-1,0:N/2-1)
AUX = Kernel_Y
Kernel_Y(0:N/2-1,0:N/2-1) = AUX(N/2:N-1,N/2:N-1)
Kernel_Y(0:N/2-1,N/2:N-1) = AUX(N/2:N-1,0:N/2-1)
Kernel_Y(N/2:N-1,0:N/2-1) = AUX(0:N/2-1,N/2:N-1)
Kernel_Y(N/2:N-1,N/2:N-1) = AUX(0:N/2-1,0:N/2-1)

; FFT Kernel
; ==========
FFT_KX = FFT(Kernel_X,-1)
FFT_KY = FFT(Kernel_Y,-1)
  





                                ; --------------------------------------------------
                                ; Make alpha

  FFTm = FFT(light,-1)

  FFT_alphaX = FFTm*FFT_KX
  FFT_alphaY = FFTm*FFT_KY

  alphaX = (FLOAT(FFT(FFT_alphaX,1)))/2.0
  alphaY = (FLOAT(FFT(FFT_alphaY,1)))/2.0

;  wset,0 & TVSCL, alphaX
;  wset,1 & TVSCL, alphaY

  writefits,'~/HiRes_Light/MACS1149BCG_alphaX_FastCompute_Smooth3_v2.fits',alphaX
  writefits,'~/HiRes_Light/MACS1149BCG_alphaY_FastCompute_Smooth3_v2.fits',alphaY


print,size(alphaX)







END
