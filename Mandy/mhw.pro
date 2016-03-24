;=======================================================
   PRO MHW, Npix_side, scale, MHW, FFT_MHW, power_MHW
;=======================================================
;
;      4 Nov. 2001
;
;      Makes a Mexican Hat Wavelet (MHW).
;
;

;   print,' MHW  :'
;   print,'-----------'
;   print,'            Scale (~ fwhm) -> ',scale,' pixels '
;   print,'            Npix_side      -> ',Npix_side,' pixels'
;   print,' ' 

   MHW     = FLTARR(Npix_side,Npix_side)
   FFT_MHW = COMPLEXARR(Npix_side,Npix_side)

;==================================
;  * Make  MHW  in  real  space  *
;==================================
 
   Distance      = DIST(Npix_side)/scale
   MHW           = (2.0 - Distance^2)*EXP(-Distance^2/2.0) 
   Normalization = TOTAL(MHW)
   MHW           = MHW/Normalization
   
;===============
;  * FFT MHW  *
;===============

   FFT_MHW = FFT(MHW,-1)

;==========================
;  *  Re-Scale  FFT MHW  *
;==========================
;
;   FFT_MHW = FFT_MHW/(ABS(FFT_MHW(0,0)))
;

;======================
;  *  Power of MHW  * 
;======================

   ff        = ROUND(DIST(Npix_side))  ; fourier k-mode filter
   nf        = FIX(SQRT(2.0*((FLOAT(Npix_side)/2.0)^2)))
   power_MHW = FLTARR(nf)

   FOR i = 1,(nf-1) DO BEGIN 

       p = WHERE(ff eq i) 

       Mm = MOMENT(ABS(FFT_MHW(p))^2) 
       power_MHW(i) = Mm(0)  ;mean value

   ENDFOR

;   print,' MHW Power Spectrum '
;   plot, power,/XLOG, /YLOG, XRANGE=[0.5,400.0] ;, YRANGE= [1.0E-21,1.0E-11]
;   print,' MHW'



END

