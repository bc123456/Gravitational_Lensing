;FUNCTION distance, z
;	RETURN, 1.0/(SQRT( 0.7 + 0.3 * (1.0+z)^(3.0)))
;END

pro mag

;CHANGE THE PARAMETERS & FILE PATHS BELOW
  Npix = 3584
  Npix0 = 512 
  z_lens = 0.348
  z_fid = 2.0
  
  alphaX_512_pix = READFITS('/home/jess/RXJ2248/WSLAP/mag/recomp_alpha_x_pix.fits')
  alphaY_512_pix = READFITS('/home/jess/RXJ2248/WSLAP/mag/recomp_alpha_y_pix.fits')
  outfile = '/home/jess/RXJ2248/WSLAP/mag/predicted_f.txt'
  Magnification_out='/home/jess/RXJ2248/WSLAP/mag/magnification.fits'

; = * = * = * = * = * = * = * = * = * = * = * = * = * = * = * 
  PRINT,' Interpolating alpha ... '
  
  alphaX0 = fltarr(Npix,Npix)
  alphaY0 = fltarr(Npix,Npix)
  X_ACS = FLOAT(Npix0)*FINDGEN(Npix)/FLOAT(Npix-1)
  Y_ACS = FLOAT(Npix0)*FINDGEN(Npix)/FLOAT(Npix-1)
  alphaX0 = INTERPOLATE(alphaX_512_pix,X_ACS,Y_ACS,/GRID)
  alphaY0 = INTERPOLATE(alphaY_512_pix,X_ACS,Y_ACS,/GRID)
  PRINT,' Interpolation done ' 

; = * = * = * = * = * = * = * = * = * = * = * = * = * = * = * 

;store the positions of images and relensed images in position.dat
  openr, 1, '/home/jess/RXJ2248/WSLAP/mag/position.dat'
  N_img=0
  readf, 1, N_img
  data=fltarr(12,N_img)
  readf,1, data
  close, 1

  openr, 2, '/home/jess/RXJ2248/WSLAP/mag/sys'
  N_sys=0
  readf, 2, N_sys
  redshift=dblarr(3,N_sys)
  readf,2, redshift
  close, 2

z_src=fltarr(N_sys)
z_src(*)=redshift(1,*)
sys=intarr(N_sys)
sys(*)=redshift(2,*)

print, z_src
print, sys

  openr, 3, '/home/jess/RXJ2248/WSLAP/mag/f_k.dat'
  N=0
  readf, 3, N
  f_k_get=dblarr(3,N)
  readf,3, f_k_get
  close, 3

  f_k_list=fltarr(N)
  f_k_list(*)=f_k_get(2,*) 

  f_k=dblarr(2,N_sys)
  f_k(0,*)=sys(*)
FOR x=0,N_sys-1 DO BEGIN 
  f_k(1,x)=f_k_get(1,WHERE(f_k_list(*) eq z_src(x))) 
ENDFOR

print, f_k


;define variables
predicted_f=fltarr(16,N_img)

flux_mag_ratio=0
l=0
FOR  y=0,N_img-1 DO BEGIN
	IF data(3,y) ne -99 and data(10,y) ne -99 and FLOOR(data(0,y)) ne 1 and FLOOR(data(0,y)) ne 2 and FLOOR(data(0,y)) ne 3 and FLOOR(data(0,y)) ne 14 then BEGIN
	  flux_mag_ratio=flux_mag_ratio+data(3,y)/(10^(-data(10,y)/2.5))
	  l=l+1
	ENDIF
ENDFOR
flux_mag_ratio=flux_mag_ratio/l

; = * = * = * = * = * = * = * = * = * = * = * = * = * = * = * 
FOR a=0,N_sys-1 DO BEGIN

;FOR y=0,N_img-1 DO BEGIN
;  z_src = data(5,y)

;  PRINT,' Scaling alpha ... '
  
;  d_z2=QROMB('distance', z_lens, z_src)/QROMB('distance', 0.0, z_src)
;  d_z1=QROMB('distance', z_lens, z_fid)/QROMB('distance', 0.0, z_fid)
;  f_k= d_z2/d_z1
  alphaX = fltarr(Npix,Npix)
  alphaY = fltarr(Npix,Npix)
  alphaX(*,*) = alphaX0(*,*)*f_k(1,a)
  alphaY(*,*) = alphaY0(*,*)*f_k(1,a)
   PRINT,' f_k ', f_k(1,a)
;   PRINT,' Sacling done '

; = * = * = * = * = * = * = * = * = * = * = * = * = * = * = * 
;  PRINT,' Calculating differentials ... '

  a_xx=fltarr(Npix,Npix)
  FOR i=0,Npix-1 DO $
  IF i eq 0 THEN a_xx(i,*)=(alphaX(i+1,*)-alphaX(i,*)) ELSE $
  IF i eq Npix-1 THEN a_xx(i,*)=(alphaX(i,*)-alphaX(i-1,*)) ELSE $
  a_xx(i,*)=(alphaX(i+1,*)-alphaX(i-1,*))/2
;  PRINT,' a_xx done '

  a_yy=fltarr(Npix,Npix)
  FOR j=0,Npix-1 DO $
  IF j eq 0 THEN a_yy(*,j)=(alphaY(*,j+1)-alphaY(*,j)) else $
  IF j eq Npix-1 THEN a_yy(*,j)=(alphaY(*,j)-alphaY(*,j-1)) else $
  a_yy(*,j)=(alphaY(*,j+1)-alphaY(*,j-1))/2
;  PRINT,' a_yy done '

  a_yx=fltarr(Npix,Npix)
  FOR k=0,Npix-1 DO $
  IF k eq 0 THEN a_yx(k,*)=(alphaY(k+1,*)-alphaY(k,*)) ELSE $
  IF k eq Npix-1 THEN a_yx(k,*)=(alphaY(k,*)-alphaY(k-1,*)) ELSE $
  a_yx(k,*)=(alphaY(k+1,*)-alphaY(k-1,*))/2
;  PRINT,' a_yx done '

  a_xy=fltarr(Npix,Npix)
  FOR l=0,Npix-1 DO $
  IF l eq 0 THEN a_xy(*,l)=(alphaX(*,l+1)-alphaX(*,l)) else $
  IF l eq Npix-1 THEN a_xy(*,l)=(alphaX(*,l)-alphaX(*,l-1)) else $
  a_xy(*,l)=(alphaX(*,l+1)-alphaX(*,l-1))/2
;  PRINT,' a_xy done '

  inverse_M=fltarr(Npix,Npix)
  inverse_M=(1-a_xx)*(1-a_yy)-(a_yx)*(a_xy)
  mag=fltarr(Npix,Npix)
  mag=1/inverse_M

WRITEFITS, Magnification_out, mag
  PRINT,' mag done '
print, where( inverse_M(*,*) eq 0.0)
print, where( mag(*,*) eq 0.0)

; = * = * = * = * = * = * = * = * = * = * = * = * = * = * = * 
;  PRINT,' Calculating predicted flux ... '

FOR y=0,N_img-1 DO BEGIN
IF FLOOR(data(0,y)) eq FLOOR(f_k(0,a)) THEN BEGIN
  predicted_f(0,y)=data(0,y)
  predicted_f(1,y)=data(5,y)
  predicted_f(2,y)=ABS(mag(ROUND(data(1,y)),ROUND(data(2,y))))
  predicted_f(3,y)=ABS(mag(ROUND(data(6,y)),ROUND(data(7,y))))
  predicted_f(4,y)=data(10,y)
  predicted_f(5,y)=data(11,y)
M_ratio= predicted_f(3,y)/predicted_f(2,y)
  predicted_f(6,y)=-2.5*ALOG10((data(8,y)+data(9,y))/flux_mag_ratio)-predicted_f(5,y)
  predicted_f(7,y)=-2.5*ALOG10((data(8,y)-data(9,y))/flux_mag_ratio)-predicted_f(5,y)
  predicted_f(8,y)= data(3,y)*M_ratio
  predicted_f(9,y)=-2.5*ALOG10(predicted_f(8,y)/flux_mag_ratio)
  predicted_f(10,y)=-2.5*ALOG10((data(3,y)+data(4,y))*M_ratio/flux_mag_ratio)-predicted_f(9,y)
  predicted_f(11,y)=-2.5*ALOG10((data(3,y)-data(4,y))*M_ratio/flux_mag_ratio)-predicted_f(9,y)
ENDIF
ENDFOR
predicted_f(12,*)=data(3,*)
predicted_f(13,*)=data(4,*)
predicted_f(14,*)=data(8,*)
predicted_f(15,*)=data(9,*)
  PRINT,' system '+ STRING(a+1)+' done ... '
; = * = * = * = * = * = * = * = * = * = * = * = * = * = * = * 
ENDFOR

PRINT,' Calculation done, start printing the predicted flux ... '

  fname=outfile
  OPENW, 1, fname 
  PRINTF, 1, "#delense	relense	M_0	M_relensed	Mag_0	Mag_obs	mag_obs_err_low	mag_obs_err_up	Flux_pre	mag_pre	mag_pre_err_low	mag_pre_err_up	flux_0	flux_err_0	flux_obs	flux_obs_err"
  PRINTF, 1, predicted_f, FORMAT='(16F)'
  CLOSE, 1

 PRINT,' predicted_f Printed in File '

END
  
