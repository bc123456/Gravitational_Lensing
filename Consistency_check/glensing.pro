;FUNCTION distance, z
;Omega=0.3
;Lambda=0.7
;Omega_k = 1.0 -Omega - Lambda
;	RETURN, 1.0/(SQRT(Lambda+Omega*(1.0+z)^(3.0)+Omega_k*(1.0+z)^2))
;END

pro glensing
 
;CHANGE THE PARAMETERS & FILE PATHS BELOW
  Npix = 3584
  Npix0 = 512 
  z_lens = DOUBLE(0.348)
  z_fid = DOUBLE(2.0)
Omega=0.3
Lambda=0.7
Omega_k = 1.0 -Omega - Lambda
Dh = 3000

  OPENR, 1, '/home/jess/RXJ2248/WSLAP/global_delens_relens/img_position.dat'
  N_img=0
  readf, 1, N_img
  position=dblarr(3,N_img)
  readf, 1, position
  close, 1

  openr, 2, '/home/jess/RXJ2248/WSLAP/global_delens_relens/sys'
  N_sys=0
  readf, 2, N_sys
  redshift=dblarr(3,N_sys)
  readf,2, redshift
  close, 2

z_src=fltarr(N_sys)
FOR x=0,N_sys-1 DO z_src(x)=redshift(1,x)
sys=intarr(N_sys)
sys(*)=redshift(2,*)

  openr, 3, '/home/jess/RXJ2248/WSLAP/global_delens_relens/f_k.dat'
  N=0
  readf, 3, N
  f_k_get=dblarr(3,N)
  readf,3, f_k_get
  close, 3
;print, f_k_get

  f_k_list=fltarr(N)
FOR x=0,N-1 DO BEGIN 
  f_k_list(x)=f_k_get(2,x) 
ENDFOR
  f_k=dblarr(2,N_sys)
  f_k(0,*)=sys(*)
FOR x=0,N_sys-1 DO BEGIN 
  f_k(1,x)=f_k_get(1,WHERE(f_k_list(*) eq z_src(x))) 
ENDFOR
;print, f_k
  
  position(1,*)=position(1,*)-709
  position(2,*)=position(2,*)-709

  alphaX_512_pix = READFITS('/home/jess/RXJ2248/WSLAP/global_delens_relens/recomp_alpha_x_pix.fits')
  alphaY_512_pix = READFITS('/home/jess/RXJ2248/WSLAP/global_delens_relens/recomp_alpha_y_pix.fits')
  theta='/home/jess/RXJ2248/WSLAP/global_delens_relens/theta.fits'
  src= '/home/jess/RXJ2248/WSLAP/global_delens_relens/beta.fits'
  outfile = '/home/jess/RXJ2248/WSLAP/global_delens_relens/ReImg.fits'

;================================================================== 
 
; = * = * = * = * = * = * = * = * = * = * = * = * = * = * = * 
  PRINT,' Interpolating alpha ... '
  X_ACS = FLOAT(Npix0)*FINDGEN(Npix)/FLOAT(Npix-1)
  Y_ACS = FLOAT(Npix0)*FINDGEN(Npix)/FLOAT(Npix-1)
  alphaX = INTERPOLATE(alphaX_512_pix,X_ACS,Y_ACS,/GRID)
  alphaY = INTERPOLATE(alphaY_512_pix,X_ACS,Y_ACS,/GRID)
  PRINT,' Interpolation done ' 

; = * = * = * = * = * = * = * = * = * = * = * = * = * = * = * 
;  PRINT,' Finding f_k ... '
  
;  z_src=dblarr(N_sys)
;  d_z2=dblarr(N_sys)
;  f_k=dblarr(2,N_sys)
;  Da_src=dblarr(N_sys)
;  Dm_src=dblarr(N_sys)
;  Da_src_lens=dblarr(N_sys)
  
;  z_src(*)=redshift(1,*)
;  Dm_fid=Dh*QROMB('distance', 0.0, z_fid)
;  Da_fid=Dm_fid/(1.0+z_fid)
;  Dm_src(*)=Dh*QROMB('distance', 0.0, z_src(*))
;  Da_src(*)=Dm_src(*)/(1.0+z_src(*))
;  Dm_lens=Dh*QROMB('distance', 0.0, z_lens)
;  Da_lens=Dm_lens/(1.0+z_lens)

;  SQRT_lens = SQRT(1.0 + Omega_k*(Dm_lens^2/Dh^2))
;  SQRT_fid = SQRT(1.0 + Omega_k*(Dm_fid^2/Dh^2))

;  Da_fid_lens=(Dm_fid*SQRT_lens-Dm_lens*SQRT_fid)/(1+z_fid)
;  Da_src_lens(*)=(Dm_src(*)*SQRT_lens-Dm_lens*SQRT(1.0+Omega_k*(Dm_src(*)^2/Dh^2)))/(1.0+z_src)

;  d_z2(*)=Da_src_lens(*)/Da_src(*)
;  d_z1=Da_fid_lens/Da_fid

;  d_z2(*)=QROMB('distance', z_lens, z_src(*))/QROMB('distance', 0.0, z_src(*))
;  d_z1=QROMB('distance', z_lens, z_fid)/QROMB('distance', 0.0, z_fid)

;  f_k(0,*)= redshift(2,*)
;  f_k(1,*)= d_z2(*)/d_z1

;  PRINT,' f_k'
;  PRINT, f_k


; = * = * = * = * = * = * = * = * = * = * = * = * = * = * = * 
  
  theta_map = FLTARR(Npix,Npix)
  Beta_map = FLTARR(Npix,Npix)
  ReImg= FLTARR(Npix,Npix)
  scaled_alphaX=FLTARR(Npix,Npix)
  scaled_alphaY=FLTARR(Npix,Npix)

  PRINT,'Delensing & Relensing ... '

  FOR i=0,N_sys-1 DO BEGIN
	  scaled_alphaX(*,*) = alphaX(*,*)*f_k(1,i)
	  scaled_alphaY(*,*) = alphaY(*,*)*f_k(1,i)
    FOR j=0, N_img-1 DO BEGIN
	BetaX=0
	BetaY=0
        IF FLOOR(position(0,j)) eq sys(i) THEN BEGIN 
	  x=ROUND(position(1,j))
	  y=ROUND(position(2,j))
	  theta_map(x,y)=position(0,j)
	  BetaX = ROUND(x-scaled_alphaX(x,y))
	  BetaY = ROUND(y-scaled_alphaY(x,y))
	  IF theta_map(x,y) gt Beta_map(BetaX, BetaY) THEN Beta_map(BetaX, BetaY)=theta_map(x,y)
          print, position(0,j), f_k(1,i)
	ENDIF
    ENDFOR
  ENDFOR

print, 'Finish Delensing, Start Relensing ...'
BetaX=0
BetaY=0
count=0
  FOR m=0,N_sys-1 DO BEGIN
	  scaled_alphaX(*,*) = alphaX(*,*)*f_k(1,m)
	  scaled_alphaY(*,*) = alphaY(*,*)*f_k(1,m)
;	print, 'sys', sys(m),  z_src(m), f_k(1,m)
    FOR k=0,Npix-1 DO BEGIN	 
      FOR l=0,Npix-1 DO BEGIN
	  BetaX = ROUND(k-scaled_alphaX(k,l))
	  BetaY = ROUND(l-scaled_alphaY(k,l))
        IF FLOOR(Beta_map(BetaX, BetaY)) eq sys(m) THEN BEGIN 
	  ReImg(k,l)=ReImg(k,l)+Beta_map(BetaX, BetaY)
	  count=count+1
;	  print, count, Beta_map(BetaX, BetaY), BetaX, BetaY
	ENDIF
      ENDFOR
;	  print,k+1
    ENDFOR
;	  print, 'count', count
  ENDFOR

  PRINT, ' Delens-Relens done, start printing maps ... '
  WRITEFITS, src, Beta_map
  WRITEFITS, theta, theta_map
  WRITEFITS,outfile,ReImg
  PRINT,' Relensed Images Printed in File '

END
  
  
   
  
  
