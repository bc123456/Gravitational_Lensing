;FUNCTION distance, z
;Omega=0.3
;Lambda=0.7
;Omega_k = 1.0 -Omega - Lambda
;	RETURN, 1.0/(SQRT(Lambda+Omega*(1.0+z)^(3.0)+Omega_k*(1.0+z)^2))
;END

pro lensing
 
;CHANGE THE PARAMETERS & FILE PATHS BELOW
  Npix = 3584
  Npix0 = 512 
  z_lens = 0.348
  z_fid = 2.0
;Omega=0.3
;Lambda=0.7
;Omega_k = 1.0 -Omega - Lambda
;Dh = 3000

  openr, 1, '/home/jess/RXJ2248/WSLAP/Delens_Relens/range.dat'
  N_range=0
  readf, 1, N_range
  heading=strarr(1)
  readf, 1, heading
  data=fltarr(6,N_range)
  readf,1, data
  close, 1

  openr, 2, '/home/jess/RXJ2248/WSLAP/Delens_Relens/sys'
  N_sys=0
  readf, 2, N_sys
  redshift=fltarr(3,N_sys)
  readf,2, redshift
  close, 2

  openr, 3, '/home/jess/RXJ2248/WSLAP/global_delens_relens/f_k.dat'
  N=0
  readf, 3, N
  f_k_get=dblarr(3,N)
  readf,3, f_k_get
  close, 3
  f_k_list=fltarr(N)
FOR x=0,N-1 DO f_k_list(x)=f_k_get(2,x)

  openr, 4, '/home/jess/RXJ2248/WSLAP/Delens_Relens/obj.dat'
  N_obj=0
  readf, 4, N_obj
  get = strarr(N_obj)
  readf, 4, get
  close, 4

  alphaX_512_pix = READFITS('/home/jess/RXJ2248/WSLAP/Delens_Relens/recomp_alpha_x_pix.fits')
  alphaY_512_pix = READFITS('/home/jess/RXJ2248/WSLAP/Delens_Relens/recomp_alpha_y_pix.fits')

FITS_OPEN, '/home/jess/RXJ2248/WSLAP/Delens_Relens/rgb_FOV.fits', fcb
FITS_READ, fcb, r, header1, exten_no=1
FITS_READ, fcb, g, header2, exten_no=2
FITS_READ, fcb, b, header3, exten_no=3
FITS_CLOSE, fcb

; = * = * = * = * = * = * = * = * = * = * = * = * = * = * = * 
  PRINT,' Interpolating alpha ... '
  X_ACS = FLOAT(Npix0)*FINDGEN(Npix)/FLOAT(Npix-1)
  Y_ACS = FLOAT(Npix0)*FINDGEN(Npix)/FLOAT(Npix-1)
  alphaX = INTERPOLATE(alphaX_512_pix,X_ACS,Y_ACS,/GRID)
  alphaY = INTERPOLATE(alphaY_512_pix,X_ACS,Y_ACS,/GRID)
  PRINT,' Interpolation done ' 

; = * = * = * = * = * = * = * = * = * = * = * = * = * = * = * 

FOR x=0,N_obj-1 DO BEGIN
obj =  STRMID(get(x),0,strlen(get(x)))
sys =  STRMID(obj, 0, strlen(obj)-1) 

  IF STRMID(obj, strlen(obj)-1, 1) eq 'a' THEN obj_no=0.1 ELSE $
  IF STRMID(obj, strlen(obj)-1, 1) eq 'b' THEN obj_no=0.2 ELSE $
  IF STRMID(obj, strlen(obj)-1, 1) eq 'c' THEN obj_no=0.3 ELSE $
  IF STRMID(obj, strlen(obj)-1, 1) eq 'd' THEN obj_no=0.4 ELSE $
  IF STRMID(obj, strlen(obj)-1, 1) eq 'e' THEN obj_no=0.5
  
  z_src=0.0
  z_src = redshift(1,WHERE(redshift(2,*) eq FIX(sys)))
  row=WHERE(data(0,*) eq FLOAT(sys)+obj_no)
;  noise_lvl = FLOAT(data(5,row(0)))
  X_min=ROUND(data(1,row)-709)
  X_max=ROUND(data(2,row)-709)
  Y_min=ROUND(data(3,row)-709)
  Y_max=ROUND(data(4,row)-709)
  X_min=FIX(X_min(0))
  X_max=FIX(X_max(0))
  Y_min=FIX(Y_min(0))
  Y_max=FIX(Y_max(0))

  theta='/home/jess/RXJ2248/WSLAP/Delens_Relens/'+sys+'/'+obj+'.fits'
  src= '/home/jess/RXJ2248/WSLAP/Delens_Relens/'+sys+'/'+'src_'+ obj+'.fits'
  outfile = '/home/jess/RXJ2248/WSLAP/Delens_Relens/'+sys+'/'+'ReImg_'+obj +'.fits'
;================================================================== 
 




  PRINT,' Scaling alpha ... '

;z_step  = 70000
;z_start = 0.005
;E_z_f_int = 0.0
;E_z_s_int = 0.0
;E_z_l_int = 0.0

;  FOR iz=0,z_step DO BEGIN
;    z_f=z_start+iz*(z_fid/z_step)
;    z_s=z_start+iz*(z_src/z_step)
;    z_l=z_start+iz*(z_lens/z_step)
;    E_z_f_int =E_z_f_int+(z_fid/z_step)/(SQRT(Omega*(1.0+z_f)^3+Omega_k*(1.0+z_f)^2+Lambda))
;    E_z_s_int =E_z_s_int+(z_src/z_step)/(SQRT(Omega*(1.0+z_s)^3+Omega_k*(1.0+z_s)^2+Lambda))
;    E_z_l_int =E_z_l_int+(z_lens/z_step)/(SQRT(Omega*(1.0+z_l)^3+Omega_k*(1.0+z_l)^2+Lambda))
;  ENDFOR

;  Dm_fid=Dh*E_z_f_int
;  Da_fid=Dm_fid/(1.0+z_fid)
;  Dm_src=Dh*E_z_s_int
;  Da_src=Dm_src/(1.0+z_src)
;  Dm_lens=Dh*E_z_l_int
;  Da_lens=Dm_lens/(1.0+z_lens)

;  SQRT_lens = SQRT(1.0 + Omega_k*(Dm_lens^2/Dh^2))
;  SQRT_fid = SQRT(1.0 + Omega_k*(Dm_fid^2/Dh^2))

;  Da_fid_lens=(Dm_fid*SQRT_lens-Dm_lens*SQRT_fid)/(1+z_fid)
;  Da_src_lens=(Dm_src*SQRT_lens-Dm_lens*SQRT(1.0+Omega_k*(Dm_src^2/Dh^2)))/(1.0+z_src)

;  d_z2=Da_src_lens/Da_src
;  d_z1=Da_fid_lens/Da_fid

;  d_z2=QROMB('distance', z_lens, z_src)/QROMB('distance', 0.0, z_src)
;  d_z1=QROMB('distance', z_lens, z_fid)/QROMB('distance', 0.0, z_fid)

  scaled_alphaX= fltarr(Npix,Npix)
  scaled_alphaY= fltarr(Npix,Npix)

  ratio=f_k_get(1,WHERE(f_k_list(*) eq z_src(0)))
  f_k=FLOAT(ratio(0))
  scaled_alphaX(*,*) = alphaX(*,*)*f_k
  scaled_alphaY(*,*) = alphaY(*,*)*f_k
   PRINT,' f_k ', f_k
   PRINT,' Scaling done '

; = * = * = * = * = * = * = * = * = * = * = * = * = * = * = * 
;Find the image to be delensed for range of i and j
  
  theta_map_r = fltarr(Npix,Npix)
  theta_map_g = fltarr(Npix,Npix)
  theta_map_b = fltarr(Npix,Npix)
  Beta_map_r = fltarr(Npix,Npix)
  Beta_map_g = fltarr(Npix,Npix)
  Beta_map_b = fltarr(Npix,Npix)

  ;HERE SHOULD JUST INCLUDE THE TEMPLATE IMAGE AREA
BetaX=0
BetaY=0

  FOR i=X_min,X_max DO BEGIN
    FOR j=Y_min,Y_max DO BEGIN
	  BetaX = ROUND(i-scaled_alphaX(i,j))
          BetaX=FIX(BetaX(0))
	  BetaY = ROUND(j-scaled_alphaY(i,j))
          BetaY=FIX(BetaY(0))
	  IF r(i,j) gt Beta_map_r(BetaX, BetaY) THEN Beta_map_r(BetaX, BetaY)=r(i,j)
	  theta_map_r(i,j)= r(i,j)
	  IF g(i,j) gt Beta_map_g(BetaX, BetaY) THEN Beta_map_g(BetaX, BetaY)=Beta_map_g(BetaX, BetaY)+g(i,j)
	  theta_map_g(i,j)= g(i,j)
	  IF b(i,j) gt Beta_map_b(BetaX, BetaY) THEN Beta_map_b(BetaX, BetaY)=Beta_map_b(BetaX, BetaY)+b(i,j)
	  theta_map_b(i,j)= b(i,j)
    ENDFOR
  ENDFOR

  PRINT, ' Delens done, start printing the beta map ... '
  FITS_OPEN, src, fcb, /write
  FITS_WRITE,fcb,Beta_map_r,extname='FLUX',EXTVER=1
  FITS_WRITE,fcb,Beta_map_g,extname='FLUX',EXTVER=2
  FITS_WRITE,fcb,Beta_map_b,extname='FLUX',EXTVER=3
  FITS_CLOSE, fcb 

  FITS_OPEN, theta, fcb, /write
  FITS_WRITE,fcb,theta_map_r,extname='FLUX',EXTVER=1
  FITS_WRITE,fcb,theta_map_g,extname='FLUX',EXTVER=2
  FITS_WRITE,fcb,theta_map_b,extname='FLUX',EXTVER=3
  FITS_CLOSE, fcb 
  PRINT,' Delensed Image Printed in File '

; = * = * = * = * = * = * = * = * = * = * = * = * = * = * = * 
;Relens
  PRINT,' Relensing ... '
  ReImg_r=FLTARR(Npix,Npix)
  ReImg_g=FLTARR(Npix,Npix)
  ReImg_b=FLTARR(Npix,Npix)

BetaX=0
BetaY=0

  FOR k=0,Npix-1 DO BEGIN	 
    FOR l=0,Npix-1 DO BEGIN
	  BetaX = ROUND(i-scaled_alphaX(i,j))
          BetaX=FIX(BetaX(0))
	  BetaY = ROUND(j-scaled_alphaY(i,j))
          BetaY=FIX(BetaY(0))
	  ReImg_r(k,l)=ReImg_r(k,l)+Beta_map_r(BetaX,BetaY)
	  ReImg_g(k,l)=ReImg_g(k,l)+Beta_map_g(BetaX,BetaY)
	  ReImg_b(k,l)=ReImg_b(k,l)+Beta_map_b(BetaX,BetaY)
    ENDFOR
  ENDFOR

 
PRINT,' Relens done, start printing the relensed image ... '
  FITS_OPEN, outfile, fcb, /write
  FITS_WRITE,fcb,ReImg_r,extname='FLUX',EXTVER=1
  FITS_WRITE,fcb,ReImg_g,extname='FLUX',EXTVER=2
  FITS_WRITE,fcb,ReImg_b,extname='FLUX',EXTVER=3
  FITS_CLOSE, fcb 
PRINT, ' ' + obj + ' Relensed Image Printed in File '

ENDFOR

END
  
  
   
  
  
