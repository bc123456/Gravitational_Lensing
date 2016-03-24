pro Abel_Transform
  
                                ;Mandy
                                ;2016

                                ;De-project 2D mass profile to 3D spherical mass distribution
                                ;with inverse Abel Transformation
  
  
G = 6.67E-8
M_sun = 1.99E33
pix_ratio = 3.08568E21*6.4*0.03   ;in cm




mass_map_2d = readfits('/Users/Mandy/FiducialMass/mass_BCG_datalight_v2_FastCompute_Smooth3.fits')
writefits,'~/HiRes_Light/KickVelocity/2DMass_FastCompute_Smooth3.fits',mass_map_2d
x0 = 799
y0 = 799      ;centre = (800,800)
  

;=========STEP 1: GET THE 2D DENSITY MAP=====================

flux = fltarr(400)                         
density_2d = fltarr(350)    
area = fltarr(350)                           
flux[0] = mass_map_2d[x0,y0]                ;This is the first
area[0] = 1.0                               ;pixel - origin
density_2d[0] = mass_map_2d[x0,y0]



for r = 1,349 do begin
  pixel = 0
  
  for x = x0-r,x0+r do begin
    for y = y0-r,y0+r do begin
      ;print,'r',r
      ;print,'x,y',x,y
      dx = float(x-x0)
      dy = float(y-y0)
      radius = SQRT(dx*dx+dy*dy)
      ;print,'radius',radius
      if radius le float(r) and radius gt float(r-1) then begin
        pixel = pixel+1
        flux[r] = flux[r] + mass_map_2d[x,y]
      endif 
      ;if r gt 300 then begin
        ;print,'r',r
        ;print,'x,y',x,y
        ;print,'x-x0,y-y0',dx,dy
        ;print,'radius',radius
      ;endif
      ;print,r
      ;print,radius 
    endfor
  endfor
  
  flux[r] = flux[r]
  area[r] = pixel
  density_2d[r] = (flux[r]/FLOAT(area[r]))
  ;print,pixel
endfor

;print,density_2d[349]
;print,'area',area[0:3]
;print,'density_2d',density_2d[0:3]
;print,'flux',flux[300:349]
;print,'area',area[300:349]
;cgPlot,indgen(350),density_2d,xrange=[0,349],yrange=[0,10]
;print,'density_2d',density_2d[300:349]







;Take rn = 349  
;
;=========STEP 2: Calculate 2D density gradient===========


density_grad_2d = fltarr(349)

density_grad_2d[0] = density_2d[1]-density_2d[0] ;the first point

for r = 1,348 do begin
  density_grad_2d[r] = (density_2d[r+1]-density_2d[r-1])/2.0
endfor

;cgplot,indgen(349),density_grad_2d,/NoErase
print,density_grad_2d

openw,lun,'~/HiRes_Light/KickVelocity/MassDensity2D_RadialGradient_FastCompute_Smooth3.dat',/get_lun
printf,lun,density_grad_2d
free_lun,lun





;=========STEP 3: Abel inverison transform===========

rn = 348
dr = -1.0
rho = dblarr(348)


for r =0,347 do begin

  for i =0,rn-r-1  do begin
    R_index = rn-i
    R_0 = float(R_index)
    rr_0 = float(r)
    a = SQRT(float(R_0^2-rr_0^2))
    z = (1.0/a)*density_grad_2d[R_index]
    f = z*dr
    rho[r] = (rho[r]+f)
    ;print,'R_index',R_index
    ;print,'a',a
    ;print,'density_grad_2d',density_grad_2d[R_index]
    ;print,'z',z
    ;print,'f',f
    ;print,'rho',rho[r]
  endfor
  
endfor

rho = rho/!pi
print,rho[290:300]

openw,lun,'~/HiRes_Light/KickVelocity/MassDensity_3D_FastCompute_Smooth3.dat',/get_lun
printf,lun,rho
free_lun,lun
print,size(density_2d)
print,size(density_grad_2d)
print,size(rho)
;cgplot,rho,/NoErase
x1=findgen(150)
cgDisplay, 1200, 550
cgPlot, x1, density_2d[0:149],$
  /xlog,/ylog,$
  xrange = [0.1,150],$
;  yrange = [min(density_2d),max(density_2d)+10],$
  Position=[0.10, 0.10, 0.35, 0.90]
cgPlot, findgen(150), density_grad_2d[0:149],Position=[0.4, 0.10, 0.65, 0.90], /NoErase
cgPlot, findgen(150)*0.03, rho[0:149], /xlog,/ylog, xrange=[0.01,150],Position=[0.70, 0.10, 0.95, 0.90], /NoErase
;cgAxis,/XLOG,/YLOG,/Window










end
