PRO PIEMD_model


                                ;Mandy
                                ;2016

                                ;This is for creating a mass map with PIEMD profile
                                ;the normalisation is arbitrary
  





  

q = 0.6     ;axis ratio
e = (1.0-q)/(1.0+q)
A=1.0
rc=30.0     ;core size
mass=fltarr(1600,1600)

theta=(150.0/180.0)*!pi   ;position angle
for i=0,1599 do begin
  for j=0,1599 do begin
    x=(i-799)*cos(theta)-(j-799)*sin(theta)
    y=(i-799)*sin(theta)+(j-799)*cos(theta)

    D=SQRT(x^2/(1+e)^2+y^2/(1-e)^2+rc^2)
    mass[i,j]=A/D
  endfor
endfor
;mass[799,799]=0.28
;you need to assign a mass value to this pixel is rc=0
writefits,'/Users/Mandy/HiRes_Light/BCG_TreuGrid_ellipticalPIEMD/BCGmassPIEMD_elliptical_v7_ellip0p6_-60deg.fits',mass

END
