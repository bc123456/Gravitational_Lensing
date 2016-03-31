PRO data_converter

;
;
; This program changes a .dat file in WSLAP to a .fits file, or reverse
; Mainly used to work on the mass_map.dat
;
;

;setting the mode of the program

mode = 2       ;mode = 1 when convert data to fits, mode = 2 when convert fits to data

;Setting initial parameters

;parameter for all two modes
x = 512.0    ;x-dimension of the file
y = 512.0    ;y-dimension of the file
data = '/home/brian/WSLAP_works/data/M0647/ver20/fiducial/layer7.dat'      ;directory and name of the data file
fits = '/home/brian/WSLAP_works/data/M0647/ver20/fiducial/layer7.fits'      ;directory and name of the fits file


;parameters for mode 1 only, change if necessary
n = 5.0    ;total number of layer
n = float(n)

;parameters for mode 2 only, change if necessary
l = 7      ;the layer you work at
m = 2.5E5  ;the mass to light ratio.

;convert data to fits

IF mode EQ 1 THEN BEGIN
  print, 'Converting .dat to .fits'
  
  openr,lun,data,/get_lun
  k = x*y*n
  datam = fltarr(4,k)
  readf,lun,datam
  
  fitsm = fltarr(x,y,n)
  
  for i = 0, ((x*y*n)-1) do begin &$
    x1 = datam(0,i) - 1
    y1 = datam(1,i) - 1 
    n1 = datam(3,i) - 1
    fitsm(x1,y1,n1) = datam(2,i)
  endfor
  writefits,fits,fitsm
ENDIF
IF mode EQ 2 THEN BEGIN
  
  print, 'Converting .fits to .dat'
  
  layer_fits = readfits(fits)
  layer_fits = layer_fits/m
  size = size(layer_fits)
  x = size[1]
  y = size[2]
  layer_dat = fltarr(4,x*y)
  
  ;construct the (x,y) coordinate column
  for i = 1,x do begin &$
    for j = 1,y do begin &$
    layer_dat(0, ((i-1)*y + (j-1))) = fix(i) &$
    layer_dat(1, ((i-1)*y + (j-1))) = fix(j) &$
    endfor
  endfor
  
  ;construct the third column (value column)
  for i = 1, (x*y) do begin &$
    layer_dat(2,(i-1)) = (layer_fits((layer_dat(0,(i-1)) - 1), (layer_dat(1,(i-1)) - 1))) &$
  endfor
  
  ;construct the fourth column (layer)
  layer_dat(3,*) = l
  
  ;write to file
  openw,lun,data,/get_lun
  printf, lun, layer_dat, format='(I,I,F,I)'
  free_lun, lun
  
ENDIF ELSE BEGIN
IF mode NE 1 AND mode NE 2 THEN $
  print, 'error, please check the mode'
ENDELSE

END
