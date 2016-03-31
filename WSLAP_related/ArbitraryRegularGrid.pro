PRO ArbitraryRegularGrid
;
;
; This program is originally designed by Chema
; In usual WSLAP configurations, we commonly start with 1024 grid points i.e. 32x32
; This program allows us to use any regular dimension of grid points
; To use this program, replace the zero_mass_grid.dat in WSLAP with the file produced by this program, and run WSLAP after Make Grid
;
;


CLOSE,1

N = 512
Ng = 38  ; Makes a NgxNg grid

CellSize = 511.0/FLOAT(Ng)

OPENW,1,'Mass_cell_multiR_ArbitraryRegGrid_38x38.dat'
CS = ROUND(CellSize)
CC = 0
Mass = 1.0E-2
X2i = 1
map = FLTARR(N,N)

FOR i=0,Ng-1 DO BEGIN
    X1i = X2i 
    X2i = X1i + CellSize
    X2j = 1
    FOR j=0,Ng-1 DO BEGIN
        CC = CC + 1
        X1j = X2j 
        X2j = X1j + CellSize
        printf,1,FORMAT='(I4,E13.2,3I5)',$
              CC,Mass,ROUND(X1i),ROUND(X1j),CS

        ix = ROUND(X1i+CS/2.)
        iy = ROUND(X1j+CS/2.)

        map(ix-3:ix+3,iy) = -1.0
        map(ix,iy-3:iy+3) = -1.0

    ENDFOR
ENDFOR
CLOSE,1

print,' Ncells ',Ng^2,CC

WINDOW,XS=N,YS=N
TVSCL, map

END
