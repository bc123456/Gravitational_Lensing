PRO ArbitraryRegularGrid_v2
;
; 
; This program is originally designed by Chema
; In usual WSLAP configurations, we commonly start with 1024 grid points i.e. 32x32
; This program creates an irregular grid with denser grid points at the middle and fewer grid points at outer region
; *Starts at the centre with a given size and works out wth increaseon sizes
; To use this program, replace the zero_mass_grid.dat in WSLAP with the file produced by this program, and run WSLAP after Make Grid
;
;


CLOSE,1,2

N=512
map = FLTARR(512,512)

Mass = 0.0

DoNextCell = 1

CellSize = 9.0  & StepC = 2.5 ; 280 cells
OPENW,1,'Mass_cell_multiR_10arcmin_V4.dat'

;OPENW,2,'XY_Mandy.dat'

;CellSize = 4.0 & StepC = 1.5  ; 576 cells. High Res
;OPENW,1,'/home/brian/WSLAP_works2/data/Mass_cell_multiR_10arcmin_V4.dat'

;CellSize = 3.0 & StepC = 1.0  ; Mannually set the parameters
;OPENW,1,'/home/brian/WSLAP_works2/data/Mass_cell_multiR_10arcmin_V4.dat'


N_rows   = 2-1
Length   = 0.0
CC       = 0 
WHILE (DoNextCell EQ 1) DO BEGIN

  CellSize = CellSize + StepC
  Length   = Length + 2.0*CellSize
  N_rows   = ROUND(Length/CellSize)
 
  CellSize_True = FLOAT(Length)/FLOAT(N_rows)

  PRINT,CellSize,CellSize_True,Length,N_rows

  IF (Length LT 515.0) THEN BEGIN
     Xi = ROUND(256.0 + (1.0+FINDGEN(N_rows) - N_rows/2.0 - 0.5)*CellSize_True)
     Yi = Xi
    
     FOR i=0,N_rows-1 DO BEGIN  ; Bottom row
         CC = CC + 1
         x=ROUND(Xi(i)-CellSize_True/2.0)
         y=ROUND(Yi(0)-CellSize_True/2.0)
         printf,1,FORMAT='(I4,E13.2,2I6,1F8.3)',$
                CC,Mass,x,y,CellSize_True
;printf,2,x-256+800,y-256+800

         x = x + CellSize_True/2.0
         y = y + CellSize_True/2.0
         map(x-3:x+3,y)=-1.0
         map(x,y-3:y+3)=-1.0
     ENDFOR
     FOR i=0,N_rows-1 DO BEGIN  ; Top row
         CC = CC + 1
         x=ROUND(Xi(i)-CellSize_True/2.0)
         y=ROUND(Yi(N_rows-1)-CellSize_True/2.0)
         printf,1,FORMAT='(I4,E13.2,2I6,1F8.3)',$
                CC,Mass,x,y,CellSize_True
;printf,2,x-256+800,y-256+800

         x = x + CellSize_True/2.0
         y = y + CellSize_True/2.0
         map(x-3:x+3,y)=-1.0
         map(x,y-3:y+3)=-1.0
     ENDFOR
     FOR i=1,N_rows-2 DO BEGIN  ; Left column
         CC = CC + 1
         x=ROUND(Xi(0)-CellSize_True/2.0)
         y=ROUND(Yi(i)-CellSize_True/2.0)
         printf,1,FORMAT='(I4,E13.2,2I6,1F8.3)',$
                CC,Mass,x,y,CellSize_True
;printf,2,x-256+800,y-256+800

         x = x + CellSize_True/2.0
         y = y + CellSize_True/2.0
         map(x-3:x+3,y)=-1.0
         map(x,y-3:y+3)=-1.0
     ENDFOR
     FOR i=1,N_rows-2 DO BEGIN  ; Right column
         CC = CC + 1
         x=ROUND(Xi(N_rows-1)-CellSize_True/2.0)
         y=ROUND(Yi(i)-CellSize_True/2.0)
         printf,1,FORMAT='(I4,E13.2,2I6,1F8.3)',$
                CC,Mass,x,y,CellSize_True
;printf,2,x-256+800,y-256+800

         x = x + CellSize_True/2.0
         y = y + CellSize_True/2.0
         map(x-3:x+3,y)=-1.0
         map(x,y-3:y+3)=-1.0
     ENDFOR
  ENDIF
  IF (Length GE 515.0) THEN DoNextCell=0
ENDWHILE
CLOSE,1

map(0:7,*) = -1
map(*,0:7) = -1
map(N-8:N-1,*) = -1
map(*,N-8:N-1) = -1

WINDOW,1,XS=512,YS=512,RETAIN=2
TVSCL, map
print,' Ncells ',CC

OPENW,1,'N_cells.dat'
PRINTF,1,CC
PRINTF,1,CC
PRINTF,1,0
PRINTF,1,-1
CLOSE,1

END
