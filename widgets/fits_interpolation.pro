PRO fits_interpolation

;
;
; Quickly interpolate a fits file to a larger file with dimension Npix
;
;


infile = '/Volumes/BRIAN/Research/WSLAP_results/initial_error/mass_sd.fits'
outfile = '/Volumes/BRIAN/Research/WSLAP_results/initial_error/mass_sd_full.fits'
;infile = '/Volumes/BRIAN/Research/WSLAP_works/data/M0647/ver22/reconstruction/critcurve_map_ver22_z3.0.fits'
;outfile = '/Volumes/BRIAN/Research/WSLAP_works/data/M0647/ver22/reconstruction/critcurve_map_ver22_z3.0_full.fits'
;infile = '/Volumes/BRIAN/Research/WSLAP_works/data/M0647/ver23/reconstruction/recomp_alpha_x_rad_ver23.fits'
;outfile = '/Volumes/BRIAN/Research/WSLAP_works/data/M0647/ver23/reconstruction/recomp_alpha_x_rad_ver23_full.fits'


Npix = 2048

original = READFITS(infile)
S = SIZE(original)
Npix0 = S(1)

X_ACS = FLOAT(Npix0)*FINDGEN(Npix)/FLOAT(Npix-1)
Y_ACS = FLOAT(Npix0)*FINDGEN(Npix)/FLOAT(Npix-1)

final = INTERPOLATE(original,X_ACS,Y_ACS,/GRID)
PRINT,' Interpolation done '

WRITEFITS,outfile,final

END
