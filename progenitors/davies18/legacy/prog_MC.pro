;; MC-simulation for estimating SN progenitors

;obsfile = 'IIPprog_obsdata_m17av.csv'
obsfile = 'IIPprog_obsdata_2019.csv'
readcol, obsfile, snname, dist, edist, dmod, edmod, mag, emag, band, $
         ebmv, eebmv, av, eav, ai, eai, notes, $
         alam, ealam, BClam, eBClam, $
         L_s09, dL_s09, L_s15, dL_s15, q1, q1, L_us, dL_us,$
         form='(a,f,f,f,f,f,f,a,  f,f,f,f,f,f,a,  f,f,f,f, f,f,f,f, f,f,f,f)', /sil, delim=','
ilfix = where(band eq 'na', nlfix)

iupper = where(notes eq 'upper lim', nupper, compl=idetect)

dum = ''
INCLUDE_UPPERS = 1
dev            = 'x'
ntrials        = 10000


if INCLUDE_UPPERS then begin
   lum_upper  = L_us[iupper]
   name_upper = snname[iupper]
endif else begin
   lum_upper  = L_us[iupper]
   name_upper = snname[iupper]

   snname = snname[idetect]
   dmod   = dmod[idetect]
   edmod  = edmod[idetect]
   mag    = mag[idetect]
   emag   = emag[idetect]
   band   = band[idetect]
   alam   = alam[idetect]
   ealam  = ealam[idetect]
   BClam  = BClam[idetect]
   eBClam = eBClam[idetect]
   L_s09  = L_s09[idetect]
   dL_s09 = dL_s09[idetect]
   L_s15  = L_s15[idetect]
   dL_s15 = dL_s15[idetect]
endelse

;pause   = lambdap("dev: dum='' & if dev eq 'ps' then epsclose else read, dum")


;------- Read-in M-Lfin relation
ML_type = 'STARS'  ;; STARS or KEPLER
case 1 of
   ML_type eq 'STARS' : MLfile = 'M-L_S09.txt'
   ML_type eq 'KEPLER': MLfile = 'M-L_KEPLER.txt'
endcase

readcol, MLfile, mass_mod, lfin_mod, /sil
;mass_mod_int = findgen(18)+9.
;lfin_mod_int = interpol(lfin_mod, mass_mod, mass_mod_int)
;mass_mod = mass_mod_int & lfin_mod = lfin_mod_int


MABS_SUN = 4.74
;LUMS_TO_USE = 'S15 Luminosities'
;LUMS_TO_USE = 'S09 Luminosities'
LUMS_TO_USE = 'updated BCs'    
MASS_09kr = 'not clust'

LUMS_TO_USE_str = strjoin(strsplit(LUMS_TO_USE,' ',/extr),'-')

nsn      = n_elements(snname)
dmod_t   = fltarr(nsn,ntrials)
BC_t     = fltarr(nsn,ntrials)
A_t      = fltarr(nsn,ntrials)
mag_t    = fltarr(nsn,ntrials)
lums_t   = fltarr(nsn,ntrials)
l_s09_t  = fltarr(nsn,ntrials)
l_s15_t  = fltarr(nsn,ntrials)
masses_t = fltarr(nsn,ntrials)


for i=0, nsn-1 do begin
   dmod_t[i,*] = dmod[i] + randomn(seed,ntrials)*edmod[i]
   BC_t[i,*]   = BClam[i] + (randomu(seed,ntrials)-0.5)*eBClam[i]*2.;flat-topped.
   A_t[i,*]    = (alam[i] + randomn(seed,ntrials)*ealam[i]) > 0.
   mag_t[i,*]  = mag[i] + randomn(seed,ntrials)*emag[i]
   l_s09_t[i,*] = l_s09[i] + randomn(seed,ntrials)*dl_s09[i]
   l_s15_t[i,*] = l_s15[i] + randomn(seed,ntrials)*dl_s15[i]
endfor
i_09kr = where(snname eq '09kr')
masses_09kr = randomu(seed,ntrials)*13.+12.


mlo_arr = fltarr(ntrials)
mhi_arr = fltarr(ntrials)


TESTING = 0
if TESTING then plot, [0], [0], xr=[2,40], yr=[-0.1,1.1]


for t=0, ntrials-1 do begin

   ;; Get randomised luminosity
   lums = (mag_t[*,t] - dmod_t[*,t] - A_t[*,t] + BC_t[*,t] - MABS_SUN)/(-2.5)
   lums[ilfix] = l_s15_t[ilfix,t] ;; use S15 values for those with good coverage
   
   ;; Which breed of L to use??
   case 1 of 
      LUMS_TO_USE eq 'S15 Luminosities': lums = l_s15_t[*,t] ;; use S15 
      LUMS_TO_USE eq 'S09 Luminosities': lums = l_s09_t[*,t] ;; use S09 
      LUMS_TO_USE eq 'updated BCs':  ;; use our new updated Ls
   endcase
   lums_t[*,t] = lums

   ;; Turn L into M, using M-L relation. 
   ;; --- If Lcur is beyond the model M-L relation, linearly
   ;;     extrapolate to higher L, assuming linear relation between
   ;;     log(L) and log(M). 
   masses = fltarr(nsn)
   for i=0, nsn-1 do begin
      ;;masses[i] = interpol(mass_mod, lfin_mod, Lums[i])
      masses[i] = 10^(interpol(alog10(mass_mod), lfin_mod, Lums[i]))
      if Lums[i] gt max(lfin_mod) or Lums[i] lt min(lfin_mod) then begin
         masses[i] = 10^(-1.2417 + 0.47868*Lums[i])
         if strpos(MLfile, 'KEPLER') ge 0 then $
            masses[i] = 10^(-1.3688718 + 0.51439852*Lums[i])
      endif
   end
   masses_t[*,t] = masses

   ;; If it's 2009kr, and we are using Justyn's cluster age rather
   ;; than a YSG, sample from a uniform distribution between 12-25Msun. 
   if MASS_09kr eq 'cluster' then masses_t[i_09kr,t] = masses_09kr[t]


   ;; Sort masses in order of increasing mass
   iorder = sort(masses)
   if INCLUDE_UPPERS then masses[iupper] = !values.f_nan
   masses = masses[iorder]

   ;; Fit mass distribution with IMF
   nn0   = findgen(nsn)/(1.*nsn-1.)   ;;Define x-axis (y-axis when plotted)
   M_lo_hi = [8.,20.]
   i_to_fit     = where(finite(masses), nfin)
   nn_tofit     = nn0[i_to_fit] 
   nn_tofit     /= max(nn_tofit)     ;; to rescale from max detected mass. 
   masses_tofit = masses[i_to_fit]

   fit_masses = curvefit(nn_tofit, masses_tofit, $
                         replicate(1.,nfin), M_lo_hi, $
                     function_name='M_hi_lo_fit', /noderiv)

   ;; Store upper/lower mass limits
   mlo_arr[t] = M_lo_hi[0]
   mhi_arr[t] = M_lo_hi[1]

   ;; Testing
   if TESTING then begin
      ydum = findgen(100)/99.
      M_hi_lo_fit, ydum, M_lo_hi, xdum
      oplot, xdum, ydum, col=coln.lblue
      oplot, masses_tofit, nn_tofit, psym=7
      oplot, fit_masses, nn_tofit, psym=5, col=255
      wait, 0.01
   endif
endfor


;; Turn trials into 2D histogram
min1 = 5.5   & max1 = 10.5
min2 = 14  & max2 = 40
bin1 = 0.25 & bin2 = 0.75
hist2d = hist_2d( mlo_arr, mhi_arr, bin1=bin1, bin2=bin2, $
                  min1=min1, max1=max1, min2=min2, max2=max2)
hist2d /= total(hist2d)
minax = findgen((max1-min1)/bin1+1)*bin1+min1
maxax = findgen((max2-min2)/bin2+1)*bin2+min2

;; Get confidence levels
iorder = reverse(sort(hist2d))
cumul = total(hist2d[iorder], /cum)
paxis = findgen(n_elements(hist2d))
pconf = round(interpol(paxis, cumul, [0.68, 0.95, 0.997]))
lconf = hist2d[iorder[pconf]]


;; Plot results
if INCLUDE_UPPERS then begin
   inctxt = '_incUP' 
   title = 'with upper limits'
endif else begin
   inctxt=''
   title = 'detections only'
endelse


setplot, col, dev=dev
loadct, 1, /sil
coln = getcols(dev,/str)
cgsun = cgsymbol('Sun')
eps_setup, 'figs/Mlo-Mhi_'+LUMS_TO_USE_str+inctxt+'_'+ML_type+'.eps', xs=20,ys=18
POSIM = [0.15,0.15,0.97,0.9]
c_col = [coln.red,coln.orange,coln.yellow]
c_col = reverse(c_col)

hist2d_byt = bytscl(hist2d, min=0, max=max(hist2d,/nan))
Xhist2d_byt = hist2d_byt*(-1) + 255
contour, hist2d, minax, maxax, lev=reverse(lconf), $
         xtit='M!dmin!n/M'+cgsun, ytit='M!dmax!n/M'+cgsun, $
         title=title, pos=POSIM,/norm, /fill, $
         c_col=c_col
contour, hist2d, minax, maxax, lev=reverse(lconf), $
         pos=POSIM,/norm, c_col=coln.white, /noerase, /over, $
         c_thick=!p.thick-1, c_line=[1,5,0]

;; error bars
i68 = where(hist2d gt lconf[0])
i95 = where(hist2d gt lconf[1])
i99 = where(hist2d gt lconf[2])
xy68 = array_indices(hist2d, i68)
xy95 = array_indices(hist2d, i95)
xy99 = array_indices(hist2d, i99)
mmin_lo_68 = min( minax[xy68[0,*]] )
mmin_hi_68 = max( minax[xy68[0,*]] )
mmax_lo_68 = min( maxax[xy68[1,*]] )
mmax_hi_68 = max( maxax[xy68[1,*]] )

;; Best fit values
mmax_prof = max(hist2d,dim=1)
mmin_prof = max(hist2d,dim=2)
int_max = lambdap('xx,yy,xmax: ' + $
                  'max = max(yy, imax) & ' +$
                  'xmax = interpol(xx[imax-1:imax+1], (deriv(yy))[imax-1:imax+1], 0.0)')
int_max, minax, mmin_prof, mmin_best
int_max, maxax, mmax_prof, mmax_best
oplot, [mmin_best], [mmax_best], psym=7, th=!p.thick+2, syms=2

legend_astro, ['M!dmin!n/M'+cgsun+' = '+string(mmin_best,'(f3.1)') + $
               '!u+' + string(mmin_hi_68-mmin_best,'(f3.1)') + $
               '!d-' + string(mmin_best-mmin_lo_68,'(f3.1)'), $
               'M!dmax!n/M'+cgsun+' = '+string(mmax_best,'(f4.1)') + $
               '!u+' + string(mmax_hi_68-mmax_best,'(f3.1)') + $
               '!d-' + string(mmax_best-mmax_lo_68,'(f3.1)')], $
              /top, /right, chars=!p.charsize*0.7, spac=2, /clear, clcol=coln.dgrey
if dev eq 'x' then read, dum else epsclose

;; ------- end of MC plot thing -----

;; Produce a table of the input 
form = '(a9, " & ", f5.2, "$\pm$", f4.2, " & ", ' + $
       'f5.2, "$\pm$", f4.2, " & ", ' + $
       'a5, " & ", ' +$
       'f5.2, "$\pm$", f4.2, " & ", ' + $
       'f5.2, "$\pm$", f4.2, " \\ ")'
row_arr = strarr(nsn)

readcol, 'IIPprog_smartt15.csv', ids15, mass_s15, mup_s15, mlo_s15, $
         delim=',', /sil, form='(a,f,f,f)'
ii_s15 = intarr(nsn)
ii_s09 = intarr(nsn)
masses_fit = fltarr(nsn,3)
lums_fit   = fltarr(nsn,2)

snname_long = strarr(nsn)
mass_final  = fltarr(nsn)
for i=0, nsn-1 do begin
   ii_s15[i] = where(ids15 eq snname[i])
   if strmid(snname[i],0,1) eq '9' then prefix = 'SN19' else prefix = 'SN20'
   namei = string(prefix+snname[i],'(a8)')
   snname_long[i] = namei
   dmodout = string(dmod[i],'(f5.2)') + "$\pm$" + string(edmod[i],'(f4.2)')
   magi  = mag[i]
   emagi = emag[i]
   magout = string(mag[i],'(f5.2)') + "$\pm$" + string(emag[i],'(f4.2)')
   bandi = band[i]
   alami = alam[i]
   ealami = ealam[i]
   aout  = string(alam[i],'(f5.2)') + "$\pm$" + string(ealam[i],'(f4.2)')
   bclami = bclam[i]
   ebclami = ebclam[i]
   bcout  = string(bclam[i],'(f5.2)') + "$\pm$" + string(ebclam[i],'(f4.2)')

   ;; ave lum
   lum_ave = median(lums_t[i,*])
   lum_err = stddev(lums_t[i,*])
   if total(i eq ilfix) ge 1 then begin
      lum_ave = L_s15[i]
      lum_err = dL_s15[i]
   endif
   lum_out = string(lum_ave,'(f4.2)') + "$\pm$" + string(lum_err,'(f4.2)')
   lums_fit[i,*] = [lum_ave, lum_err]

   ;; ave mass
   masses_sorted = masses_t[i,sort(masses_t[i,*])]
   mass_ave = masses_sorted[(1.*ntrials*0.5)]
   mass_ave = median(masses_t[i,*])
   mass_final[i] = mass_ave
   mass_ep  = masses_sorted[(1.*ntrials*0.84)] - mass_ave
   mass_em  = mass_ave - masses_sorted[(1.*ntrials*0.16)]
   mass_out = string(mass_ave,'(f4.1)') + $
              "$^{+" + string(mass_ep,'(f3.1)') + '}' + $
              "_{-" + string(mass_em,'(f3.1)') + '}$'

   if total(i eq iupper) gt 0 then begin
      magout   = '$>$'+string(mag[i],'(f5.2)')
      mass_out = '$<$'+string(mass_ave,'(f4.1)') + $
                 "$^{+" + string(mass_ep,'(f3.1)') + '}$'
      lum_out  = '$<$'+string(lum_ave,'(f4.2)') + $
                 "$^{+" + string(lum_err,'(f4.2)') + '}$'
   endif


   masses_fit[i,*] = [mass_ave,mass_ep,mass_em]

   if bandi eq 'na' then begin
      bandi = 'multi'
      aout   = ' -   '
      magout = ' -   '
      bcout  = ' -   '
   endif
   bandi = string(bandi,'(a10)')
   aout   = string(aout,'(a14)')
   magout = string(magout,'(a14)')
   bcout  = string(bcout,'(a14)')


   rowout = strjoin([namei,dmodout,magout,bandi,aout,bcout,lum_out,mass_out], $
                    ' & ') + ' \\ '
   row_arr[i] = rowout

   ;;if snname[i] eq '08bk' then stop

endfor
row_order = sort(strtrim(snname_long,2))
row_arr = row_arr[row_order]
openw, 1, 'figs/prog_MC_tab.tex'
for i=0, n_elements(row_arr)-1 do printf, 1, row_arr[i]
close, 1
spawn, 'cat figs/prog_MC_tab.tex'

;;----------- Lum vs lum
eps_setup, 'figs/lum_comp'+inctxt+'.eps', xs=20, ys=18
plot, [0], [0], xr=[4,6], yr=[4,6], xtit='L/L'+cgsun, $
      ytit='L!dS!n/L'+cgsun, /iso
oplot, !x.crange, !y.crange, li=1
oploterror, L_us, L_s15, $
            dL_us, dL_s15, psym=8, syms=1.5
;oploterror, L_us[*,0], L_s15[ii_s15], $
;            L_us[*,2], L_s15[ii_s15], psym=8, /lo, syms=1.5
if dev eq 'x' then read, dum else epsclose


;;----------- Mass vs Mass
eps_setup, 'figs/mass_comp'+inctxt+'.eps', xs=20, ys=18
plot, [0], [0], xr=[4,26], yr=[4,24], xtit='M!dinit!n/M'+cgsun, $
      ytit='M!dS15!n/M'+cgsun, /iso
oplot, !x.crange, !x.crange, li=1
;oploterror, masses_fit[*,0], mass_s15[ii_s15], $
;            masses_fit[*,1], mup_s15[ii_s15], psym=8, /hi, syms=1.5
;oploterror, masses_fit[*,0], mass_s15[ii_s15], $
;            masses_fit[*,2], mlo_s15[ii_s15], psym=8, /lo, syms=1.5
oploterror, masses_fit[idetect,0], mass_s15[ii_s15[idetect]], $
            masses_fit[idetect,1], mup_s15[ii_s15[idetect]], psym=8, /hi, syms=1.8
oploterror, masses_fit[idetect,0], mass_s15[ii_s15[idetect]], $
            masses_fit[idetect,2], mlo_s15[ii_s15[idetect]], psym=8, /lo, syms=1.8

if INCLUDE_UPPERS then begin
   plotsym, 1, th=!p.thick
   mup_s15[ii_s15[iupper]] = 0.
   oploterror, masses_fit[iupper,0], mass_s15[ii_s15[iupper]], $
               masses_fit[iupper,1], mup_s15[ii_s15[iupper]], psym=8, /hi, syms=2.5, $
               col=col[1], errcol=col[1]
   
   plotsym, 6, th=!p.thick
   oploterror, masses_fit[iupper,0], mass_s15[ii_s15[iupper]], $
               masses_fit[iupper,1], mup_s15[ii_s15[iupper]], psym=8, /hi, syms=2.5, $
               col=col[1], errcol=col[1]
   plotsym, 0, th=!p.thick
endif

;pause, dev
dum = ''
if dev eq 'x' then read, dum

;;----------- Mass vs Jerk Mass
eps_setup, 'figs/mass_comp_jerk'+inctxt+'.eps', xs=20, ys=18
ii_jerk = intarr(nsn)
readcol, 'prog_masses.csv', snjerk, q1, q2, q3, m_j, mp_j, mn_j, $
         form='(a,f,f,f,f,f,f)', /sil
for i=0, nsn-1 do begin
   ii_jerk[i] = where(snjerk eq snname[i])
endfor
plot, [0], [0], xr=[5,25], yr=[5,25], xtit='M!dinit!n/M'+cgsun, $
      ytit='M!dneb!n/M'+cgsun, /iso
oplot, !x.crange, !y.crange, li=1
oploterror, masses_fit[*,0], m_j[ii_jerk], $
            masses_fit[*,1], mp_j[ii_jerk], psym=8, /hi, syms=1.5
oploterror, masses_fit[*,0], m_j[ii_jerk], $
            masses_fit[*,2], mn_j[ii_jerk], psym=8, /lo, syms=1.5
if dev eq 'x' then read, dum else epsclose
;pause, dev


;;---------- Mass spec plot
eps_setup, 'figs/mass_spec'+inctxt+'.eps', xs=24, ys=18

iuse = where(masses_fit[*,0] le max(masses_fit[idetect,0]),nsn)

stop
masses_fit = masses_fit[iuse,*]
snname     = snname[iuse]
notes      = notes[iuse]
iupper = where(notes eq 'upper lim', nupper, compl=idetect)

masses_upp = masses_fit
mso = sort(masses_fit[*,0])

if INCLUDE_UPPERS then begin
   masses_fit[iupper,*]  = !values.f_nan
   masses_upp[idetect,*] = !values.f_nan
endif

plot, [0], [0], yr=[-0.5,nsn-0.5], xr=[3,33], ytickv=indgen(nsn), $
      ytickn=snname[mso], yticks=nsn-1, $
      xtit='M!dinit!n/M'+cgsun, ychars=0.6, $
      ytit='SN progenitor', pos=[0.14,0.14,0.97,0.98]
yax  = findgen(nsn)
yax0 = yax / (nsn-1.)
ey  = replicate(0.,nsn)

;; 3-sig fits
xy_sig = xy99
n_sig  = n_elements(xy_sig[0,*])
x_sig_arr = fltarr(nsn,n_sig)
for n=0, n_sig-1 do begin
   mlo_c = minax[xy_sig[0,n]]
   mhi_c = maxax[xy_sig[1,n]]
   M_hi_lo_fit, yax0, [mlo_c,mhi_c], xax
   x_sig_arr[*,n] = xax
endfor
x_sig_min = min(x_sig_arr,dim=2)
x_sig_max = max(x_sig_arr,dim=2)
polyfill, [x_sig_max,reverse(x_sig_min),x_sig_max[0]], $
          [yax,reverse(yax),yax[0]], col=coln.yellow

;; 2-sig fits
n95  = n_elements(xy95[0,*])
x95_arr = fltarr(nsn,n95)
for n=0, n95-1 do begin
   mlo_c = minax[xy95[0,n]]
   mhi_c = maxax[xy95[1,n]]
   M_hi_lo_fit, yax0, [mlo_c,mhi_c], xax
   x95_arr[*,n] = xax
endfor
x95_min = min(x95_arr,dim=2)
x95_max = max(x95_arr,dim=2)
polyfill, [x95_max,reverse(x95_min),x95_max[0]], $
          [yax,reverse(yax),yax[0]], col=coln.orange

;; 1-sig fits
n68  = n_elements(xy68[0,*])
x68_arr = fltarr(nsn,n68)
for n=0, n68-1 do begin
   mlo_c = minax[xy68[0,n]]
   mhi_c = maxax[xy68[1,n]]
   M_hi_lo_fit, yax0, [mlo_c,mhi_c], xax
   x68_arr[*,n] = xax
endfor
x68_min = min(x68_arr,dim=2)
x68_max = max(x68_arr,dim=2)
polyfill, [x68_max,reverse(x68_min),x68_max[0]], $
          [yax,reverse(yax),yax[0]], col=coln.red

;; Best fit
M_hi_lo_fit, yax0, [mmin_best,mmax_best], xax
oplot, xax, yax, col=coln.blue, th=!p.thick+2, li=0

;; Fit when Mup=30Msun
M_hi_lo_fit, yax0, [mmin_best,30.], xax
oplot, xax, yax, col=coln.blue, th=!p.thick+2, li=5

;;Data
plotsym, 0, th=!p.thick
oploterror, masses_fit[mso,0], yax, masses_fit[mso,1], ey, psym=8, /hi, syms=1.5
oploterror, masses_fit[mso,0], yax, masses_fit[mso,2], ey, psym=8, /lo, syms=1.5

if INCLUDE_UPPERS then begin
   plotsym, 6, th=!p.thick
   oploterror, masses_upp[mso,0], yax, masses_upp[mso,1], ey, psym=8, /hi, syms=1.8
   oplot, masses_upp[mso,0], yax, psym=1
endif
if dev eq 'x' then read, dum else epsclose
;--------------------------


;;---------- Generic mass spec plot ------
eps_setup, 'figs/mass_spec_koch.eps', xs=22, ys=18

plot, [0], [0], yr=[-0.5,nsn-0.5], xr=[3,25], ytickv=indgen(nsn), $
      ytickn=snname[mso], yticks=nsn-1, $
      xtit='M!dinit!n/M'+cgsun, ychars=0.6, $
      ytit='SN progenitor', pos=[0.14,0.14,0.97,0.98]

;;Data
plotsym, 0, th=!p.thick
oploterror, masses_fit[mso,0], yax, masses_fit[mso,1], ey, psym=8, /hi, syms=1.5
oploterror, masses_fit[mso,0], yax, masses_fit[mso,2], ey, psym=8, /lo, syms=1.5

if INCLUDE_UPPERS then begin
   plotsym, 6, th=!p.thick
   oploterror, masses_upp[mso,0], yax, masses_upp[mso,1], ey, psym=8, /hi, syms=1.8
   oplot, masses_upp[mso,0], yax, psym=1
endif

nfin = float(total(finite(masses_fit[mso,0])))

;; Kochanek's fit 
M_hi_lo_fit, yax0, [6.3,15.7], xax_koch
oplot, xax_koch, yax, col=coln.orange, th=!p.thick+2, li=2

;; chisq of Kochanek's fit
diffs_koch = (xax_koch - masses_fit[mso,0])/masses_fit[mso,2]
chisq_koch = total(diffs_koch^2,/nan)
prob_koch  = exp(-0.5*chisq_koch/nfin)

;; chisq of our fit
M_hi_lo_fit, yax0, [7.5,19.0], xax_us
diffs = (xax_us - masses_fit[mso,0])/masses_fit[mso,2]
diffs_poz = (xax_us - masses_fit[mso,0])/masses_fit[mso,1]
ipoz = where(diffs gt 0.)
diffs[ipoz] = diffs_poz[ipoz]
chisq_us = total(diffs^2,/nan)
oplot, xax_us, yax, col=coln.blue, th=!p.thick+2, li=5
prob_us  = exp(-0.5*chisq_us/nfin)

print, 'Koch prob: ',(1.-prob_koch)*100.

;; Us, 2020
M_hi_lo_fit, yax0, [7.1,17.9], xax_us20
oplot, xax_us20, yax, col=coln.blue, th=!p.thick+2, li=1



legend_astro, ['Kochanek 2020','Davies & Beasor 2018','Davies & Beasor 2020'], $
              li=[2,5,1], col=[coln.orange,coln.blue,coln.blue], $
              /bot, /right, chars=!p.charsize*0.7, spac=2




prob_diff = exp(-0.5* (chisq_koch-chisq_us))


;; MC the errors
nmc = 10000
xax_t = fltarr(nsn,nmc)
masses_comb = masses_fit
masses_comb[where(~finite(masses_fit))] = masses_upp[where(~finite(masses_fit))]
for i=0, nmc-1 do begin
   for n=0, nsn-1 do begin
      nn = mso[n]
      xax_t[nn,i] = randomn(seed,1)*masses_comb[nn,1] + xax[nn]
      
   endfor
   xax_t[*,i] = xax_t[sort(xax_t[*,i]),i]
endfor
;oplot, median(xax_t,dim=2), yax
stop

if dev eq 'x' then read, dum else epsclose

;;-----------------------------

;;--- Upper limits
mass_upper = fltarr(nupper)
out_str_arr = strarr(nupper)
for i=0, nupper-1 do begin
   mass_upper[i] = interpol(mass_mod, lfin_mod, Lum_upper[i])
   up_name = 'SN ' + name_upper[i]
   up_mass = '<' + string(mass_upper[i],'(f4.1)')
   out_str = strjoin([up_name, up_mass], ' & ') + '\\'
   out_str_arr[i] = out_str
end
uo = sort(mass_upper)
out_str_arr = out_str_arr[uo]
print, ' '
print, out_str_arr

for i=0, nsn-1 do begin
   print, snname[i], masses_fit[i,*], f='(a6, 3(2x,f4.1))'
endfor
for i=0, nupper-1 do begin
   print, name_upper[i], mass_upper[i], 0.0, 0.0, f='(a6, 3(2x,f4.1))'
endfor



end
