;;
;;------ v3:
;;       --: removed all the 'cumulative' stuff
;;       --: changed to chisq from P values in tri-plots.
;;       --: adding in synthetic obs, for testing purposes
;;       --: removed all plotting code

function wmean, vals, weights

   dw = total(1D * weights)
   sum = total(1D * vals * weights) / dw

   return, sum

end

function opthist, var, binmin, binmax, binsize

  outhist = histogram(var, binsize=binsize,min=binmin,max=binmax)

  return, outhist

end

function outxhist, binmin, binmax, binsize

  xarray = findgen((binmax-binmin)/binsize)/((binmax-binmin)/binsize) * $
             (binmax-binmin)+binmin + binsize/2.

  return, xarray

end

pro create_synth_obs, NTRIALS=NTRIALS, LSTEP=LSTEP, LMIN=LMIN, LMAX=LMAX, $
                     LPOBS_SYNTH=LPOBS_SYNTH, NSN_OBS=NSN_OBS, OBS_STR=OBS_STR, $
                     LLO_IN=LLO_IN, LHI_IN=LHI_IN, GAM_IN=GAM_IN, NSN_OUT=NSN_OUT

  if ~keyword_Set(NSN_OBS) then NSN_OBS = 24
  if ~keyword_Set(NSN_OUT) then NSN_OUT = 24
  if ~keyword_Set(NTRIALS) then NTRIALS = 10000
  if ~keyword_Set(LSTEP)   then LSTEP   = 0.02
  if ~keyword_Set(LMIN)   then LMIN   = 4.0
  if ~keyword_Set(LMAX)   then LMAX   = 5.8
  if ~keyword_Set(NSN_OBS) then NSN_OBS = 24
  if ~keyword_Set(LLO_IN) then LLO_IN = 4.3
  if ~keyword_Set(LHI_IN) then LHI_IN = 5.2
  if ~keyword_Set(GAM_IN) then GAM_IN = -1.0

  message, 'Model params {' + $
           strjoin(string([LLO_IN,LHI_IN,GAM_IN],'(f5.2)'),',') + '}', $
           /info

  nsn = NSN_OBS

  lhist = outxhist(lmin,lmax,lstep)
  nl    = n_elements(lhist)
  parr  = dblarr(nl,nsn)
  tarr  = fltarr(NTRIALS,nsn)

  nsn_real = 1.*(n_elements(obs_str))

  for t=0, NTRIALS-1 do begin

     ;; --- method 2: take real data, add in more if necessary
     if nsn le nsn_real then begin
        randi = round(randomu(seed,nsn) * (nsn_real-1.))
        lum_t = obs_str[randi].L_us
        lum_t = lum_t + randomn(seed,nsn)*obs_str[randi].DL_us
        tarr[t,*] = lum_t[sort(lum_t)]
     endif else begin
        lum_t = fltarr(nsn)
        nsn_ex = nsn - nsn_real
        lum_t[0:(nsn_real-1)] = obs_str.L_us + $
                                randomn(seed,nsn_real)*obs_str.DL_us
        randomp, lum_t2, GAM_IN, nsn_ex, range=[10^LLO_IN, 10^LHI_IN]
        lum_t2 = alog10(lum_t2)
        lum_t2 = lum_t2[sort(lum_t2)]
        dL = interpol(obs_str.DL_us, findgen(nsn_real)/nsn_real, $
                      findgen(nsn_ex)/(1.*nsn_ex) )
        Li = lum_t2 + randomn(seed,nsn_ex)*dL
        Li_sort = Li[sort(Li)]
        lum_t[nsn_real:*] = Li_sort
        tarr[t,*] = lum_t[sort(lum_t)]
     endelse
  endfor

  for i=0, nsn-1 do begin
     nhist = double(opthist(tarr[*,i],lmin,lmax-lstep,lstep))
     nhist /= 1.*total(nhist)
     parr[*,i] = nhist
  endfor

  LPOBS_SYNTH = parr

  yin    = findgen(NSN_OBS )
  yout   = findgen(NSN_REAL)/float(NSN_REAL)*float(NSN_OBS)

end

;=======================================================================
;=======================================================================


pro create_Lobs, NTRIALS=NTRIALS, LSTEP=LSTEP, LMIN=LMIN, LMAX=LMAX, $
                 LPOBS=LPOBS, NSN_OBS=NSN_OBS, OBS_STR=OBS_STR

;; Create the observed P-distribution of Ls

obsfile = 'IIPprog_obsdata_2022.csv'
readcol, obsfile, snname, dist, edist, dmod, edmod, mag, emag, band, $
         ebmv, eebmv, av, eav, ai, eai, notes, $
         alam, ealam, BClam, eBClam, $
         L_s09, dL_s09, L_s15, dL_s15, q1, q1, L_us, dL_us,$
         form='(a,f,f,f,f,f,f,a,  f,f,f,f,f,f,a,  f,f,f,f, f,f,f,f, f,f,f,f)', /sil, delim=','
nsn = n_elements(snname)

iupper = where(notes eq 'upper lim', nupper, compl=idetect)
iuse = where(L_us le max(L_us[idetect]), nuse)
nsn = nuse
L_us = L_us[iuse]
dL_us = dL_us[iuse]
notes = notes[iuse]
snname = snname[iuse]
dmod = dmod[iuse]
edmod = edmod[iuse]
mag = mag[iuse]
emag = emag[iuse]
alam = alam[iuse]
ealam = ealam[iuse]
BClam = BClam[iuse]
eBClam = eBClam[iuse]
iupper = where(notes eq 'upper lim', nupper, compl=idetect)

s       = {snname:'', dmod:0., edmod:0., mag:0., emag:0., $
           alam:0., ealam:0., BClam:0., eBClam:0., $
           L_us:0., dL_us:0., notes:''}
obs_str = replicate(s,nsn)
obs_str.snname = snname
obs_str.dmod   = dmod
obs_str.edmod  = edmod
obs_str.mag    = mag
obs_str.emag   = emag
obs_str.alam   = alam
obs_str.ealam  = ealam
obs_str.BClam  = BClam
obs_str.eBClam = eBClam
obs_str.L_us   = L_us
obs_str.dL_us  = dL_us
obs_str.notes  = notes

if ~keyword_Set(NTRIALS) then NTRIALS = 10000
if ~keyword_Set(LSTEP)   then LSTEP   = 0.01
if ~keyword_Set(LMIN)   then LMIN   = 4.0
if ~keyword_Set(LMAX)   then LMAX   = 5.8

lhist = outxhist(lmin,lmax,lstep) ;;; ?
nl    = n_elements(lhist)
parr  = dblarr(nl,nsn)
tarr  = fltarr(NTRIALS,nsn)

for t=0, NTRIALS-1 do begin
   Li = fltarr(nsn)
   for i=0, nsn-1 do begin
      Li[i] = ( (mag[i]+randomn(seed,1)*emag[i]) - $
                (dmod[i]+randomn(seed,1)*edmod[i]) - $
                ( (alam[i]+randomn(seed,1)*ealam[i]) > 0.) + $
                (BClam[i]+randomn(seed,1)*eBClam[i]) - $
                4.74) / (-2.5)
      if mag[i] lt 0. then Li[i] = L_us[i] + randomn(seed,1)*dL_us[i]
      Li[i] = L_us[i] + randomn(seed,1)*dL_us[i] ;;testing
   endfor
   Li_nan  = Li
   Li_nan[iupper] = !values.f_nan
   Li_sort = Li_nan[sort(Li)]
   tarr[t,*] = Li_sort
endfor

print, lmin, lmax, lstep

for i=0, nsn-1 do begin
   nhist = double(opthist(tarr[*,i],lmin,lmax-lstep,lstep))  ;;; ?
   nhist /= 1.*total(nhist)
   parr[*,i] = nhist
endfor

parr_n = parr
for i=0, nsn-1 do begin
   prof = total(parr[*,i],/cum)
   prof /= max(prof)
   parr_n[*,i] = parr[*,i] / max(parr[*,i])
   lmed = interpol(lhist,prof,0.50)
   elo = lmed - interpol(lhist,prof,0.16)
   ehi = interpol(lhist,prof,0.16) - lmed
endfor

LPOBS = parr
NSN_OBS = nsn

for i=0, nsn-1 do begin
   ii = Li[i]
   yval = [i+1.1]/float(nsn)
endfor

save, LPOBS, file='LPOBS.sav'

end

;==============================================

pro create_modelgrid,  LSTEP=LSTEP, LMIN=LMIN, LMAX=LMAX, $
                       NSN_OBS=NSN_OBS, LPMOD=LPMOD, NTRIALS=NTRIALS, $
                       NLLO=NLLO, NLHI=NLHI, NGAM=NGAM, OBS_STR=OBS_STR, $
                       LLO_RANGE=LLO_RANGE, LHI_RANGE=LHI_RANGE, GAM_RANGE=GAM_RANGE, $
                       NO_RESULTS=NO_RESULTS, $
                       TAPER_LLO=TAPER_LLO, TAPER_LHI=TAPER_LHI, $
                       OUTNAME=OUTNAME

  if ~keyword_Set(NTRIALS) then NTRIALS = 10000
  if ~keyword_Set(LSTEP)   then LSTEP   = 0.01
  if ~keyword_Set(LMIN)    then LMIN    = 4.45
  if ~keyword_Set(LMAX)    then LMAX    = 5.28
  if ~keyword_Set(NSN_OBS) then NSN_OBS = 24
  if ~keyword_Set(NLLO)    then NLLO    = 20
  if ~keyword_Set(NLHI)    then NLHI    = 19
  if ~keyword_Set(NGAM)    then NGAM    = 10
  if ~keyword_Set(TAPER_LLO) then TAPER_LLO    = 0
  if ~keyword_Set(TAPER_LHI) then TAPER_LHI    = 0
  if ~keyword_Set(OUTNAME)   then OUTNAME      = 'LPMOD'

  if ~keyword_set(LLO_RANGE) then begin
     llo_min = 4.1 & llo_max = 4.5
  endif else begin
     llo_min = LLO_RANGE[0] & llo_max = LLO_RANGE[1]
  endelse

  if ~keyword_set(LHI_RANGE) then begin
     lhi_min = 5.1 & lhi_max = 5.7
  endif else begin
     lhi_min = LHI_RANGE[0] & lhi_max = LHI_RANGE[1]
  endelse

  if ~keyword_set(GAM_RANGE) then begin
     gam_min = 0.0 & gam_max = 2.0
  endif else begin
     gam_min = GAM_RANGE[0] & gam_max = GAM_RANGE[1]
  endelse

  nsn_real = n_elements(obs_str)
  nsn = nsn_obs
  nllo = NLLO & nlhi = NLHI & npl = NGAM
  n_array_elems = nllo*nlhi*npl


  llo = findgen(nllo)/(nllo-1.) * (llo_max-llo_min) + llo_min
  lhi = findgen(nlhi)/(nlhi-1.) * (lhi_max-lhi_min) + lhi_min
  pl  = findgen(npl)/(npl-1.) * (gam_max-gam_min) + gam_min

  if npl eq 1 then pl = [( - 1.7)]
  loglik = fltarr(nllo,nlhi) * !values.f_nan
  chisqarr  = fltarr(nllo,nlhi) * !values.f_nan
  cumspecs  = fltarr(nllo,nlhi,nsn)

  lhmin = LMIN & lhmax=LMAX & lhstep=LSTEP
  lxhist = outxhist(lhmin,lhmax,lhstep)
  nl     = n_elements(lxhist)

  ;; Random Lums
  results = ptrarr(nllo,nlhi,npl)
  counter = 0
  pc_done = 0

  if keyword_set(no_results) then begin
     message, 'NO_RESULTS is set - skipping computation of cumulative L', /info
     sav_extension = '_void'
     goto, SKIPRESULTS
  endif else sav_extension = ''

  if TAPER_LLO ne 0. then sav_extension += '_TAPER'

  message, 'Array elements: '+strn(n_array_elems), /info

  print, llo
  print, lhi
  print, pl

  time0 = systime(/sec)
  calctime = 1

  for xx=0,nllo-1 do begin
     for yy=0,nlhi-1 do begin
        for pp=0,npl-1 do begin

           LRAN_MIN = 10^llo[xx]
           LRAN_MAX = 10^lhi[yy]
           PLCUR    = pl[pp]
           if TAPER_LLO ne 0. then LRAN_MIN = 10^(llo[xx]-TAPER_LLO)
           if TAPER_LHI ne 0. then LRAN_MAX = 10^(lhi[yy]+TAPER_LHI)
           randomp, lums_mast, PLCUR, 1e6, range=[LRAN_MIN,LRAN_MAX]

           if TAPER_LLO ne 0. then begin
              iTAPER = where(lums_mast lt 10^(llo[xx]+TAPER_LLO), nTAPER, comp=iptaper)
              if nTAPER gt 0 then begin
                 tprob  = (alog10(lums_mast[iTAPER]) - alog10(LRAN_MIN) ) / (2.*TAPER_LLO)
                 trand  = randomu(seed,nTAPER)
                 ikeep  = where(tprob gt trand)
                 lout   = [lums_mast[iptaper], lums_mast[iTAPER[ikeep]] ]
                 lums_mast = lout
              endif
           endif

           if TAPER_LHI ne 0. then begin
              iTAPER = where(lums_mast gt 10^(lhi[yy]-TAPER_LHI), nTAPER, comp=iptaper)
              if nTAPER gt 0 then begin
                 tprob  = (alog10(lums_mast[iTAPER]) - alog10(LRAN_MAX) ) / (-2.*TAPER_LHI)
                 trand  = randomu(seed,nTAPER)
                 ikeep  = where(tprob gt trand)
                 lout   = [lums_mast[iptaper], lums_mast[iTAPER[ikeep]] ]
                 lums_mast = lout
              endif
           endif

           lums_mast = alog10(lums_mast)
           larray = dblarr(nsn,NTRIALS)
           for t=0, NTRIALS-1 do begin
              it = round(randomu(seed,nsn)*NTRIALS)
              lums_t = lums_mast[it]
              lums_t = lums_t[sort(lums_t)]
              dL = interpol(obs_str.dL_us, findgen(nsn_real)/nsn_real, $
                            findgen(nsn)/(1.*nsn) )
              lums_t = lums_t + randomn(seed,nsn)*dL
              ;lums_t += randomn(seed,nsn)*obs_str.dL_us ; add errors
              lums_t = lums_t[sort(lums_t)]
              larray[*,t] = lums_t


           endfor

           ;; Calculate probs
           lprob = dblarr(nsn)
           chisq = dblarr(nsn)
           parr  = dblarr(nl,nsn)
           for i=0, nsn-1 do begin
              cur_arr = larray[i,*]
              lyhist  = opthist(cur_arr,lhmin,lhmax-lhstep,lhstep)
              lyhist /= double(NTRIALS)
              parr[*,i] = lyhist
           endfor
           results[xx,yy,pp] = ptr_new(parr)

           counter++
           if 1.*counter/(1.*n_array_elems) gt 0.1 then begin
              if calctime then begin
                 telapsed = systime(/sec) - time0
                 tremain  = telapsed*9.
                 message, 'Time now: ' + strn(systime()), /info
                 message, 'Time remaining/min: '+strn(tremain/60.), /info
                 calctime = 0
              endif

              pc_done += 10
              print, strn(pc_done)+'%... ', f='(a7,$)'
              counter = 0
              wait, 0.05

           endif
        endfor
     endfor

  endfor
  SKIPRESULTS:
  print, ' '

  LPMOD = {Llo:llo, Lhi:lhi, GAMMA_L:pl, results:results}
  save, LPMOD, file=OUTNAME+sav_extension+'.sav'
end

;=============================================================

pro compare_obs_mod, LPOBS=LPOBS, LPMOD=LPMOD, $
                     LSTEP=LSTEP, LMIN=LMIN, LMAX=LMAX, NSN_OBS=NSN_OBS, $
                     SILENT=SILENT, BESTFITPARS=BESTFITPARS, $
                     PDNAME=PDNAME

  LLo  = LPMOD.Llo
  LHi  = LPMOD.LHi
  GAML = LPMOD.GAMMA_L
  mod_probs = LPMOD.results
  lhist = outxhist(lmin,lmax,lstep)
  nl    = n_elements(lhist)
  nhist = findgen(NSN_OBS)
  dum = ''
  nsn = nsn_obs

  nlo  = n_elements(Llo)
  nhi  = n_elements(Lhi)
  ngam = n_elements(GAML)
  loglikarr = dblarr(nlo,nhi,ngam)
  probarr   = dblarr(nlo,nhi,ngam)
  chisqarr  = dblarr(nlo,nhi,ngam)
  prob_i_best = -99.
  for xx=0, nlo-1 do begin
     for yy=0, nhi-1 do begin
        for gg=0, ngam-1 do begin
           parr_mod = *(mod_probs[xx,yy,gg])
           loglik_i = total(alog(parr_mod)+alog(LPOBS),/nan)
           loglikarr[xx,yy,gg] = loglik_i

           if total(size(lpobs,/dim) eq size(parr_mod,/dim)) ne 2 then begin
              message, 'Error - dim(LPOBS) .ne. dim(PARR_MOD)'
           endif

           prob_prof = (total(lpobs*parr_mod, nan=1, 1));[0:-1]
           ifin = where(finite(prob_prof))
           prob_i2 = total(alog(prob_prof[ifin]))

           prob_i = prob_i2
           probarr[xx,yy,gg] = prob_i

           if prob_i gt prob_i_best then begin
              prob_i_best = prob_i
              if ~keyword_set(SILENT) then print, llo[xx], lhi[yy], gamL[gg], -2.*prob_i
              llo_best = llo[xx]
              lhi_best = lhi[yy]
              gam_best = gamL[gg]
              void = max(parr_mod,dim=1,imax)
              lfit = lhist[(array_indices(parr_mod,imax))[0,*]]
           endif
        endfor
     endfor
  endfor

  ;; dimension reduction, if needed.
  chisqarr = probarr * (-2.)
  probarr  = exp(probarr-max(probarr))
  probarr /= (total(probarr,/nan)*1.2)

  ;; best-fit values
  minchisq = min(chisqarr, imax,/nan)
  xyg  = array_indices(chisqarr,imax)
  llo_best = xyg[0]
  lhi_best = xyg[1]
  gam_best = xyg[2]

  prob2d = max(probarr,dim=3,/nan)
  prob2d_hig = max(probarr,dim=1,/nan)
  prob2d_log = max(probarr,dim=2,/nan)

  chisq2d = min(chisqarr,dim=3,/nan)
  chisq2d_hig = min(chisqarr,dim=1,/nan)
  chisq2d_log = min(chisqarr,dim=2,/nan)

  ;; probability values
  llcum = probarr[sort(probarr)]/max(probarr,/nan)
  pax  = findgen(n_elements(probarr))
  isig = round(interpol(pax, llcum, [0.003,0.05,0.32]))
  lsig = probarr[(sort(probarr))[isig]]

  lsig = minchisq + [3.53, 7.81, 14.16]

  ;; best fit vals
  i95 = where(chisqarr lt lsig[0], n95)
  xyg_95 = array_indices(chisqarr,i95)
  llo_lo = min(llo[xyg_95[0,*]])
  llo_hi = max(llo[xyg_95[0,*]])

  llo_m = llo # replicate(1.,nhi)
  lhi_m = replicate(1.,nlo) # lhi
  weights = exp(-0.5*chisq2d)
  i1sig = where(chisq2d lt lsig[0])
  llo_best_w = wmean(llo_m[i1sig], weights[i1sig])
  lhi_best_w = wmean(lhi_m[i1sig], weights[i1sig])

  lhi_lo = min(lhi[xyg_95[1,*]])
  lhi_hi = max(lhi[xyg_95[1,*]])

  gam_lo = min(gaml[xyg_95[2,*]])
  gam_hi = max(gaml[xyg_95[2,*]])

  llo_bf = [llo_best_w, llo_hi-llo_best_w, llo_best_w-llo_lo]
  lhi_bf = [lhi_best_w, lhi_hi-lhi_best_w, lhi_best_w-lhi_lo]
  gam_bf = [gaml[gam_best], gam_hi-gaml[gam_best], gaml[gam_best]-gam_lo]

  BESTFITPARS = {llo:llo_bf, lhi:lhi_bf, gam:gam_bf}

  openw, 1, 'LFUNC_results.tex'
  vars = ['$\log(L_{\rm lo}/L_\odot)$','$\log(L_{\rm hi}/L_\odot)$','$\Gamma_{\rm L}$']
  for i=0, 2 do begin
     case i of
        0: var = llo_bf
        1: var = lhi_bf
        2: var = gam_bf
     endcase
     printf, 1, vars[i], var, form='(a26, " & = & $", f5.2, "^{+",f4.2,"}_{-",f4.2,"}$ \\" )'
  endfor
  close, 1
  print, ' '


  ;; Save distribution to output files
  im_c  = (-1.*chisq2d)+max(chisq2d,/nan)
  lev_c = reverse((-1.*lsig)+max(chisq2d,/nan))
  save, im_c, llo, lhi, lev_c, file='LUM_CONTOUR.sav'


  ;; best fitting profile ---------------------
  maxp = min(chisqarr,imax)
  xyg_best = array_indices(chisqarr,imax)
  parr_mod = *(mod_probs[xyg_best[0],xyg_best[1],xyg_best[2]])

  POSIM = [0.12,0.15,0.97,0.97]

  parr_n = parr_mod
  parr_sig = parr_mod
  for i=0, nsn-1 do begin
     parr_n[*,i] /= max(parr_n[*,i])
     pm_sort = (parr_mod[*,i])[sort(parr_mod[*,i])]
     ind = findgen(n_elements(pm_sort))
     isig = interpol(ind,total(pm_sort,/cum),[0.001,0.003,0.05,0.68])
     lsig = parr_mod[(sort(parr_mod[*,i]))[isig],0 ]
  endfor

  prange = [0,1]
  bytim = round(1.*bytscl(parr_n, min=prange[0], max=prange[1], top=ncols)/32.)*32.

  for i=0, nsn-1 do begin
     prof = total(lpobs[*,i],/cum)
     prof /= max(prof)
     lmed = interpol(lhist,prof,0.50)
     elo = lmed - interpol(lhist,prof,0.16)
     ehi = interpol(lhist,prof,0.16) - lmed
     icur = (i+1.)/float(nsn)
  endfor

  ;; If Gamma_L forced to be = -1.675 ----------------------------
  gamL_pre = -1.675
  void = min(abs(gamL-gamL_pre),ig)
  prob2D_G = probarr[*,*,ig]
  prob2D_G /= total(prob2D_G,/nan)
  chisq2d_G = chisqarr[*,*,ig]
  max = max(prob2D_G,imax,/nan)
  xybest = array_indices(prob2D_G,imax)
  llo_best_g = llo[xybest[0]]
  lhi_best_g = lhi[xybest[1]]

  lsig = min(chisq2d_G,/nan) + [2.3,6.18,11.83]

  POSIM = [0.15,0.15,0.97,0.97]
  pwid = 0.325
  im_c  = (-1.*chisq2d_G)+max(chisq2d_G,/nan)
  lev_c = reverse((-1.*lsig)+max(chisq2d_G,/nan))

  ;; trim white space
  i999 = where(chisq2d_G lt (min(chisq2d_G,/nan) + 20.))
  xy999 = array_indices(chisq2d_G, i999)
  llo_minmax = [min(llo[xy999[0,*]]), max(llo)]

  ;; M-L relations ------------------------------
  MLnames = ['Geneva-nr','Geneva-r','KEPLER','STARS','MIST']
  AL = [2.67, 2.92, 2.68, 2.67, 2.84]
  GL = [1.93, 1.82, 1.93, 2.02, 1.83]


  ;; COnvert to Mass plane ---------------------------
  mlo_ax = 10^((llo-AL[3])/GL[3]) & mlo_ax_range = [min(mlo_ax),max(mlo_ax)]
  mhi_ax = 10^((lhi-AL[3])/GL[3]) & mhi_ax_range = [min(mhi_ax),max(mhi_ax)]

  mlo_ax2 = findgen(nlo)/(nlo-1.)*(mlo_ax_range[1]-mlo_ax_range[0]) + mlo_ax_range[0]
  mhi_ax2 = findgen(nhi)/(nhi-1.)*(mhi_ax_range[1]-mhi_ax_range[0]) + mhi_ax_range[0]

  mlo_ax_ind = interpol(findgen(nlo),mlo_ax,mlo_ax2)
  mhi_ax_ind = interpol(findgen(nhi),mhi_ax,mhi_ax2)

  im_mass = interpolate(im_c, mlo_ax_ind, mhi_ax_ind, /grid, cub=-0.5)

  mass_chisq = interpolate(chisq2d_G, mlo_ax_ind, mhi_ax_ind, /grid, cub=-0.5)
  koch_xy = [interpol(findgen(nlo), mlo_ax2, 6.3), $
             interpol(findgen(nhi), mhi_ax2, 15.8) ]
  koch_chisq = interpolate(chisq2d_G, koch_xy[0], koch_xy[1])

  ;; Masses ------------------------------
  MLnames = ['Geneva-nr','Geneva-r','KEPLER','STARS','MIST']
  AL = [2.67, 2.92, 2.68, 2.67, 2.84]
  GL = [1.93, 1.82, 1.93, 2.02, 1.83]

  openu, 1, 'LFUNC_results.tex', /append
  printf, 1, ' '
  for i=0, n_elements(MLnames)-1 do begin
     mlo_best = 10^(( abs(llo_bf+[0,llo_bf[0],-llo_bf[0]])-AL[i])/GL[i])
     mlo_best = abs(mlo_best-[0,mlo_best[0],mlo_best[0]])

     mhi_best = 10^(( abs(lhi_bf+[0,lhi_bf[0],-lhi_bf[0]])-AL[i])/GL[i])
     mhi_best = abs(mhi_best-[0,mhi_best[0],mhi_best[0]])

     printf, 1, MLnames[i], mlo_best, mhi_best, $
             form='(a12," & ",f4.1,"$^{+",f3.1,"}_{-",f3.1,"}$ & ", '+$
             'f4.1,"$^{+",f3.1,"}_{-",f3.1,"}$ \\")'
  endfor
  close, 1
  if ~keyword_set(silent) then spawn, 'cat LFUNC_results.tex'

end

;=============================================================
;=============================================================
;=============================================================

pro mainproc

  NTRIALS_OBS = 100000
  NTRIALS_MOD = 100000
  LSTEP   = 0.01
  LMIN    = 4.0
  LMAX    = 5.8

  NSN_OBS=25

  NLLO    = 30
  NLHI    = 29
  NGAM    = 23

  LLO_RANGE=[4.0,4.6]
  LHI_RANGE=[5.0,5.7]
  GAM_RANGE=[0.5,-2.2]

  message, 'Creating Lobs...', /info
  create_Lobs, NTRIALS=NTRIALS_OBS, LSTEP=LSTEP, LMIN=LMIN, LMAX=LMAX, $
               LPOBS=LPOBS, NSN_OBS=NSN_OBS, OBS_STR=OBS_STR
  restore, 'LPOBS.sav'

  message, 'Creating Lmod...', /info
  NO_RESULTS = 1
  TAPER_LLO = 0. ;0.1
  TAPER_LHI = 0. ;0.1


  message, 'Doing synth tests...', /info
  LLO_IN = 4.45
  LHI_IN = 5.28
  GAM_IN = -1.675
  NSN_SYNTH = 31
  create_synth_obs, NTRIALS=NTRIALS_OBS, LSTEP=LSTEP, LMIN=LMIN, LMAX=LMAX, $
                    LPOBS_SYNTH=LPOBS_SYNTH, NSN_OBS=NSN_SYNTH, OBS_STR=OBS_STR, $
                    LLO_IN=LLO_IN, LHI_IN=LHI_IN, GAM_IN=GAM_IN

  LPNAME = 'LPMOD_SYNTH_N'+strn(NSN_SYNTH)+$
           '_LL'+strtrim(string(LLO_IN,'(f4.2)'),2) + $
           '_LH'+strtrim(string(LHI_IN,'(f4.2)'),2) + $
           '_G'+strtrim(string(GAM_IN,'(f6.3)'),2)
  message, 'filename = '+LPNAME, /info

  if file_test(LPNAME+'.sav') then begin
     message, 'Restoring '+LPNAME, /info
     restore, LPNAME+'.sav'
  endif else begin
     message, 'Creating '+LPNAME, /info
     create_modelgrid,  LSTEP=LSTEP, LMIN=LMIN, LMAX=LMAX, $
                        NLLO=NLLO, NLHI=NLHI, NGAM=NGAM, $
                        NSN_OBS=NSN_SYNTH, LPMOD=LPMOD, NTRIALS=NTRIALS_MOD, $
                        OBS_STR=OBS_STR, NO_RESULTS=0, $
                        LLO_RANGE=LLO_RANGE, LHI_RANGE=LHI_RANGE, GAM_RANGE=GAM_RANGE, $
                        OUTNAME=LPNAME
  endelse

  PDNAME = LPNAME + '_gfix-plotdata.sav'
  compare_obs_mod, LPOBS=LPOBS_SYNTH, LPMOD=LPMOD, $
                   LSTEP=LSTEP, LMIN=LMIN, LMAX=LMAX, NSN_OBS=NSN_SYNTH, $
                   SILENT=SILENT, $
                   PDNAME=PDNAME
  
end


