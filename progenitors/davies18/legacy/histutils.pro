;*****************************************************
;*****************************************************
;*****************************************************

function opthist, var, binmin, binmax, binsize

  outhist = histogram(var, binsize=binsize,min=binmin,max=binmax)

  return, outhist

end

;*****************************************************
;*****************************************************
;*****************************************************

function outxhist, binmin, binmax, binsize

  xarray = findgen((binmax-binmin)/binsize)/((binmax-binmin)/binsize) * $
             (binmax-binmin)+binmin + binsize/2.

  return, xarray

end

