pro eps_setup, filename, xs=xs, ys=ys, _extra=_extra_keywords

  if !d.name eq 'PS' then begin

     device, /encap, /col, xs=xs, ys=ys, $
             file=filename, _extra=_extra_keywords

  endif

end
