pro setplot, col, dev=dev, font=font, thick=thick, $
             file=file, xs=xs, ys=ys

if n_params() eq 0 then begin
   print, " "
   print, " SETPLOT, col, dev=dev, font=font"
   return
endif

if not keyword_set(dev) then dev = 'x'
if not keyword_set(font) then font = -2*(dev eq 'x') + 1
if not keyword_set(thick) then thick = (dev eq 'ps')*4+1

set_plot, dev
cleanplot, /sil
!x.style = 1
!y.style = 1
!p.font = font
!p.multi = [0, 1, 1]
if dev eq 'ps' then begin
    charsize = 1.8
    TVLCT, [0, 255,   0,   0,   0, 255, 170, 255, 255], $
           [0,   0, 155,   0, 205, 152, 170, 255,   0], $
           [0,   0,  50, 255, 205,   0, 170, 255, 255]
    col = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    ss  = 2
    device, /encapsulated, xs=25, ys=18, file='psfile.eps',$
            /color
endif else begin
    if not keyword_set(thick) then thick = 1.
    charsize = 1.8
    device, decomp=1
    cube = indgen(256,256,256, /long)
    col  = [ cube[255,255,255], $ ; white
             cube[255,0,0], $     ; red
             cube[0,255,0], $     ; green
             cube[0,0,255], $     ; blue
             cube[0,255,255], $   ; cyan
             cube[255,255,0], $   ; yellow
             cube[255,125,0], $   ; orange
             cube[255,0,255], $   ; magenta
             cube[120,120,120], $ ; grey
             cube[255,125,125], $  ; pink 
             cube[125,255,125], $  ; light green 
             cube[125,125,255] $  ; light blue 
           ]
    ss = 2
endelse
circ = dindgen(34)*!pi*2/32
usersym,sin(circ),cos(circ), thick=thick
!p.thick = 1.*thick
!x.thick = 1.*thick
!y.thick = 1.*thick
!p.charthick = 1.*thick/2.
!p.charsize  = charsize
!p.symsize = ss

if dev eq 'ps' then device, /encap, /col, xs=xs, ys=ys, $
                            file=file

end

