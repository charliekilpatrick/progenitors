function getcols, dev, struct=struct


;if n_params() eq 0 then begin
;   print, " "
;   print, " col = getcols(dev, struct=struct)"
;   return, 0
;endif

if n_params() eq 0 then dev = !d.name


if dev eq 'ps' then begin
;.......... Bl, Rd, Grn, Blu, Cyn, Yel, Org, Gry, wht, Mg, pnk, lgr, lbl
    TVLCT, [0, 255,   0,   0,   0, 245, 255, 120, 255, 255, 255, 125, 125, 210], $
           [0,   0, 225,  20, 205, 245, 175, 120, 255,   0, 125, 255, 125, 210], $
           [0,   0,   0, 255, 205,   0,   0, 120, 255, 255, 125, 125, 255, 210]
    col_ind =  [0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,13]
    col_str = { white:0, $
             black:8, $
             red:1, $
             green:2, $ 
             blue:3, $
             cyan:4, $ 
             yellow:5, $  
             orange:6, $  
             magenta:9, $   
             grey:7, $ 
             pink:10, $ 
             lgreen:11, $ 
             dgrey:13, $ 
             lblue:12 $
           }
endif else begin
    cube = indgen(256,256,256, /long)
    col_ind  = [ cube[255,255,255], $ ; white
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
    col_str = { white:cube[255,255,255], $
             red:cube[255,0,0], $
             green:cube[0,255,0], $ 
             blue:cube[0,0,255], $
             cyan:cube[0,255,255], $ 
             yellow:cube[255,255,0], $  
             orange:cube[255,125,0], $  
             magenta:cube[255,0,255], $   
             grey:cube[120,120,120], $ 
             pink:cube[255,125,125], $ 
             lgreen:cube[125,255,125], $ 
             black:cube[0,0,0], $ 
             dgrey:cube[50,50,50], $ 
             lblue:cube[125,125,255] $
           }
endelse

if keyword_set(struct) then col = col_str else col = col_ind

return, col

end

