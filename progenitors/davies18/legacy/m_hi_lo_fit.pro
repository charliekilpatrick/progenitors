PRO M_hi_lo_fit, X, A, F, pder

  GAMMA = -1.35
  norm  = 1./(A[1]^GAMMA - A[0]^GAMMA)
  expo  = 1./GAMMA * alog10(X/norm + A[0]^GAMMA)
  F     = 10.^expo

  IF N_PARAMS() GE 4 THEN $
     pder = [ [x], [x] ]
END
