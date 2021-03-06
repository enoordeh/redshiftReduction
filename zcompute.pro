;+
; NAME:
;   zcompute
;
; PURPOSE:
;   Compute relative redshift of object(s) vs. eigen-templates.
;
; CALLING SEQUENCE:
;   zans = zcompute(objflux, objivar, starflux, [starmask, nfind=, $
;    poffset=, pspace=, pmin=, pmax=, mindof=, width=, minsep=, $
;    plottitle=, /doplot, /debug ]
;
; INPUTS:
;   objflux    - Object fluxes [NPIXOBJ,NOBJ]
;   objivar    - Object inverse variances [NPIXOBJ,NOBJ]
;   starflux   - Eigen-template fluxes [NPIXSTAR,NTEMPLATE]
;
; OPTIONAL INPUTS:
;   starmask   - Eigen-template mask; 0=bad, 1=good [NPIXSTAR]
;   nfind      - Number of solutions to find per object; default to 1.
;   poffset    - Offset between all objects and templates, in pixels.
;                A value of 10 indicates that STARFLUX begins ten pixels
;                after OBJFLUX, i.e. OBJFLUX[i+10] = STARFLUX[i] for the
;                case when the relative redshift should be zero.  If the
;                wavelength coverage of the templates is larger, then the
;                value of ZOFFSET will always be negative.
;                [Scalar or vector of size NOBJ]
;   pspace     - The spacing in redshifts to consider; default to 1 [pixels];
;                [Scalar or vector of size NOBJ]
;   pmin       - The smallest redshift to consider [pixels].
;                [Scalar or vector of size NOBJ]
;   pmax       - The largest redshift to consider [pixels].
;                [Scalar or vector of size NOBJ]
;   mindof     - Minimum number of degrees of freedom for a valid fit;
;                default to 10.
;   width      - Parameter for FIND_NMINIMA(); default to 3 * PSPACE.
;   minsep     - Parameter for FIND_NMINIMA(); default to the same as WIDTH.
;   plottitle  - ???
;   doplot     - ???
;   debug      - ???
;
; OUTPUTS:
;   zans       - Output structure [NOBJECT,NFIND] with the following elements:
;                z : The relative redshift.
;                z_err : Error in the redshift, based upon the local quadratic
;                        fit to the chi^2 minimum. 
;                chi2 : Fit value for the best (minimum) chi^2
;                dof : Number of degrees of freedom, equal to the number of
;                      pixels in common between the object and templates
;                      minus the number of templates.
;                theta : Mixing angles [NTEMPLATE].  These are computed at the
;                        nearest integral redshift, e.g. at ROUND(ZOFFSET).
;
; COMMENTS:
;   Fits are done to chi^2/DOF, not to chi^2.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/TEMPLATEFILES
;
; PROCEDURES CALLED:
;   find_nminima()
;   splog
;
; INTERNAL SUPPORT ROUTINES:
;   create_zans()
;
; REVISION HISTORY:
;   10-Jul-2000  Written by D. Schlegel, Princeton
;   23-Sep-2002  adapted by MD for use by DEIMOS, DEEP2 survey
;------------------------------------------------------------------------------
; Create output structure
function create_zans, nstar, nfind

   zans1 = create_struct( $
    name = 'ZANS'+strtrim(string(nstar),1), $
    'z'     , 0.0, $
    'z_err' , 0.0, $
    'chi2'  , 0.0, $
    'dof'   ,  0L, $
    'theta' , fltarr(nstar) )

   return, replicate(zans1, nfind)
end

;------------------------------------------------------------------------------
function zcompute, objflux, objivar, starflux, starmask, nfind=nfind, $
                   poffset=poffset, pspace=pspace, pmin=pmin, pmax=pmax, $
                   mindof=mindof, width=width, minsep=minsep, $
                   plottitle=plottitle, doplot=doplot1, debug=debug, fname=fname, tclass=tclass, physcheck=physcheck

   if (NOT keyword_set(nfind)) then nfind = 1
   if (NOT keyword_set(pspace)) then pspace = 1 ; *Emil* pspace passed through _EXTRA
   if (NOT keyword_set(width)) then width = 3 * pspace ; 
   if (NOT keyword_set(plottitle)) then plottitle = ''

;   on_error, 0 ;stop at error
   ; print, "debug = ", debug
   ; Plot if either /DOPLOT or /DEBUG is set.
   if (keyword_set(doplot1)) then doplot = doplot1
   if (keyword_set(debug)) then doplot = 1
   ; print, "doplot = ", doplot

   ;---------------------------------------------------------------------------
   ; Check dimensions of object vectors

   ndim = size(objflux, /n_dimen)
   dims = size(objflux, /dimens)
   npixobj = dims[0]
   if (ndim EQ 1) then nobj = 1 $
    else if (ndim EQ 2) then nobj = dims[1] $
    else message, 'OBJFLUX is neither 1-D or 2-D'

   if total(abs(size(objflux, /dimens)-size(objivar, /dimens))) NE 0 $
    OR size(objflux, /n_dimen) NE size(objivar, /n_dimen) THEN  $
    message, 'Dimensions of OBJFLUX and OBJIVAR do not match'

   ;---------------------------------------------------------------------------
   ; If multiple object vectors exist, then call this routine recursively.

   if (nobj GT 1) then begin
      t0 = systime(1)
      for iobj=0L, nobj-1 do begin
         if (n_elements(poffset) EQ 1) then poffset1 = poffset $
          else poffset1 = poffset[iobj]
         if (n_elements(pspace) EQ 1) then pspace1 = pspace $
          else pspace1 = pspace[iobj]
         if (n_elements(pmin) EQ 1) then pmin1 = pmin $
          else pmin1 = pmin[iobj]
         if (n_elements(pmax) EQ 1) then pmax1 = pmax $
          else pmax1 = pmax[iobj]

         zans1 = zcompute(objflux[*,iobj], objivar[*,iobj], $
          starflux, starmask, nfind=nfind, $
          poffset=poffset1, pspace=pspace1, pmin=pmin1, pmax=pmax1, $
          mindof=mindof, width=width, minsep=minsep, $
          plottitle=plottitle+ ' -- Object #'+strtrim(string(iobj+1),2), $
          doplot=doplot, debug=debug)
         if (iobj EQ 0) then zans = zans1 $
          else zans = [[zans], [zans1]]
         splog, 'Object #', iobj, '  Elap time=', systime(1)-t0, $
          ' (sec)  z=', zans1[0].z, ' (pix)'
      endfor
      return, zans
   endif

;---------------------------------------------------------------------------

   if (NOT keyword_set(mindof)) then mindof = 10
   if (NOT keyword_set(width)) then width = 3 * pspace
   if (NOT keyword_set(minsep)) then minsep = width

   ndim = size(starflux, /n_dimen)
   dims = size(starflux, /dimens)
   npixstar = dims[0]
   if (ndim EQ 1) then nstar = 1 $
    else if (ndim EQ 2) then nstar = dims[1] $
    else message, 'STARFLUX is neither 1-D or 2-D'

   if (NOT keyword_set(starmask)) then begin
      starmask = bytarr(npixstar) + 1
   endif else begin
      if (n_elements(starmask) NE npixstar) then $
       message, 'Dimensions of STARFLUX and STARMASK do not match'
   endelse

   if (n_elements(poffset) EQ 0) then poffset = 0
   if (n_elements(poffset) GT 1) then $
    message, 'ZOFFSET must be a scalar'
   pixoffset = round(poffset)

   if (n_elements(pmin) EQ 0) then $
    pmin = -2 * ((npixobj < npixstar) + 1) + pixoffset
   if (n_elements(pmax) EQ 0) then $
    pmax = pmin + 2 * ((npixobj < npixstar) - 1)

   if (n_elements(pmin) GT 1) then $
    message, 'PMIN must be a scalar'
   if (n_elements(pmax) GT 1) then $
    message, 'PMAX must be a scalar'
   if (pmin GT pmax) then $
    message, 'PMAX must be greater than or equal to PMIN'

    ; print, "pmax", pmax
    ; print, "pmin", pmin
    ; print, "pspace", pspace
   nlag = ((pmax - pmin + 1) / pspace) > 1
   lags = - lindgen(nlag) * pspace + pixoffset - long(pmin) ; must be integers

   chi2arr = fltarr(nlag)
   debugarr = fltarr(7,nlag)
   dofarr = fltarr(nlag)
   ccorarr = fltarr(nlag)
   normarr = fltarr(nlag)
   thetaarr = fltarr(nstar,nlag)
   physicalarr = intarr(nlag)
   zans = create_zans(nstar, nfind)

   ;---------------------------------------------------------------------------

   sqivar = sqrt(objivar > 0) ;force this to be positive
   objmask = objivar NE 0

   for ilag=0L, nlag-1 do begin
      j1 = lags[ilag] < (npixstar-1L)
;;; FOR THE NEGATIVE LAGS, JUST SHIFT THE OBJECT SPECTRUM.
      if (j1 LT 0) then i1 = -j1 $
;;; OTHERWISE SHIFT THE TEMPLATE SPECTRUM (SET i1=0). 
       else i1 = 0L
;;; IF j1 IS NEGATIVE, THEN SET j1=0 (SHIFT OBJECT SPECTRUM
;;; INSTEAD).
      j1 = j1 > 0L
;;; SET j2 SO THAT WE DON'T EXCEED THE LENGTH OF THE TEMPLATE.
      j2 = npixstar-1 < (npixobj+j1-i1-1L)
;;; SIMILARLY 
      i2 = i1 + j2 - j1


IF abs(j2-j1) LT 100 THEN BEGIN
    chi2arr[ilag] = 1E50
    dofarr[ilag] = 1
    thetaarr[*,ilag] = 1
    debugarr[ilag] = 1
ENDIF ELSE BEGIN
      ; print, "objflux size: ", n_elements(objflux[i1:i2])
      ; print, "starflux size", n_elements(starmask[j1:j2])
      chi2arr[ilag] = computechi2( objflux[i1:i2], $
                       sqivar[i1:i2] * starmask[j1:j2], starflux[j1:j2,*], $
                            acoeff=acoeff, dof=dof)
      dofarr[ilag] = dof
      thetaarr[*,ilag] = acoeff
      debugarr[*,ilag] = [ilag, lags[ilag], n_elements(objflux[i1:i2]), n_elements(starmask[j1:j2]), mean(objflux[i1:i2]), mean(starflux[j1:j2,*]), chi2arr[ilag]/dofarr[ilag]]
;


; *** "physical solution" check below commented out by Emil

; is this a physical solution or not??- either both positive, or one +
; ;                                       and the other modest in strength

; ** Emil adjustment to physarr check:
; If star, don't apply the check
; If CV, apply check only up to acoeff[2] (CV only has 3 eigenspectra)
; If gal/AGN, apply up to acoeff[3]

; print, "n_elements(acoeff) = ", n_elements(acoeff)
; print, "acoeff = ", acoeff

; ** EMIL's PHYSARR Check

if physcheck eq 1 then begin
  if tclass eq 'STAR' then begin ;
    ; print, 'TCLASS = STAR' 
    if (acoeff[0] ge 0.) then physicalarr[ilag] = 1 $
      else physicalarr[ilag] = 0
  endif

  if tclass eq 'CV' then begin
    ; print, 'TCLASS = CV'
    coeff_adjust = [acoeff[0],acoeff[1],acoeff[2]] 
    coeffcopy = coeff_adjust(sort(coeff_adjust))
    if (coeffcopy[0] GT 0.) OR $ ; i.e. all a_coeff are positive
      ((coeffcopy[2] GT 0) AND (abs(coeffcopy[0]) lt 0.5*coeffcopy[2])) $ ; i.e. largest coeff is atleast 2x smallest coefficient
    then begin 
      physicalarr[ilag] = 1 
    endif else begin 
      if coeffcopy[2] GT 0 then physicalarr[ilag] = 0 $
        else physicalarr[ilag] = -1 
    endelse
  endif

  if (tclass eq 'AGN') OR (tclass eq 'GALAXY') then begin
    ; print, 'TCLASS = AGN or GAL'
    coeff_adjust = [acoeff[0],acoeff[1],acoeff[2],acoeff[3]] ; **** EMIL ADDED acoeff[3] &  all [3] below used to be [2]
    coeffcopy = coeff_adjust(sort(coeff_adjust))
    if (coeffcopy[0] GT 0.) OR $ 
      ((coeffcopy[3] GT 0) AND (abs(coeffcopy[0]) lt 0.5*coeffcopy[3])) $
    then begin 
      physicalarr[ilag] = 1 
    endif else begin 
      if coeffcopy[3] GT 0 then physicalarr[ilag] = 0 $
        else physicalarr[ilag] = -1 
    endelse
  endif
endif

if physcheck eq 2 then begin
  if tclass eq 'STAR' then begin ;
    ; print, 'TCLASS = STAR' 
    if (acoeff[0] ge 0.) then physicalarr[ilag] = 1 $
      else physicalarr[ilag] = 0
  endif

  if tclass eq 'CV' then begin
    ; print, 'TCLASS = CV'
    coeff_adjust = [acoeff[0],acoeff[1],acoeff[2]] 
    coeffcopy = coeff_adjust(sort(abs(coeff_adjust)))
    if coeffcopy[2] GT 0 then physicalarr[ilag] = 1 $
      else physicalarr[ilag] = 0
  endif

  if (tclass eq 'AGN') OR (tclass eq 'GALAXY') then begin
    ; print, 'TCLASS = AGN or GAL'
    coeff_adjust = [acoeff[0],acoeff[1],acoeff[2],acoeff[3]] 
    coeffcopy = coeff_adjust(sort(abs(coeff_adjust)))
    if coeffcopy[3] GT 0 then physicalarr[ilag] = 1 $
      else physicalarr[ilag] = 0
  endif
endif

if physcheck eq 3 then begin
  if (acoeff[0] ge 0.) then begin
    physicalarr[ilag] = 1 
    
  endif else begin
    physicalarr[ilag] = 0
    ; print, '********PHYSCHECK=3 ACOEFF0 < 0**************'
  endelse
endif

; OLD PHYSARR Check

;       if n_elements(acoeff) eq 1 then begin ; never the case (even for stars) when we have npoly
;           if (acoeff[0] ge 0.) then physicalarr[ilag] = 1 $
;           else physicalarr[ilag] = 0
;       endif else begin
; ;          if ((acoeff[0] GE 0. AND acoeff[1] GE 0.) OR $ ;both positive
; ;              (acoeff[0] GT 0. AND abs(acoeff[1]) lt .5*acoeff[0])  OR $ ;1 +
; ;              (acoeff[1] GT 0. AND abs(acoeff[0]) lt .5*acoeff[1])) then $
; ;            physicalarr[ilag] = 1 ELSE physicalarr[ilag] = 0 
; ;          coeff_adjust = [acoeff[0]/3.,acoeff[1]*14./3,acoeff[2]]
;           coeff_adjust = [acoeff[0],acoeff[1],acoeff[2],acoeff[3]] ; **** EMIL ADDED acoeff[3] &  all [3] below used to be [2]

;           coeffcopy = coeff_adjust(sort(coeff_adjust))
; 	  if ( coeffcopy[0] GT 0. OR $ 
; 	       (coeffcopy[3] GT 0 AND abs(coeffcopy[0]) lt 0.5*coeffcopy[3]) ) $
; 	      then begin 
; 	         physicalarr[ilag] = 1 
;           endif else begin 
; 	      if coeffcopy[3] GT 0 then physicalarr[ilag] = 0 $
; 		                   else physicalarr[ilag] = -1 
;           endelse
;       endelse





; EXTRA JUNK, used in testing
;      ccorarr[ilag] = total(objflux[i1:i2]*(objivar[i1:i2] > 0.)* $
;           starmask[j1:j2]*starflux[j1:j2] )/(i2-i1)
;      normarr[ilag] = sqrt(variance(objflux[i1:i2]*(objivar[i1:i2] > 0.))* $
;            variance(starmask[j1:j2]*starflux[j1:j2]) )
ENDELSE


      ; Burles counter of lag number...
      if(ilag mod 100 eq 50) then $ 
        print, format='("Lag ",i5," of ",i5,a1,$)', $
         ilag, nlag, string(13b)

   endfor
; print, "chi2arr:"
; print, chi2arr

if (keyword_set(doplot)) then begin
print, "************ PLOTTING PRE-SMOOTHING CHI2ARR/DOF vs LAG ************"
cgplot, lags, chi2arr/dofarr, color=FSC_color('green'), XTicklen=1.0, YTicklen=1.0, AxisColor='white',$
    YTitle = textoidl('\chi^2/DOF'), XTitle = 'Lag', xrange=[5050,9450]
pause

; print, "***********MINCHI2/DOF***********:", min(chi2arr/dofarr)
; print, debugarr
  ; print, "nlag"
  ; print, nlag
endif

   if nlag gt 1000 then begin ;remove long range trends
      medchi2 = median(chi2arr)
      diffchi2arr =  djs_median(chi2arr,  width= 200*5/pspace,  boundary='reflect') ; *** Emil changed width from 200 to 200*5/pspace so the width in redshift is independent of pspace
      chi2arr = chi2arr-diffchi2arr + medchi2 ;redefine chisqr array
      print, "chi2arr smoothed - long trends removed"
   endif

if (keyword_set(doplot)) then begin
; print, "***********MINCHI2/DOF (POST SMOOTHING)***********:", min(chi2arr/dofarr)
   ;-----
   ; Fit this chi2/DOF minimum

; WRITE_CSV, fname, debugarr

    lengthtest=n_elements(chi2arr)
    xaxistest=indgen(lengthtest)
    ; print, "xaxistest", xaxistest
    ; print, "chi2arr", chi2arr
    print, "min chi2arr/dof: ", min(chi2arr/dof)
    print, "index of minchi2arr/dof: ", where(chi2arr/dof EQ min(chi2arr/dof))
    print, "************ PLOTTING POST-SMOOTHING CHI2ARR/DOF vs LAG  ************"
     cgoplot, lags, chi2arr/dofarr, color=FSC_color('red')
    print, "post smoothing min chi2arr/dof: ", min(chi2arr/dof)
     pause
    ; oplot, xaxistest, dofarr, color=FSC_color('green')
    ; wait, 2
endif
  ;   indx = where(dofarr GT mindof, ct)
  ; wait, 1
  ;   oplot, xaxistest[indx], chi2arr[indx], color=FSC_color('red')

   indx = where(dofarr GT mindof, ct)
   ; print, "ct: ", ct
   ; print, "mindof: ", mindof
   ; print, "indx: ", indx
   ; print, "width: ", width
   ; print, "ct = ", ct
   ; print, 'width = ', width
   if (ct GT width) then begin
;      xpeak1 = find_npeaks(-chi2arr[indx]/dofarr[indx], lags[indx], $
;       nfind=nfind, minsep=minsep, width=width, $
;       ypeak=ypeak1, xerr=xerr1, npeak=npeak)
;      zans[0:npeak-1].chi2 = -ypeak1
      ; print, "lags[indx] ", lags[indx]
      ; print, [Transpose(lags[indx]),Transpose(chi2arr[indx]/dofarr[indx])]
      ; print, chi2arr[indx]
      xpeak1 = find_nminima(chi2arr[indx], lags[indx], $ 
       dofarr=dofarr[indx], nfind=nfind, minsep=minsep, width=width, $
       ypeak=ypeak1, xerr=xerr1, npeak=npeak, errcode=errcode, $;/flatnoise, $
       physicalarr=physicalarr,  plottitle=plottitle, doplot=doplot)
      ; print, "xpeak1 = ", xpeak1
      ; print, "ypeak1 = ", ypeak1
      zans[0:npeak-1].z = poffset - xpeak1
      ; Set Z_ERR equal to the error-code if it is non-zero
      zans[0:npeak-1].z_err = xerr1 * (errcode EQ 0) + errcode
      ; print, "ypeak1 = ", ypeak1 ; *** Emil This is where negative chi^2 is coming from
      zans[0:npeak-1].chi2 = ypeak1
      for ipeak=0L, npeak-1 do begin
         junk = min(abs(lags-xpeak1[ipeak]), ilag)
         zans[ipeak].dof = dofarr[ilag]
         zans[ipeak].theta = thetaarr[*,ilag]
      endfor
      zans[0:npeak-1].chi2 = zans[0:npeak-1].chi2 * zans[0:npeak-1].dof

      ; Wait for a keystroke...
      if (keyword_set(debug)) then begin
         print, 'Press any key...'
         cc = strupcase(get_kbrd(1))
      endif

   endif else if (ct GE 1) then begin
      zans[0].chi2 = -min(-chi2arr[indx]/dofarr[indx], ilag)
      zans[0].z = poffset - lags[indx[ilag]]
      zans[0].z_err = 0
      zans[0].dof = dofarr[indx[ilag]]
      zans[0].chi2 = zans[0].chi2 * zans[0].dof
      zans[0].theta = thetaarr[*,indx[ilag]]
   endif

   return, zans
end
;------------------------------------------------------------------------------


