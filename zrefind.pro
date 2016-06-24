;+
; NAME:
;   zrefind
;
; PURPOSE:
;   Re-fit redshifts from a previous call to ZFIND(), but doing a local
;   fit around the previous answers.
;
; CALLING SEQUENCE:
;   result = zrefind( ss1d, objflux, objivar, pwidth=, pspace=, width=, zold= ]
;
; INPUTS:
;   ss1d       - structure of 1d spectra (not used, but passed to zfind)
;   objflux    - Object fluxes [NPIXOBJ,NOBJ]
;   objivar    - Object inverse variances [NPIXOBJ,NOBJ]
;
; REQUIRED KEYWORDS:
;   zold       - Z structure from an initial call to ZFIND().
;
; OPTIONAL KEYWORDS:
;   pwidth     - Search width in pixels about the intial redshift; default to 5
;   pspace     - Keyword for ZCOMPUTE().
;   width      - Keyword for ZCOMPUTE().
;   doplot     - Keyword for ZCOMPUTE
;
; OUTPUTS:
;   result     - Structure with redshift-fit information, modified from ZOLD.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   zfind()
;
; REVISION HISTORY:
;   17-Aug-2000  Written by D. Schlegel, Princeton
;   30-sep-02    hacked by MD for DEIMOS pipeline
;------------------------------------------------------------------------------
function zrefind, ss1d, objflux, objivar, hdr=hdr, loglam=loglam, $
 pwidth=pwidth, pspace=pspace, width=width, zold=zold, doplot=doplot, fname=fname, $
      _EXTRA=EXTRA

  print, "******* ZREFIND *******"
   if (NOT keyword_set(pwidth)) then pwidth = 5

   if n_elements(doplot) eq 0 then doplot=0

   ;----------
   ; Copy results

   result = zold
   nres = n_elements(result)

   if (size(zold,/n_dim) EQ 1) then nfind = 1 $
    else nfind = (size(zold,/dimens))[0]

   ;----------
   ; Identify which redshift measurements are valid and will be re-measured

   qvalid = result.tfile NE '' AND result.tcolumn[0] NE -1 $ ; result.tfile stores an array of the eigenspectra used for top 5 chi^2 results
    AND result.dof NE 0
   ivalid = where(qvalid)
   ; print, "ivalid = ", ivalid
   if (ivalid[0] EQ -1) then return, result

   ;----------
   ; Find all eigen-spectra files being used
   ; The ISELECT array will give indices for unique (and valid) sort strings.

   ntcol = n_elements(result[0].tcolumn)
   tcolstring = strarr(nres)
   for ires=0, nres-1 do $
    tcolstring = string(result[ires].tcolumn, format='('+string(ntcol)+'i)')
   sortstring = result[*].tfile + tcolstring + string(result[*].npoly)
   ; print, "sortstring is:", sortstring

   isort = ivalid[ sort(sortstring[ivalid]) ]
   iselect = isort[ uniq(sortstring[isort]) ]

   for ifile=0, n_elements(iselect)-1 do begin

      i0 = iselect[ifile] ; Index of result specifying this TFILE,TCOLUMN,NPOLY

      ; Identify which redshifts are measured with this template file,
      ; TCOLUMN, and NPOLY.
      indx = where(sortstring EQ sortstring[i0])

      ; Find object numbers corresponding to these indices
      iobj = indx / nfind

      ; Re-fit the redshifts using the specified PWIDTH,PSPACE
      ncol = (reverse(where(result[i0].tcolumn NE -1)))[0] + 1
      
      ; TEST PARAMETERS
      ; zmin = -0.0001
      ; ; zmin = 0.01
      ; zmax = max(ss1d.lambda)/3727. -1.  ;allows for variable end
      ; ; zmax = 0.4727
      ; pspace = 5            ; was 1
      ; width = 5*pspace
      ; nfind = 5
      ; npoly = 0             ; 3
      ; END TEST PARAMETERS


      for jj=0, n_elements(iobj)-1 do begin ;call one at a time
        eigenfile = result[i0].tfile
        if (not strmatch(eigenfile, '*QSO*',/fold)) then $ 
        ; TEST with same parameters (zmin, zmax) as initial zfind
          ; res1 = zfind(ss1d, hdr=hdr,objflux=objflux, objivar=objivar, $
          ; eigenfile=result[i0].tfile, columns=result[i0].tcolumn[0:ncol-1], $
          ; npoly=result[i0].npoly, zmin=zmin, zmax=zmax, loglam=loglam, $
          ; pspace=pspace, nfind=1, width=width, $
          ;  doplot=doplot*(jj eq 0), /nosubtract, _EXTRA=EXTRA)
        ; END TEST
          ; ; Original
          res1 = zfind(ss1d, hdr=hdr,objflux=objflux, objivar=objivar, $
          eigenfile=result[i0].tfile, columns=result[i0].tcolumn[0:ncol-1], $
          npoly=result[i0].npoly, zguess=result[indx[jj]].z, loglam=loglam, $
          pwidth=pwidth, pspace=pspace, nfind=1, width=width, $
           doplot=doplot*(jj eq 0), /nosubtract, fname=fname, _EXTRA=EXTRA) else $ 
;slightly different version of zfind for qso spectra
;
         res1 = zfind_qso(ss1d, hdr=hdr,objflux=objflux, objivar=objivar, $
          eigenfile=result[i0].tfile, columns=result[i0].tcolumn[0:ncol-1], $
          npoly=result[i0].npoly, zguess=result[indx[jj]].z, loglam=loglam, $
          pwidth=pwidth, pspace=pspace, nfind=1, width=width, $
           doplot=doplot*(jj eq 0), /nosubtract, _EXTRA=EXTRA)

      ; Copy the results into the output structure
         ; print, "res1 redshift: ", res1[*].z

        
         result[indx[jj]].z = res1[*].z
         result[indx[jj]].z_err = res1[*].z_err
         result[indx[jj]].rchi2 = res1[*].rchi2
         result[indx[jj]].dof = res1[*].dof
         result[indx[jj]].theta = res1[*].theta


         ; print, "columns=result[i0].tcolumn[0:ncol-1]:", result[i0].tcolumn[0:ncol-1]
         ; print, "res1 redshift: ", res1[*].z  
         ; print, "result redshift [jj]: ", result[indx[jj]].z

         ; print, "res1 rchi2: ", res1[*].rchi2
         ; print, "result rchi2 [jj]: ", result[indx[jj]].rchi2
 
      endfor

   endfor

   ; print, result
   return, result
end
;------------------------------------------------------------------------------