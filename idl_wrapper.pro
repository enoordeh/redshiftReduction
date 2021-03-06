pro idl_wrapper

;+
; PURPOSE
;      A wrapper routine to organize fit redshifts.
;      Reads in input spec, err spec and mask file
;      Uses zfind.pro to search for good redshifts
;      Uses zrefind.pro to refit the redshifts on finer gric to determine local minima.
;
; SYNTAX
;      idl_wrapper, [ result, files=files, /doplot, /debug, allresult= ]
;
; INPUTS
;      files = an optional parameter giving a vector of spec1d file
;              names. 
;
; KEYWORDS
;      /TODAY -- 
;
; OUTPUTS
;     result -- output structure giving result of last object in call
;               (or single object)
;   
; PROCEDURES CALLED 
;      splog -
;      zfind -
;      zrefind -
;      vdispfit
;      fill_gap
;      cpbackup
;      dfpsplot
;      concat_dir
;      remove_telluric
;      airtovac - find this...
;
; EXAMPLES
;
;
; COMMENTS
;
;
; HISTORY
;      Created March 8th 2016
;
;-

; files needed - main file and mask and errors

; input_dir = ''
; output_dir = ''
; eigen_dir = ''

; file1 = ''
; ext1 = 
; file2 = ''
; ext2 = 
; file3 = ''
; ext3 = 

; hdr_ext =

; ; input the files and find array size
; im=MRDFITS(file1,ext1,hdr1)
;   imsize=size(im)
;   xsize1=imsize[1]
; er=MRDFITS(file2,ext2,hdr2)
;   imsize=size(er)
;   xsize2=imsize[1]
; mask=MRDFITS(file3,ext3,hdr3)
;   imsize=size(mask)
;   xsize3=imsize[1]

; check the files are all there - TO DO!!!
; check files are same size
; is header info different?

; ; find the wavelength information
;   hdr2 = headfits(file, ext=hdr_ext, /silent)
;   start_w = strcompress( sxpar(hdr2, 'CRVAL1'), /rem)
;   increment_w = strcompress( sxpar(hdr2, 'CD1_1'), /rem)
;   ; make lambda array
;   lambda=(FINDGEN(xsize1, increment=1.0)*increment_w)+start_w

; ; Eigen dir directories and file names
;   eigenfile_gal = ''
;   eigenfile_star = ''
;   eigenfile_qso = ''
;  eigenfile = 'spEigenGal-*.fits'
;  eigenfile = 'emi_line.fits'

; write info to log file.
;  splog, 'Compute GALAXY redshifts: ', $
;    ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
; determine airmass value for the header of the 0th spec1d file.  - TO DO!!!
;  airmass = sxpar(hdr, 'AIRMASS')
;  splog, 'Mean airmass for this mask: ', airmass

; correct for telluric absorption bands (do this before shifting to - TO DO!!!
; the vacuum wavelengths!).
;          remove_telluric, ss1d, airmass
;	  ;fix_response,ss1d		CHANGED BL 7/22, fix_response or fix_response_blue is already done in fill_gap, no reason to do it again
          
; convert lambda values from air to vacuum. 
;          airtovac, lambda    
;          ss1d.lambda = float(lambda) ;replace, reconvert from double


debug=0
debug_gal = 0
debug_star = 0
debug_qso = 0
pspace_qso = 1   ; Normally run with 5. pspace_qso bug fixed now, can run at 1 (See XAGN doc for details)
physcheck=0 ; Check if linear combos are 'physical'? (0= no, 1 = full physcheck, 2 = light physcheck)
limitz=1 ; Apply secondary limitation on maximum redshift for each template fit to force spectrum to always be fully within the template bounds
nosubtract = 0
npolyall = 0 ; Add an (npolya ll - 1)th order polynomial to the template fitting procedure ; *** This appears to be causing issues with computechi2 (details in XAGN doc)
cont_sub_width = 200
trim=0 ; Trim data?
wvmin=10^3.58 ; Trim all input data to these wavelengths
wvmax=10^3.99
print, "wvmin = ", wvmin
print, "wvmax = ", wvmax
padding = -10000 ; -10000 for SDSS, -100 for PC

; Input files
; spec_file='/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/eigentest/objflux.csv'
; ivar_file='/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/eigentest/objivar.csv'
; wav_file='/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/eigentest/wavelengths.csv'

; SDSS Test Data Set
; spec_file='/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/objflux.csv'
; ivar_file='/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/objivar.csv'
; wav_file='/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/wavelengths.csv'

; ex='noneg_npoly0/run2/ex4'

; exs=['noneg_npoly0/run2/Telluric_Mask/ex5','noneg_npoly0/run2/Telluric_Mask/ex6','noneg_npoly0/run2/Telluric_Mask/ex7','noneg_npoly0/run2/Telluric_Mask/ex8','noneg_npoly0/run2/Telluric_Mask/ex9']
 
; exs = 'noneg_npoly0/ex4_pspace_qso_1_fixed'
; exs = '090A0958B'
; exs = ['090A0958B','092A0405A','094A0557A']
exs = 'all_formatted'
for i=0,N_ELEMENTS(exs)-1 do begin
ex=exs[i]
print, 'CURRENT folder: ', ex

; VIMOS SPECTRA
; spec_file='/scratch/emiln/VIMOS_AGN/'+ex+'/reflex_end_products/all_reductions/formatted/objflux.csv'
; ivar_file='/scratch/emiln/VIMOS_AGN/'+ex+'/reflex_end_products/all_reductions/formatted/objivar.csv'
; wav_file='/scratch/emiln/VIMOS_AGN/'+ex+'/reflex_end_products/all_reductions/formatted/wavelengths.csv'

; ALL VIMOS SPECTRA
spec_file='/scratch/emiln/VIMOS_AGN/'+ex+'/objflux.csv'
ivar_file='/scratch/emiln/VIMOS_AGN/'+ex+'/objivar.csv'
wav_file='/scratch/emiln/VIMOS_AGN/'+ex+'/wavelengths.csv'

; SDSS NON STACKED SPECTRA
; spec_file='/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/NotStacked/'+ex+'/objflux.csv'
; ivar_file='/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/NotStacked/'+ex+'/objivar.csv'
; wav_file='/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/NotStacked/'+ex+'/wavelengths.csv'

; Single SDSS Test spectrum that was poorly reduced
; spec_file='/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/bad_objflux.csv'
; ivar_file='/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/bad_objivar.csv'
; wav_file='/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/bad_wavelengths.csv'


; output database name


; VIMOS OUTPUTS
out_db='/scratch/emiln/VIMOS_AGN/'+ex+'/cont_width_500/outdb.csv' ; Best fit results
out_db2='/scratch/emiln/VIMOS_AGN/'+ex+'/cont_width_500/outdb2.csv' ; Holds smallest 2 rchi^2 for each class
a_out= '/scratch/emiln/VIMOS_AGN/'+ex+'/cont_width_500/a_coeff.csv' ; Holds the best fit [npolyall, a_coeff] including polynomial terms
a_out_AGN= '/scratch/emiln/VIMOS_AGN/'+ex+'/cont_width_500/a_coeff_AGN.csv' ; Holds the best fit AGN [npolyall, a_coeff]
a_out_GAL= '/scratch/emiln/VIMOS_AGN/'+ex+'/cont_width_500/a_coeff_GAL.csv' ; Holds the best fit GAL [npolyall, a_coeff]
a_out_STAR= '/scratch/emiln/VIMOS_AGN/'+ex+'/cont_width_500/a_coeff_STAR.csv' ; Holds the best fit STAR [npolyall, a_coeff]
z_out = '/scratch/emiln/VIMOS_AGN/'+ex+'/cont_width_500/z.csv' ; holds the top 8 redshifts corresopnding to top 8 classifications below
class_out = '/scratch/emiln/VIMOS_AGN/'+ex+'/cont_width_500/classes.csv' ; holds the top 8 classifications
tfile_out = '/scratch/emiln/VIMOS_AGN/'+ex+'/cont_width_500/tfile.csv' ; holds the eigentemplates associated with top 8 classifications


; VIMOS OUTPUTS
; out_db='/scratch/emiln/VIMOS_AGN/'+ex+'/reflex_end_products/all_reductions/formatted/outdb.csv' ; Best fit results
; out_db2='/scratch/emiln/VIMOS_AGN/'+ex+'/reflex_end_products/all_reductions/formatted/outdb2.csv' ; Holds smallest 2 rchi^2 for each class
; a_out= '/scratch/emiln/VIMOS_AGN/'+ex+'/reflex_end_products/all_reductions/formatted/a_coeff.csv' ; Holds the best fit [npolyall, a_coeff] including polynomial terms
; a_out_AGN= '/scratch/emiln/VIMOS_AGN/'+ex+'/reflex_end_products/all_reductions/formatted/a_coeff_AGN.csv' ; Holds the best fit AGN [npolyall, a_coeff]
; a_out_GAL= '/scratch/emiln/VIMOS_AGN/'+ex+'/reflex_end_products/all_reductions/formatted/a_coeff_GAL.csv' ; Holds the best fit GAL [npolyall, a_coeff]
; a_out_STAR= '/scratch/emiln/VIMOS_AGN/'+ex+'/reflex_end_products/all_reductions/formatted/a_coeff_STAR.csv' ; Holds the best fit STAR [npolyall, a_coeff]
; z_out = '/scratch/emiln/VIMOS_AGN/'+ex+'/reflex_end_products/all_reductions/formatted/z.csv' ; holds the top 8 redshifts corresopnding to top 8 classifications below
; class_out = '/scratch/emiln/VIMOS_AGN/'+ex+'/reflex_end_products/all_reductions/formatted/classes.csv' ; holds the top 8 classifications
; tfile_out = '/scratch/emiln/VIMOS_AGN/'+ex+'/reflex_end_products/all_reductions/formatted/tfile.csv' ; holds the eigentemplates associated with top 8 classifications


; Single SDSS Test spectrum that was poorly reduced
; out_db='/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/bad_outdb.csv' ; Best fit results
; out_db2='/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/bad_outdb2.csv' ; Holds smallest 2 rchi^2 for each class
; a_out= '/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/bad_a_coeff.csv' ; Holds the best fit [npolyall, a_coeff] including polynomial terms
; z_out = '/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/bad_z.csv' ; holds the top 8 redshifts corresopnding to top 8 classifications below
; class_out = '/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/bad_classes.csv' ; holds the top 8 classifications
; tfile_out = '/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/bad_tfile.csv' ; holds the eigentemplates associated with top 8 classifications

; SDSS Test Data set outputs
; out_db='/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/NotStacked/'+ex+'/outdb.csv' ; Best fit results
; out_db2='/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/NotStacked/'+ex+'/outdb2.csv' ; Holds smallest 2 rchi^2 for each class
; a_out= '/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/NotStacked/'+ex+'/a_coeff.csv' ; Holds the best fit [npolyall, a_coeff] including polynomial terms
; a_out_AGN= '/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/NotStacked/'+ex+'/a_coeff_AGN.csv' ; Holds the best fit AGN [npolyall, a_coeff]
; a_out_GAL= '/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/NotStacked/'+ex+'/a_coeff_GAL.csv' ; Holds the best fit GAL [npolyall, a_coeff]
; a_out_STAR= '/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/NotStacked/'+ex+'/a_coeff_STAR.csv' ; Holds the best fit STAR [npolyall, a_coeff]
; z_out = '/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/NotStacked/'+ex+'/z.csv' ; holds the top 8 redshifts corresopnding to top 8 classifications below
; class_out = '/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/NotStacked/'+ex+'/classes.csv' ; holds the top 8 classifications
; tfile_out = '/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/testdataset/SDSS_TEST_DATA/NotStacked/'+ex+'/tfile.csv' ; holds the eigentemplates associated with top 8 classifications

;mask_file='/scratch/emiln/VIMOS_AGN/090A0958B/reflex_end_products/test/mask.csv'
; readcol, spec_file, spec, format='F'
; readcol, ivar_file, ivar, format='F'
; readcol, wav_file, wav, format='F'
; ;readcol, mask_file, mask, format='F'
; mask=FLTARR(N_ELEMENTS(spec))

; length=N_ELEMENTS(spec)
;no_spec=
;ss1d = replicate({spec:FLTARR(length), ivar:FLTARR(length), lambda:FLTARR(length), mask:FLTARR(length)}, no_spec)


allobjflux = READ_CSV(spec_file)
allobjivar = READ_CSV(ivar_file)
allobjwav = READ_CSV(wav_file)

; HELP, allobjflux, /STRUCTURES
; print, N_TAGS(allobjflux)
print, "Counting # of Spectra..."
ncols = N_TAGS(allobjflux) ; number of spectra to iterate over
; ncols = 15 ; Only look at first column
; print, SIZE(allobjflux.(0))

;number of spec - integer probably found using nspec=N_ELEMENTS(some_file_we_read_in)
n_spec=ncols
print, "Number of spectra: ", n_spec

; storage structures
; gal_struct = replicate({zf_z:0.0, zr_z:0.0, zr_chi2:0.0, zr_dof:0.0, zr_zerr:0.0, type:''}, n_spec)
; qso_struct = replicate({zf_z:0.0, zr_z:0.0, zr_chi2:0.0, zr_dof:0.0, zr_zerr:0.0, type:''}, n_spec)
; star_struct = replicate({zf_z:0.0, zr_z:0.0, zr_chi2:0.0, zr_dof:0.0, zr_zerr:0.0, type:''}, n_spec)
gal_struct = replicate({z_1:0.0, rchi2_1:0.0, z_2:0.0, rchi2_2:0.0, class:''}, n_spec) ; Store the top 2 redshifts and their reduced chi^2 values
star_struct = replicate({z_1:0.0, rchi2_1:0.0, z_2:0.0, rchi2_2:0.0, class:''}, n_spec) ; Store the top 2 redshifts and their reduced chi^2 values
CV_struct = replicate({z_1:0.0, rchi2_1:0.0, z_2:0.0, rchi2_2:0.0, class:''}, n_spec) ; Store the top 2 redshifts and their reduced chi^2 values
qso_struct = replicate({z_1:0.0, rchi2_1:0.0, z_2:0.0, rchi2_2:0.0, class:''}, n_spec) ; Store the top 2 redshifts and their reduced chi^2 values
a_struct = replicate({npoly:0.0, a1:0.0, a2:0.0, a3:0.0, a4:0.0, a5:0.0, a6:0.0, a7:0.0, a8:0.0}, n_spec)
a_struct_AGN = replicate({z:0.0, a1:0.0, a2:0.0, a3:0.0, a4:0.0, a5:0.0, a6:0.0, a7:0.0, a8:0.0}, n_spec)
a_struct_GAL = replicate({z:0.0, a1:0.0, a2:0.0, a3:0.0, a4:0.0, a5:0.0, a6:0.0, a7:0.0, a8:0.0}, n_spec)
a_struct_STAR = replicate({z:0.0, a1:0.0, a2:0.0, a3:0.0, a4:0.0, a5:0.0, a6:0.0, a7:0.0, a8:0.0}, n_spec)
result_struct = replicate({z:0.0, rchi2:0.0, dof:0.0, z_err:0.0, class:'', tfile:'', istar:0.0}, n_spec)
class_struct = replicate({c1:'', c2:'', c3:'', c4:'', c5:'', c6:'', c7:'', c8:''}, n_spec) ; Store top 8 classes
z_struct = replicate({z1:0.0, z2:0.0, z3:0.0, z4:0.0, z5:0.0, z6:0.0, z7:0.0, z8:0.0}, n_spec) ; Store top 8 redshifts
tfile_struct = replicate({t1:'', t2:'', t3:'', t4:'', t5:'', t6:'', t7:'', t8:''}, n_spec) ; Store eigentemplates of top 8 fits


for ii=0, ncols-1 do begin
print, '#####################################################'
print, 'Processing Spectrum ', ii, ' out of ', n_spec
print, '#####################################################'
; ii=444
; some unique identifier for the spec
spec_id=ii

spec = allobjflux.(ii)
ivar = allobjivar.(ii)
wav = allobjwav.(ii)

; Filter out zero padding (-10000 padding) from test dataset
spec = spec[WHERE(spec GT padding, /NULL)]
; print, "objflux: ", spec
ivar = ivar[WHERE(ivar GT padding, /NULL)]
; print, "ivar: ", ivar
wav = wav[WHERE(wav GT padding, /NULL)]
; print, "wavelengths: ", wav




if trim gt 0 then begin

plot, wav, spec, color=FSC_color('white')
; pause

  diff = abs(wav - wvmin)
  minval = min(diff,minsub)
  minpix = (minsub - 1) > 0
  diff = abs(wav - wvmax)
  minval = min(diff,minsub)
  maxpix = (minsub - 1) > 0

  wav = wav[minpix:maxpix]
  spec = spec[minpix:maxpix]
  ivar = ivar[minpix:maxpix]


oplot, wav, spec, color=FSC_color('red')
; pause
endif



; print, "spec: ", spec
mask=FLTARR(N_ELEMENTS(spec))
length=N_ELEMENTS(spec)

print, "length(spec): ", N_ELEMENTS(spec)
print, "length(ivar): ", N_ELEMENTS(ivar)
print, "length(wav): ", N_ELEMENTS(wav)

ss1d = {spec:FLTARR(length), ivar:FLTARR(length), lambda:FLTARR(length), mask:FLTARR(length)}

ss1d.spec=spec
ss1d.ivar=ivar
ss1d.lambda=wav
ss1d.mask=mask

; print, wav
plot, ss1d.lambda, ss1d.spec

; estimate signal-to-noise in 1-d spectrum.
;          s2n = mean(ss1d.spec*sqrt(ss1d.ivar))
;          splog, ' '
;          splog, 'mean S2N for frame: ', files[ii], ' : ', s2n

;;; SMOOTH THE 1d SPECTRUM - needed to ease interpolation.
;   should we do this before we fill_gap?
;  WE were getting rchi2<1 with this smoothing.  Eliminate it!!
;    ss1d.spec = ivarsmooth(ss1d.spec, ss1d.ivar, 3, ss1d.ivar)
          
; check for NaNs and Infs...interpolate over them.
          nfn = where(finite(ss1d.spec) eq 0, nfnum, compl=fn)
          if nfnum gt 0 then begin 
              ss1d.spec[nfn] = interpol(ss1d.spec[fn], ss1d.lambda[fn], $
                                        ss1d.lambda[nfn])
          endif


; ------------------------------
; determine galaxy redshifts.
; first set the parameters for call to zfind.

; SDSS Test 2012
      eigenfile_gal='spEigenGal-55740.fits'
      eigenfile_CV = 'spEigenCVstar-55734.fits'
      eigenfile_star = 'spEigenStar-55734.fits'
      eigenfile_qso = 'spEigenQSO-55732.fits'

; print, "********", eigenfile_gal

; SDSS Test 2008
      ; eigenfile_gal='spEigenGal-53724.fits'
      ; eigenfile_CV = 'spEigenCVstar-55734.fits'
      ; eigenfile_star = 'spEigenStar-55734.fits' ; Should be spEigenStar-54474.fits but we don't have it. Doesn't make a big difference. Only 5 stars in dataset.
      ; eigenfile_qso = 'spEigenQSO-53724.fits'

; Big BOSS
		  ; eigenfile_gal='spEigenGal-51838.fits'
    ;   eigenfile_CV = 'spEigenCVstar-55734.fits'
    ;   eigenfile_star = 'spEigenStar-51845.fits'
    ;   eigenfile_qso = 'spEigenQSO-52223.fits'
		  eigendir='/scratch/bcanning/VIMOS_AGN/sdss12_templates/software/idlspec2d-v5_8_0/templates/'
          zmin = -0.0001
          zmax = max(ss1d.lambda)/3727. -1.  ;allows for variable end
          ; zmax = 12
          pspace = 5            ; was 1
          width = 5*pspace
          nfind = 5
          npoly = npolyall            ; 3  *** Changed by Emil (was 0)
          tclass = 'GALAXY'

          ; columns = Column numbers of the eigenspectra image to use in the
;                PCA fit; default to all columns.
; *** Should these not be rows? Eigentemplate has 4 rows of top eigenspectra, not columns...
;*          columns = [0,  1, 2]     ;include absorption and emission templates in all cases ;*** Why only include 3 eigenspectra of the 4?
          
; determine the redshift of object by comapring to the galaxy
; templates.
print, "zmax = ", zmax
fname='debug/zfind.dat' ; Debug file name

          resulti = zfind(ss1d, /linear_lambda, eigenfile=eigenfile_gal, $
                         eigendir=eigendir, npoly=npoly, $
                         zmin=zmin, zmax=zmax, nfind=nfind, pspace=pspace, $
                         width=width, plottitle=plottitle, doplot=doplot, $
                         debug=debug_gal, objflux=objflux, objivar=objivar, cont_sub_width=cont_sub_width, $
                         loglam=loglam, columns=columns, fname=fname, limitz=limitz, tclass=tclass, physcheck=physcheck, nosubtract=nosubtract)
          ; fname=fname, wvmin=wvmin, wvmax=wvmax)

print, "******************************************************************************************************"
print, "******************************************************************************************************"
print, "Initial z values: ", resulti.z   
print, "Initial chi^2 values: ", resulti.rchi2 ;*** These are per dof
print, "Initial DoF values: ", resulti.dof 
; print, "Initial z errors: ", result.z_err ; *** what do these mean?
; print, "Result.tfile: ", result.tfile ;*** array of template spectra used in top 5 chi^2 cases
; print, "Result.tcolumn: ", result.tcolumn 
print, "******************************************************************************************************"
print, "******************************************************************************************************"

; s

;print, result          
; redo the analysis at higher sampling, just around 5 minima chisqr regions
           pspace = 1 ; New sampling width in pixels for higher resolution fit
            splog, 'Locally re-fitting GALAXY redshifts'
            fname='debug/zrefind.dat' ; Debug file name
            res_gal = zrefind(ss1d, objflux, objivar, hdr=hdr, pwidth=15, $ ; look 15 pixels around each old peak for a better fit
                              pspace=pspace, width=3.*pspace, zold=resulti, $
                             loglam=loglam, plottitle=plottitle,  $
                             doplot=doplot, debug=debug_gal, columns=columns, fname=fname, tclass=tclass, physcheck=physcheck, nosubtract=nosubtract)


; ; calculate the delta chi^2 values from minimum to minimum.
          res_gal.class = 'GALAXY' 
           res_gal.subclass = ' '
           result = res_gal
           ; delta_chisqr = (result[1].rchi2-result[0].rchi2)*.5* $
           ;   (result[1].dof + result[0].dof)
           ; splog, 'Minimum chisqr is better than next min. by: ', delta_chisqr
           ; result[0].rchi2diff = result[1].rchi2-result[0].rchi2
           

print, "******************************************************************************************************"
print, "******************************************************************************************************"
print, "Zrefind z values: ", res_gal.z   ; *** nothing changing?
; print, "Zrefind z values (res_gal): ", res_gal.z   ; *** nothing changing?
print, "Zrefind chi^2 values: ", res_gal.rchi2 ; *** Running zrefind gives same z values but increases chi2? why...?
; print, "Zrefind chi^2 values (res_gal): ", res_gal.rchi2 ; *** Running zrefind gives same z values but increases chi2? why...?
print, "Zrefind DoF values: ", res_gal.dof 
; print, "Zrefind z errors: ", result.z_err ; *** what do these mean?
; print, "Zrefind object class: ", result.class ; *** These were just stated to be GALAXY above...
print, "******************************************************************************************************"
print, "******************************************************************************************************"
          


; Find CV STAR redshifts - do all stellar templates with one call.
; *** when running galaxy templates as a test, star template is giving best fit sometimes...

          npoly = npolyall             ; With only 1 eigen-template, fit 3 poly terms as well. *** Changed by Emil (was 0)
          zmin = -0.004         ; -1200 km/sec
          zmax = 0.004          ; +1200 km/sec
          pspace = 1
          nfind = 2
          tclass = 'CV'


resulti_CV = zfind_CV(ss1d, /linear_lambda, eigenfile=eigenfile_CV, $
                         eigendir=eigendir, npoly=npoly, $
                         zmin=zmin, zmax=zmax, nfind=nfind, pspace=pspace, $
                         width=width, plottitle=plottitle, doplot=doplot, $
                         debug=debug_star, objflux=objflux, objivar=objivar, cont_sub_width=cont_sub_width, $
                         loglam=loglam, columns=columns, fname=fname, tclass=tclass, physcheck=physcheck, nosubtract=nosubtract)
          ; fname=fname, wvmin=wvmin, wvmax=wvmax)

print, "******************************************************************************************************"
print, "******************************************************************************************************"
print, "Initial z values (CV): ", resulti_CV.z   
print, "Initial chi^2 values (CV): ", resulti_CV.rchi2 ;*** These are per dof
print, "Initial DoF values (CV): ", resulti_CV.dof 
; print, "Initial z errors: ", result.z_err ; *** what do these mean?
; print, "Result.tfile: ", result.tfile ;*** array of template spectra used in top 5 chi^2 cases
; print, "Result.tcolumn: ", result.tcolumn 
print, "******************************************************************************************************"
print, "******************************************************************************************************"

; s

;print, result          
; redo the analysis at higher sampling, just around 5 minima chisqr regions
           pspace = 1 ; New sampling width in pixels for higher resolution fit
            splog, 'Locally re-fitting GALAXY redshifts'
            fname='debug/zrefind.dat' ; Debug file name
            res_CV = zrefind_CV(ss1d, objflux, objivar, hdr=hdr, pwidth=15, $ ; look 15 pixels around each old peak for a better fit
                              pspace=pspace, width=5.*pspace, zold=resulti_CV, $
                             loglam=loglam, plottitle=plottitle,  $
                             doplot=doplot, debug=debug_star, columns=columns, fname=fname, tclass=tclass, physcheck=physcheck, nosubtract=nosubtract)


; ; calculate the delta chi^2 values from minimum to minimum.
          res_CV.class = 'CV' ; *** This is just setting the class to galaxy? what? b/c this rerun is only for eigenfile_gal
           res_CV.subclass = ' '
           result = [result, res_CV]
           ; delta_chisqr = (result[1].rchi2-result[0].rchi2)*.5* $
           ;   (result[1].dof + result[0].dof)
           ; splog, 'Minimum chisqr is better than next min. by: ', delta_chisqr
           ; result[0].rchi2diff = result[1].rchi2-result[0].rchi2
           

print, "******************************************************************************************************"
print, "******************************************************************************************************"
print, "Zrefind z values (CV): ", res_CV.z   ; *** nothing changing?
; print, "Zrefind z values (res_gal): ", res_gal.z   ; *** nothing changing?
print, "Zrefind chi^2 values (CV): ", res_CV.rchi2 ; *** Running zrefind gives same z values but increases chi2? why...?
; print, "Zrefind chi^2 values (res_gal): ", res_gal.rchi2 ; *** Running zrefind gives same z values but increases chi2? why...?
print, "Zrefind DoF values (CV): ", res_CV.dof 
; print, "Zrefind z errors: ", result.z_err ; *** what do these mean?
; print, "Zrefind object class: ", result.class ; *** These were just stated to be GALAXY above...
print, "******************************************************************************************************"
print, "******************************************************************************************************"







     
; now find the velocity dispersions aacording to the galaxy models.
;            splog, 'Find velocity dispersions for galaxies'
;          vdisp_wrapper, files[ii], zresult=result, $
;            nfind=nfind, airmass=airmass, nlsky=nlsky
;         
;            splog, 'CPU time to fit GALAXY velocity dispersions = ', systime(1)-t0
          


; Find STAR redshifts - do all stellar templates with one call.
; *** when running galaxy templates as a test, star template is giving best fit sometimes...

          npoly = npolyall             ; With only 1 eigen-template, fit 3 poly terms as well. *** Changed by Emil (was 0)
          zmin = -0.004         ; -1200 km/sec
          zmax = 0.004          ; +1200 km/sec
          pspace = 1
          nfind = 1
          tclass = 'STAR'
; check the system time.
          ts0 = systime(1)
            
;            splog, 'Compute STAR (' + subclass + ') redshifts:', $
;              ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
          res_star = zfind_star(ss1d, /linear_lambda, $
                                eigenfile=eigenfile_star,eigendir=eigendir, npoly=npoly, $
                                zmin=zmin, zmax=zmax, pspace=1, $
                                nfind=nfind, width=5*pspace, cont_sub_width=cont_sub_width, $
                                doplot=0, debug=debug_star, tclass=tclass, physcheck=physcheck, nosubtract=nosubtract) ;doplot=doplot
          
          res_star.class = 'STAR'
          res_star = res_star[sort(res_star.rchi2)] ;sort by rchi2 
          
          print, "res_star.z :", res_star.z
          result = [result, res_star[0:2]] ; Append result of top 3 star
          istars = res_star.tcolumn[0]

print, "******************************************************************************************************"
print, "******************************************************************************************************"
print, "Zfind z values for galaxies and stars: ", result.z   ;  
print, "Zfind chi^2 values for galaxies and stars ", result.rchi2 
print, 'istar : ', istars[0]
; print, 'istar = ', res_star.istar[0:2]
print, "******************************************************************************************************"
print, "******************************************************************************************************" 



 ;----------
   ; Find QSO redshifts


          npoly = npolyall ; *** Changed by Emil (was 0)
          zmin = 0.0033         ; +1000 km/sec
;          zmax = max(ss1d.lambda)/1215. -1  ; Max range we can expect to see
          zmax = 5.
          pspace = pspace_qso
          nfind = 2             ;find best QSO candidate z
          plottitle = 'QSO Redshift'
          tclass='AGN'

          

          ; splog, 'Compute QSO redshifts:', $
          ;        ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
          ; t0 = systime(1)
          
          ; print, "eigenfile QSO: ", eigenfile_qso
          ; wait, 5
          res_qso =  zfind_qso(ss1d, /linear_lambda, $
                      eigenfile = eigenfile_qso, eigendir=eigendir, npoly = npoly, zmin = zmin, $
                      zmax = zmax, pspace = pspace, loglam = loglam, $
                      nfind = nfind, width = 25, objflux = objflux, $ ; width changed from 5*pspace to 25
                      objivar=objivar, plottitle = plottitle, cont_sub_width=cont_sub_width, $
                       doplot = doplot, debug = debug_qso, tclass=tclass, physcheck=physcheck, nosubtract=nosubtract)


          ; splog, 'CPU time to compute QSO redshifts = ', systime(1)-t0

          res_qso.class = 'AGN'
          res_qso.subclass = ' '


          result = [result, res_qso[0:1]] ; Append result of top 2 QSO *Emil

print, "******************************************************************************************************"
print, "******************************************************************************************************"
print, "Zfind QSO z values: ", res_qso.z   ; *** nothing changing?
print, "Zfind QSO chi^2 values: ", res_qso.rchi2 ; *** Running zrefind gives same z values but increases chi2? why...?
print, "Zfind QSO DoF values: ", res_qso.dof 
print, "Zfind QSO theta ", res_qso.theta
print, "******************************************************************************************************"
print, "******************************************************************************************************"


;skip refitting QSO', as not necessary
;
;          splog, 'Locally re-fitting QSO redshifts'
;         t0 = systime(1)

;          reres_qso = zrefind(ss1d, objflux, objivar, hdr=hdr, pwidth=91, $ ; Emil - why pwidth=91?
;                            pspace=1, width=3.*pspace, zold=res_qso, $
;                            loglam=loglam, plottitle=plottitle, $
;                            doplot=doplot, debug=debug  )
; ;
;
;          splog, 'CPU time to re-fit QSO redshifts = ', systime(1)-t0

          ; reres_qso.class = 'AGN'
          ; reres_qso.subclass = ' '

; don't append QSO results here. we don't want them to be sorted by
; chi2 value. instead, we will sort the GALAXY and STAR z values and
; then append the 2 QSO redshifts.
;          result = [result, res_qso] ; Append results

; print, "******************************************************************************************************"
; print, "******************************************************************************************************"
; print, "Zrefind QSO z values: ", reres_qso.z   ; *** nothing changing?
; print, "Zrefind QSO chi^2 values: ", reres_qso.rchi2 ; *** Running zrefind gives same z values but increases chi2? why...?
; print, "Zrefind QSO DoF values: ", reres_qso.dof 
; print, "******************************************************************************************************"
; print, "******************************************************************************************************"

; Sort resulti, res_gal, res_star, and res_qso independetly for data storage purposes

; now decide which result gives best rchisqr
          
;----------
          
; Individually sort the results for galaxies, stars, and quasars

isort_gal = sort(res_gal.rchi2)
res_gal_sorted = res_gal[isort_gal]

isort_star = sort(res_star.rchi2)
res_star_sorted = res_star[isort_star]

isort_qso = sort(res_qso.rchi2)
res_qso_sorted = res_qso[isort_qso]

isort_CV = sort(res_CV.rchi2)
res_CV_sorted = res_CV[isort_CV]



   ;----------
   ; Sort results for each object by ascending order in chi^2/DOF,
   ; but putting any results with zero degrees-of-freedom at the end.

          minvdiff = 1000.0     ; km/s
          cspeed = 2.99792458e5
          rchi2 = result.rchi2
          
          isort = sort(rchi2 + (result.dof EQ 0)*max(rchi2))
          result = result[isort]

; append the QSO results here! -- so that they aren't sorted by chi2! *** Why?
          ; result = [result, res_qso] ; *Emil - res_qso stored in result and sorted with others for us
          nper = (size(result,/dimens))[0]
; Find the difference in reduced chi^2 between each result and the next
          rchi2 = result.rchi2
          for ia=0,nper-2 do begin
              inext = (where(abs(result[ia+1:nper-1].z - result[ia].z) GT $
                        minvdiff/cspeed AND result[ia+1:nper-1].dof GT 0))[0]
              if (inext ne -1) then $
                result[ia].rchi2diff = rchi2[ia+1+inext] - rchi2[ia]
          endfor


; define how many results to keep.
      nres = (n_elements(result) < 10) - 1          

; form result structure for output-- will evolve as Z's become real
      if not(keyword_set(cresult)) then begin
          cresult = result[0]   ;if first in list
          allresult = result[0:nres] ;save best nres fits
      endif else begin
          cresult = [[cresult], [result[0]]] ;if not first, append
          allresult = [[allresult], [result[0:nres]]]
      endelse

; print, cresult
; print, allresult

; Fill in data structures
print, "result.z", result.z
print, "result.rchi2", result.rchi2
print, "result.dof", result.dof
print, "result.z_err", result.z_err
print, "result.class", result.class
print, "result.theta", result.theta

; Save everything -- This needs to be modified to have many more columns to hold above results
; gal_struct[ii].zf_z=result.z ; here just tag the results that are returned to us in the result structures
; gal_struct[ii].zr_z=res_gal.z
; gal_struct[ii].zr_chi2=res_gal.chi2
; gal_struct[ii].zr_dof=res_gal.dof
; gal_struct[ii].zr_zerr=res_gal.zerr
; gal_struct[ii].type='GALAXY' ; or maybe you want name of eigentemplate or something?

; qso_struct[ii].zf_z=res_qso.z
; qso_struct[ii].zf_chi2=res_qso.chi2
; qso_struct[ii].zf_dof=res_qso.dof
; qso_struct[ii].zf_zerr=res_qso.zerr
; qso_struct[ii].type='QSO' ; or maybe you want name of eigentemplate or something?

; star_struct[ii].zf_z=res_star.z
; star_struct[ii].zf_chi2=res_star.chi2
; star_struct[ii].zf_dof=res_star.dof
; star_struct[ii].zf_zerr=res_star.zerr
; star_struct[ii].type='STAR' ; or maybe you want name of eigentemplate or something?

; Just save best fit result after sorting gal, star, qso by chi^2
print, 'result[0].z: ', result[0].z
print, 'result[0].theta: ', result[0].theta
print, 'result[0].tfile: ', result[0].tfile
print, 'npolyall=', npolyall
result_struct[ii].z=result[0].z ; here just tag the results that are returned to us in the result structures
result_struct[ii].rchi2=result[0].rchi2
result_struct[ii].dof=result[0].dof
result_struct[ii].z_err=result[0].z_err
result_struct[ii].class=result[0].class ; or maybe you want name of eigentemplate or something?
result_struct[ii].tfile = result[0].tfile
result_struct[ii].istar = istars[0] ; Stores the (row# - 1) of best fit star template

; Also save the top 2 rchi^2 and respective redshifts for each class
gal_struct[ii].z_1 = res_gal_sorted[0].z
gal_struct[ii].rchi2_1 = res_gal_sorted[0].rchi2
gal_struct[ii].z_2 = res_gal_sorted[1].z
gal_struct[ii].rchi2_2 = res_gal_sorted[1].rchi2

star_struct[ii].z_1 = res_star_sorted[0].z
star_struct[ii].rchi2_1 = res_star_sorted[0].rchi2
star_struct[ii].z_2 = res_star_sorted[1].z
star_struct[ii].rchi2_2 = res_star_sorted[1].rchi2

qso_struct[ii].z_1 = res_qso_sorted[0].z
qso_struct[ii].rchi2_1 = res_qso_sorted[0].rchi2
qso_struct[ii].z_2 = res_qso_sorted[1].z
qso_struct[ii].rchi2_2 = res_qso_sorted[1].rchi2

CV_struct[ii].z_1 = res_CV_sorted[0].z
CV_struct[ii].rchi2_1 = res_CV_sorted[0].rchi2
CV_struct[ii].z_2 = res_CV_sorted[1].z
CV_struct[ii].rchi2_2 = res_CV_sorted[1].rchi2

a_struct[ii].npoly = npolyall
a_struct[ii].a1 = result[0].theta[0]
a_struct[ii].a2 = result[0].theta[1]
a_struct[ii].a3 = result[0].theta[2]
a_struct[ii].a4 = result[0].theta[3]
a_struct[ii].a5 = result[0].theta[4]
a_struct[ii].a6 = result[0].theta[5]
a_struct[ii].a7 = result[0].theta[6]
a_struct[ii].a8 = result[0].theta[7]

a_struct_AGN[ii].z = res_qso_sorted[0].z
a_struct_AGN[ii].a1 = res_qso_sorted[0].theta[0]
a_struct_AGN[ii].a2 = res_qso_sorted[0].theta[1]
a_struct_AGN[ii].a3 = res_qso_sorted[0].theta[2]
a_struct_AGN[ii].a4 = res_qso_sorted[0].theta[3]
a_struct_AGN[ii].a5 = res_qso_sorted[0].theta[4]
a_struct_AGN[ii].a6 = res_qso_sorted[0].theta[5]
a_struct_AGN[ii].a7 = res_qso_sorted[0].theta[6]
a_struct_AGN[ii].a8 = res_qso_sorted[0].theta[7]

a_struct_GAL[ii].z = res_gal_sorted[0].z
a_struct_GAL[ii].a1 = res_gal_sorted[0].theta[0]
a_struct_GAL[ii].a2 = res_gal_sorted[0].theta[1]
a_struct_GAL[ii].a3 = res_gal_sorted[0].theta[2]
a_struct_GAL[ii].a4 = res_gal_sorted[0].theta[3]
a_struct_GAL[ii].a5 = res_gal_sorted[0].theta[4]
a_struct_GAL[ii].a6 = res_gal_sorted[0].theta[5]
a_struct_GAL[ii].a7 = res_gal_sorted[0].theta[6]
a_struct_GAL[ii].a8 = res_gal_sorted[0].theta[7]

a_struct_STAR[ii].z = res_star_sorted[0].z
a_struct_STAR[ii].a1 = res_star_sorted[0].theta[0]
a_struct_STAR[ii].a2 = res_star_sorted[0].theta[1]
a_struct_STAR[ii].a3 = res_star_sorted[0].theta[2]
a_struct_STAR[ii].a4 = res_star_sorted[0].theta[3]
a_struct_STAR[ii].a5 = res_star_sorted[0].theta[4]
a_struct_STAR[ii].a6 = res_star_sorted[0].theta[5]
a_struct_STAR[ii].a7 = res_star_sorted[0].theta[6]
a_struct_STAR[ii].a8 = res_star_sorted[0].theta[7]

class_struct[ii].c1 = result[0].class
class_struct[ii].c2 = result[1].class
class_struct[ii].c3 = result[2].class
class_struct[ii].c4 = result[3].class
class_struct[ii].c5 = result[4].class
class_struct[ii].c6 = result[5].class
class_struct[ii].c7 = result[6].class
class_struct[ii].c8 = result[7].class

z_struct[ii].z1 = result[0].z
z_struct[ii].z2 = result[1].z
z_struct[ii].z3 = result[2].z
z_struct[ii].z4 = result[3].z
z_struct[ii].z5 = result[4].z
z_struct[ii].z6 = result[5].z
z_struct[ii].z7 = result[6].z
z_struct[ii].z8 = result[7].z

tfile_struct[ii].t1 = result[0].tfile
tfile_struct[ii].t2 = result[1].tfile
tfile_struct[ii].t3 = result[2].tfile
tfile_struct[ii].t4 = result[3].tfile
tfile_struct[ii].t5 = result[4].tfile
tfile_struct[ii].t6 = result[5].tfile
tfile_struct[ii].t7 = result[6].tfile
tfile_struct[ii].t8 = result[7].tfile

endfor
; print, "result_struct.z", result_struct.z
; print, "result_struct.rchi2", result_struct.rchi2
; print, "result_struct.dof", result_struct.dof
; print, "result_struct.z_err", result_struct.z_err
; print, "result_struct.class", result_struct.class
; print ALL to file
; writecol, out_db, spec_id, $
; gal_struct.zf_z, gal_struct.zr_z, gal_struct.zr_chi2, gal_struct.zr_dof, gal_struct.zr_zerr, gal_struct.type, $
; star_struct.zr_z, star_struct.zr_chi2, star_struct.zr_dof, star_struct.zr_zerr, star_struct.type, $
; qso_struct.zr_z, qso_struct.zr_chi2, qso_struct.zr_dof, qso_struct.zr_zerr, qso_struct.type, $
; fmt=('f,f,f,f,f,f,s,f,f,f,f,s,f,f,f,f,s'), delimit=','

; print best fit to file
; writecol, out_db, $
; result_struct.z, result_struct.rchi2, result_struct.dof, result_struct.z_err, result_struct.class, fmt='(f,f,f,f,a)'

; WRITE_CSV takes a maximum of 8 columns of data,

; output_data = [result_struct.z, result_struct.rchi2, result_struct.dof, result_struct.z_err, result_struct.class, gal_struct.rchi2_1, $
;  gal_struct.rchi2_2, star_struct.rchi2_1, star_struct.rchi2_2, qso_struct.rchi2_1, qso_struct.rchi2_2] ; Not working - can't have a 2D array with floats AND strings?
; WRITE_CSV, out_db, output_data
WRITE_CSV, out_db, result_struct
; WRITE_CSV, out_db, result_struct.z, result_struct.rchi2, result_struct.dof, result_struct.z_err, result_struct.class
WRITE_CSV, out_db2, gal_struct.rchi2_1, gal_struct.rchi2_2, star_struct.rchi2_1, star_struct.rchi2_2, qso_struct.rchi2_1, qso_struct.rchi2_2, CV_struct.rchi2_1, CV_struct.rchi2_2
WRITE_CSV, a_out, a_struct
WRITE_CSV, a_out_AGN, a_struct_AGN
WRITE_CSV, a_out_GAL, a_struct_GAL
WRITE_CSV, a_out_STAR, a_struct_STAR
WRITE_CSV, z_out, z_struct
WRITE_CSV, class_out, class_struct
WRITE_CSV, tfile_out, tfile_struct

endfor
end

 