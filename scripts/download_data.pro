;;==============================================================================
;; Author: Michael Palumbo
;; Date: July 2017
;; Purpose: Download SDO data from VSO
;;==============================================================================
PRO DOWNLOAD_SDO, rawdir
    time1 = '2014/01/01 00:00:00'
    time2 = '2014/12/31 24:00:00'

    ;;==============================================================================
    ;; *** Ensure that rawdir is empty ***
    ;;==============================================================================
    IF FILE_TEST(rawdir + '*.fits') EQ 1 THEN BEGIN
        SPAWN, 'rm ' + rawdir + '*.fits'
    ENDIF

    ;;==============================================================================
    ;; *** Query VSO for files and download  ***
    ;;==============================================================================
    PRINT, "Downloading files"
    dop = vso_search(time1, time2, instr='hmi', provider='JSOC', physobs='LOS_velocity', sample=14400)
    a = vso_get(dop, out_dir=rawdir, /QUIET, /NOWAIT)

    mag = vso_search(time1, time2, instr='hmi', provider='JSOC', physobs='LOS_magnetic_field', sample=14400)
    b = vso_get(mag, out_dir=rawdir, /QUIET, /NOWAIT)

    int = vso_search(time1, time2, instr='hmi', provider='JSOC', physobs='intensity', sample=14400)
    c = vso_get(int, out_dir=rawdir, /QUIET, /NOWAIT)
    
    aia = vso_search(time1, time2, instr='aia', provider='JSOC', physobs='intensity', wave='1700', sample=14400)
    d = vso_get(aia, out_dir=rawdir, /QUIET, /NOWAIT)
    
    ;;==============================================================================
    ;; *** End of program ***
    ;;==============================================================================
END
