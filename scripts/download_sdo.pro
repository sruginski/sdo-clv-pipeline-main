;;==============================================================================
;; Author: Michael Palumbo
;; Date: July 2017
;; Purpose: Download SDO data from VSO
;;==============================================================================
PRO DOWNLOAD_SDO, rawdir
    time1 = '2017/01/01 00:00:00'
    time2 = '2017/01/01 12:00:00'

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
    dop = vso_search(time1, time2, instr='hmi', provider='JSOC', physobs='LOS_velocity')
    a = vso_get(dop, out_dir=rawdir, /QUIET, /NOWAIT)

    mag = vso_search(time1, time2, instr='hmi', provider='JSOC', physobs='LOS_magnetic_field')
    b = vso_get(mag, out_dir=rawdir, /QUIET, /NOWAIT)

    int = vso_search(time1, time2, instr='hmi', provider='JSOC', physobs='intensity')
    c = vso_get(int, out_dir=rawdir, /QUIET, /NOWAIT)
    
    aia = vso_search(time1, time2, instr='aia', provider='JSOC', wave='1700', physobs='intensity', sample=45)
    d = vso_get(aia, out_dir=rawdir, /QUIET, /NOWAIT)
    
    ;;==============================================================================
    ;; *** End of program ***
    ;;==============================================================================
END
