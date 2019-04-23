!dis !dis    Open Source License/Disclaimer, Forecast Systems Laboratory
!dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
!dis
!dis    This software is distributed under the Open Source Definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis
!dis    In particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis
!dis    - Redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis
!dis    - Redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis
!dis    - All modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis
!dis    - If significant modifications or enhancements are made to this
!dis    software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis
!dis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
!dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
!dis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
!dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
!dis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
!dis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
!dis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
!dis
!dis

!WRF:PACKAGE:IO

     MODULE EMISSIONS_WPS_V4

!
! Version 2.0 of the emissions standard initialization routine. This
!  program is for use with the chemistry code in WRF V3
!
! This program takes formatted output of em11v1 NEI-2011 v1 inventory and
! grids to a different grid.  Simple grid dumping done here so it is
! not to be used with domains using grid spacing less than 12 km.
!
!  Stu McKeen 6/20/04
!  Steven Peckham 1/25/05
!
! compile with:
!  pgf90 -w -byteswapio -Mfree -Mlfs emiss_v03_wps.F
!
!or
!
!  mpif90 -axP -free -convert big_endian emiss_v03_wps.F
!-----------------------------------------------------------------------
!
!  Fields to set before running:
!
!  zfa                    elevation at grid cell top (m)
!  ix2                    x-dimension of output  data (ix2=nx-1 for WRF domain)
!  jx2                    y-dimension of output  data (jx2=ny-1 for WRF domain)
!  kx                     z-dimension of output  data (kx =kemit for WRF domain)
!  IPROJ                  =1 if Lambert Conformal, =2 if Polar Stereographic
!  DX                     horizontal grid spacing (m)
!  REKM                   Earth radius (km)
!  XLATC                  Center latitude of mother doamin projection
!  XLONC                  Center longitude of mother domain projection (-180->180E)
!  CLAT1                  Northern most reference latitude of mother domain projection
!  CLAT2                  Southern most ref. lat. of mother domain for Lambert Conformal, not used for Polar Stereo (TRUELAT1 > TRUELAT2)
!  INEST1                 =0 if emissions in mother domain, =1 if a nest within mother domain
!  XNESSTR                X value of southwestern most (dot) point of nested domain in mother (=1. if mother)
!  YNESSTR                Y value of southwestern most (dot) point of nested domain in mother (=1. if mother)
!  DXBIGDO                horizontal grid spacing (m) of mother domain (=DX if emissions in mother)
!  IL                     X (west-east) stagger dimension of mother (=IX2+1 if emissions in mother)
!  JL                     Y (south-north) stagger dimension of mother (=JX2+1 if emissions in mother)
!  HEMI                   1 for Northern Hemisphere, -1 for Southern Hemisphere
!  ISTART                 1 for starting at 00z , 12 for starting at 12z
!  ZFA                    grid level height (m) of 4-D emission file specified at the top of the WRF grid cell (w data point).
!  MAXHR                  Number of hours the emissions data is to be generated for (e.g., maxhr=3 for 3 hours of data)
!  POINDIR                directories of point emissions input data
!  AREADIR                directories of area emissions input data
!  KWIN                   Number of vertical levels in the horizontal wind profile
!  REFWZ                  interval grid level height (m) of wind-data (climatology) used in momentum lift calcs.
!  WSPD                   wind speed at level height (m) specified by REFWZ
!-----------------------------------------------------------------------
!
! - Input emissions data dimensions.  This is fixed for the 4-km grid
!    that the NEI-05 emissions is located on.
!  Note the grid numbering convention.  The lower left corner is numbered (1,1)
!    and the lower left cell center is numbered (1.5,1.5)
!
     PARAMETER(IX=1332,JX=1008,IP=IX+1,JP=JX+1,IM=IX-1,JM=JX-1)
!
! - Output emissions data dimensions (user specified)
!
     PARAMETER(IX2=99,JX2=99,KX=20,KP=KX+1)
!
!
! Number of release points
!
!    PARAMETER(IPOINT=168516)  ! NEI2005 version 2
     PARAMETER(IPOINT=224097) !  NEI-11 US, version1, and NEI05 Mexico and Canada 11/1/13
!
! Number of Primary (IPRIM), VOC (IVOC), and PM2.5 (IPM25) species in NEI-11 emission files
!
     PARAMETER(IPRIM=7, IVOC=50, IPM25=19)
     PARAMETER(NAL2DO=56, NRADM=31, NAMFILE2=64)
! 
! Number of grid levels, and interval levels for wind climatology data used in vertical momentum lift calcs.
! 
     PARAMETER(KWIN=20,KWINP=KWIN+1)
! 
     CHARACTER (len=80)                       ::  INFIL,OUFILU,FIL,FIL24
     CHARACTER (len= 9), DIMENSION(NAMFILE2)  ::  NAM, NAMRAD
     CHARACTER (len= 9), DIMENSION(NAL2DO)    ::  NAM2EM
     CHARACTER (len= 9), DIMENSION(NRADM)     ::  ENAME
     CHARACTER (len=10), DIMENSION(IPRIM)     ::  NAMINF
     CHARACTER (len=10), DIMENSION(IVOC)      ::  NAMVOC
     CHARACTER (len=10), DIMENSION(IPM25)     ::  NAMPM2
     CHARACTER (len=10)                       ::  INARG,CHSPEC,CHSCRT
     CHARACTER (len= 4)                       ::  HRNAM
     CHARACTER (len=128)                      ::  POINDIR,AREADIR

     INTEGER                                  ::  I, J, K, N
     INTEGER                                  ::  IDEX
     INTEGER                                  ::  IEMAX, JEMAX
     INTEGER, DIMENSION(NAMFILE2)             ::  MWT

     REAL, DIMENSION(IX,JX)     :: EMT
     REAL, DIMENSION(IPOINT)    :: EMP
     REAL, DIMENSION(IX2,JX2)   :: EM2D
     REAL, DIMENSION(24,KWIN)     :: WSPD

     REAL, DIMENSION(IX2,KX,JX2,NRADM)   :: EM3RD
     REAL, DIMENSION(IX2,KX,JX2)         :: EM3RS
     REAL, DIMENSION(NAMFILE2)           :: FAC
!
!     Output grid map information
!
     REAL, DIMENSION(IP,JP)    :: XLATD,XLOND
     REAL, DIMENSION(KP)       :: ZFA(KP)  ! interval grid level height (m) of 4-D emission file
     REAL, DIMENSION(KWINP)    :: REFWZ(KWINP)  ! interval grid level height (m) of wind-data (climatology) used in momentum lift calcs.
     REAL                      :: ETOT, EMAX, EMIN
     REAL                      :: ZTOP

!     Constants
     REAL       :: rad_per_deg = 0.017453293
     REAL       :: deg_per_rad = 57.29577951

!
! Map information for output grid
!  IPROJ                  =1 if Lambert Conformal, =2 if Polar Stereographic
!  DX                     horizontal grid spacing (m)
!  REKM                   Earth radius (km)
!  XLATC                  Center latitude of mother doamin projection
!  XLONC                  Center longitude of mother domain projection (-180->180E)
!  CLAT1                  Northern most reference latitude of mother domain projection
!  CLAT2                  Southern most ref. lat. of mother domain for Lambert Conformal, not used for Polar Stereo (TRUELAT1 > TRUELAT2)
!  INEST1                 =0 if emissions in mother domain, =1 if a nest within mother domain
!  XNESSTR                X value of southwestern most (dot) point of nested domain in mother (=1. if mother)
!  YNESSTR                Y value of southwestern most (dot) point of nested domain in mother (=1. if mother)
!  DXBIGDO                horizontal grid spacing (m) of mother domain (=DX if emissions in mother)
!  IL                     X (west-east) stagger dimension of mother (=IX2+1 if emissions in mother)
!  JL                     Y (south-north) stagger dimension of mother (=JX2+1 if emissions in mother)
!  HEMI                   1 for Northern Hemisphere, -1 for Southern Hemisphere
!  STARTHR                1 for starting at 00z , 12 for starting at 12z
!  MAXHR                  Number of hours the emissions data is to be generated for (e.g., maxhr=3 for 3 hours of data)
!
     INTEGER    :: hemi = +1  ! SAM 7/9/08 - Only northern hemisphere Lambert Conf. and Polar Stereograph tested
! Following for WRF/Chem 60 km domain (Lambert Conformal), no nest, eastern U.S. (November, 2003)
!     INTEGER    :: iproj = 1
!     REAL       :: rekm = 6371.
!     REAL       :: dx = 60.E3
!     REAL       :: xlatc =  40.00
!     REAL       :: xlonc = -115.00
!     REAL       :: clat1 = 45.00
!     REAL       :: clat2 = 35.00
!     INTEGER    :: inest1 = 0
!     REAL       :: xnesstr = 1.00
!     REAL       :: ynesstr = 1.00
!     REAL       :: dxbigdo = 60.E3
!     INTEGER    :: il = 134
!     INTEGER    :: jl = 110
! Following for WRF/Chem 4 km domain (Lambert Conformal), no nest, continental U.S. (Jun 2018)
     INTEGER    :: iproj = 1
     REAL       :: rekm = 6370.
     REAL       :: dx = 4.E3
     REAL       :: dxbigdo = 4.E3
     REAL       :: xlatc =  33.75
     REAL       :: xlonc = -88.75
     REAL       :: clat1 = 60.00
     REAL       :: clat2 = 30.00
     INTEGER    :: inest1 = 0
     REAL       :: xnesstr = 1.00
     REAL       :: ynesstr = 1.00
     INTEGER    :: il = 100
     INTEGER    :: jl = 100
     INTEGER    :: starthr = 01
     INTEGER    :: maxhr = 24
     INTEGER    :: endhr 
!
! Names of emission variables
!
     DATA ENAME /    &
      'e_so2 ','e_no  ','e_ald ','e_hcho','e_ora2',                     &
      'e_nh3 ','e_hc3 ','e_hc5 ','e_hc8 ',                              &
      'e_eth ','e_co  ','e_ol2 ','e_olt ','e_oli ','e_tol ','e_xyl ',   &
      'e_ket ','e_csl ','e_iso ','e_ch4','e_pm25i','e_pm25j',           &
      'e_so4i','e_so4j','e_no3i','e_no3j','e_orgi','e_orgj','e_eci',    &
      'e_ecj','e_pm10'/
! WARNING-WARNING-WARNING - Order sensitive in gas-phase species emission processing: 9/27/10
! Any new gas-phase species (i.e. hc16 or ch4) needs to be included in binary
! file before PM2.5 emissions. Order of binary file output is set-up in
! variable ENAME.  Notice convert_emiss also needs to modified, and order of
! binary file reads is significant in that subroutine as well.

     DATA NAM2EM /  &
       'CO', 'NOX', 'SO2', 'NH3','HC01','HC02','HC03','HC04','HC05',    &
     'HC06','HC07','HC08','HC09','HC10','HC12','HC13','HC14','HC15',    &
     'HC16','HC17','HC18','HC19','HC20','HC21','HC22','HC23','HC24',    &
     'HC25','HC26','HC27','HC28','HC29','HC31','HC32','HC33','HC34',    &
     'HC36','HC37','HC38','HC39','HC40','HC41','HC42','HC43','HC44',    &
     'HC45','HC46','HC47','HC48','HC49','PM01','PM02','PM03','PM04',    &
     'PM05','PM10-PRI'/


     DATA NAMINF /'VOC       ','NOX       ','CO        ','SO2       ',  &
       'PM10-PRI  ','PM25-PRI  ','NH3       '/

     DATA NAMVOC /   &
       ' METHANE  ',' ALKANE1  ',' ALKANE2  ',' ALKANE3  ',' ALKANE4  ',   &
       ' ALKANE5  ',' ETHYLENE ',' OLEFIN1  ',' OLEFIN2  ',' ISOPRENE ',   &
       ' TERPENES ','AROMATIC1 ','AROMATIC2 ','   CH2O   ',' CH3CHO   ',   &
       'HI_ALDEHY ','BENZALDHY ',' ACETONE  ','   MEK    ','  PRD2    ',   &
       '   MEOH   ',' GLYOXAL  ','METHGLYOX ','   BACL   ',' PHENOLS  ',   &
       ' CRESOLS  ','  MACR    ','  MVK     ','   IPRD   ',' HCOOH    ',   &
       ' CH3COOH  ',' RCOOH    ',' XYLENEOLS',' CATECHOLS',' NONVOLATL',   &
       ' PROPYLENE',' ACETYLENE',' BENZENE  ',' BUTANES  ',' PENTANES ',   &
       ' TOLUENE  ',' m-XYLENE ',' o-XYLENE ',' p-XYLENE ',' PROPANE  ',   &
       '  DIENES  ',' STYRENES ',' ETHANOL  ','ETHE_GLYC ',' UNKNOWN  '/

     DATA NAMPM2 /'PMFINE','PSO4','PNO3','POC','PEC','PNCOM','PNH4','PAL', &
       'PCA','PFE','PH2O','PMG','PMOTHR','PK','PMN','PCL','PNA','PTI','PSI'/

      DATA POINDIR / '/met1/WRF/DATA/NEI/point/' /
      DATA AREADIR / '/met1/WRF/DATA/NEI/' /

!
!
! A total kludge placed here for the plume rise.  Give some vertical wind profile climatology for the domain.  Ignoring the buoyancy effects in
!  the plume rise.
!
!
! Elevation at grid cell top (m)
!
!     DATA ZFA/    0.,   50.,  100.,  160.,  230.,  300.,  380.,  470.,  590.,  720.,             &
!                850., 1000., 1200., 1400., 1600., 1800., 2000., 2250., 2500., 2750.              /

    ! from Pablo's settings for Rim fire work
     DATA ZFA/    0.,   47.5,  95.3,  143.,  200.,  256.,  354.,  453., 553.,  653.,             &
                806.,  1006., 1209., 1416., 1628., 1845., 2106., 2374., 2649., 2997., 3320.  /

     DATA REFWZ/0.,16.8,50.5,84.3,127.,170.,256.,343.,476.,610.,791., &
      883.,975.,1160.,1350.,1644.,1943.,2252.,2677.,3010.,3350./

!      K_ht hr1   hr2   hr3   hr4   hr5   hr6   hr7   hr8   hr9  hr10  hr11  hr12  hr13  hr14  hr15  hr16  hr17  hr18  hr19  hr20  hr21  hr22  hr23  hr24
!
!        1         2         3         4         5         6         7         8         9         0         1         2         3 
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
    DATA WSPD /         &
       4.27, 4.01, 4.16, 4.27, 4.30, 4.26, 4.20, 4.12, 4.11, 4.11, 4.15, 4.25, 4.44, 3.43, 3.61, 3.86, 4.12, 4.36, 4.62, 4.94, &
       5.12, 5.21, 5.24, 5.12, 5.30, 5.62, 5.79, 5.96, 6.02, 5.96, 5.87, 5.74, 5.73, 5.70, 5.75, 5.85, 5.96, 4.51, 4.67, 4.93, &
       5.26, 5.58, 5.95, 6.40, 6.70, 6.88, 6.99, 6.91, 6.14, 6.82, 7.07, 7.27, 7.36, 7.28, 7.19, 7.02, 6.99, 6.92, 6.93, 7.02, &
       7.08, 5.46, 5.43, 5.64, 5.97, 6.31, 6.72, 7.24, 7.60, 7.88, 8.09, 8.12, 6.87, 7.73, 8.08, 8.36, 8.45, 8.37, 8.28, 8.09, &
       8.03, 7.92, 7.87, 7.91, 7.91, 6.34, 6.15, 6.22, 6.50, 6.82, 7.25, 7.80, 8.22, 8.58, 8.88, 9.04, 7.53, 8.46, 8.87, 9.22, &
       9.34, 9.26, 9.16, 8.95, 8.86, 8.70, 8.60, 8.57, 8.51, 6.94, 6.75, 6.74, 6.94, 7.21, 7.61, 8.16, 8.62, 9.05, 9.43, 9.73, &
       8.17, 9.20, 9.72,10.14,10.30,10.26,10.16, 9.92, 9.77, 9.55, 9.36, 9.25, 9.11, 7.42, 7.28, 7.38, 7.52, 7.65, 7.95, 8.44, &
       8.92, 9.42, 9.92,10.38, 8.55, 9.65,10.28,10.81,11.03,11.02,10.92,10.65,10.44,10.17, 9.93, 9.77, 9.59, 7.72, 7.72, 7.95, &
       8.29, 8.37, 8.35, 8.64, 9.07, 9.57,10.12,10.71, 8.74, 9.82,10.49,11.08,11.36,11.39,11.32,11.05,10.83,10.54,10.25,10.12, &
      10.02, 7.89, 7.96, 8.21, 8.69, 8.89, 8.78, 8.88, 9.24, 9.71,10.26,10.89, 8.82, 9.84,10.50,11.11,11.41,11.43,11.33,11.02, &
      10.75,10.44,10.18,10.14,10.25, 7.92, 8.09, 8.39, 8.93, 9.32, 9.34, 9.32, 9.59, 9.96,10.45,11.06, 8.86, 9.81,10.42,10.98, &
      11.20,11.11,10.92,10.54,10.23, 9.91, 9.71, 9.77,10.02, 7.95, 8.20, 8.62, 9.23, 9.72, 9.93, 9.92,10.03,10.31,10.70,11.19, &
       8.89, 9.74,10.28,10.76,10.88,10.72,10.46,10.07, 9.80, 9.55, 9.39, 9.50, 9.81, 8.09, 8.36, 8.83, 9.48,10.02,10.33,10.34, &
      10.40,10.60,10.90,11.26, 8.90, 9.67,10.14,10.54,10.61,10.42,10.16, 9.79, 9.56, 9.34, 9.24, 9.40, 9.73, 8.19, 8.45, 8.94, &
       9.63,10.21,10.56,10.58,10.65,10.79,11.02,11.31, 8.96, 9.56, 9.91,10.21,10.18, 9.96, 9.71, 9.39, 9.25, 9.12, 9.12, 9.33, &
       9.69, 8.38, 8.63, 9.17, 9.91,10.54,10.91,11.01,11.03,11.12,11.22,11.40, 9.10, 9.50, 9.66, 9.84, 9.77, 9.58, 9.42, 9.23, &
       9.27, 9.31, 9.36, 9.53, 9.83, 8.71, 8.98, 9.53,10.30,10.93,11.33,11.47,11.45,11.48,11.51,11.52, 9.45, 9.64, 9.60, 9.62, &
       9.54, 9.44, 9.42, 9.39, 9.58, 9.74, 9.79, 9.89,10.08, 9.24, 9.49, 9.96,10.70,11.36,11.77,11.96,11.95,11.90,11.83,11.68, &
      10.05,10.08, 9.90, 9.80, 9.72, 9.74, 9.90,10.02,10.28,10.45,10.47,10.52,10.63,10.05,10.23,10.54,11.15,11.76,12.14,12.28, &
      12.32,12.30,12.21,11.97,10.66,10.66,10.50,10.37,10.32,10.42,10.66,10.81,11.05,11.22,11.25,11.29,11.35,10.99,11.08,11.20, &
      11.62,12.10,12.37,12.44,12.52,12.62,12.62,12.43,11.40,11.53,11.45,11.42,11.43,11.50,11.68,11.77,11.98,12.21,12.33,12.40, &
      12.46,12.17,12.18,12.08,12.23,12.44,12.54,12.61,12.84,13.13,13.33,13.34,12.21,12.47,12.43,12.47,12.51,12.53,12.63,12.70, &
      12.95,13.26,13.44,13.56,13.64,13.37,13.36,13.07,12.91,12.82,12.71,12.74,13.07,13.52,13.93,14.19,12.76,13.18,13.21,13.32, &
      13.40,13.43,13.54,13.64,13.97,14.32,14.53,14.68,14.72,14.37,14.37,14.05,13.78,13.56,13.31,13.17,13.42,13.87,14.33,14.73/
!
!  Conversion table.  Table lists the field, the emissions name, the weight factor to apply, the molecular weight (not used) and the emissions
!          name  - some with units.  The fields are basically SAPRAC emissions that are being converted to RADM2 emissions input.
!
!
! Speciation profile for AL emissions into RADM2, SAM and GF 6/17/04
! CO        e_co                      1.00            28         1
! NOX       e_no                      1.00            46         1
! SO2       e_so2                     1.00            64         1
! NH3       e_nh3                     1.00            17         1
! HC01      e_ch4                     1.00            00         Methane
! HC02      e_eth                     1.00            00         Ethane  kOH<500 /ppm/min
! HC03      e_hc3                     1.00            00         Alkane 500<kOH<2500 exclude(C3H8,C2H2,ethanol,acids)
! HC04      e_hc3                     1.11            00         Alkane 2500<kOH<5000 exlude(butanes)
! HC05      e_hc5                     0.97            00         Alkane 5000<kOH<10000 exlude(pentanes)
! HC06      e_hc8                     1.00            00         Alkane kOH>10000 exclude(ethylene glycol)
! HC07      e_ol2                     1.00            00         Ethylene
! HC08      e_olt                     1.00            00         Alkene kOH <20000 /ppm/min
! HC09      e_oli                     1.00            00         Alkene kOH >20000 /ppm/min exclude(dienes,styrenes)
! HC10      e_iso                     1.00            00         Isoprene
! HC12      e_tol                     1.00            00         Aromatic kOH <20000 /ppm/min exclude(benzene and toluene)
! HC13      e_xyl                     1.00            00         Aromatic kOH >20000 /ppm/min exclude(xylenes)
! HC14      e_hcho                    1.00            00         Formaldehyde
! HC15      e_ald                     1.00            00         Acetaldehyde
! HC16      e_ald                     1.00            00         Higher aldehydes
! HC17      e_ald                     1.00            00         Benzaldehyde
! HC18      e_ket                     0.33            00         Acetone
! HC19      e_ket                     1.61            00         Methylethyl ketone
! HC20      e_ket                     1.61            00         PRD2 SAPRAC species (aromatic ketones)
! HC21      e_hc3                     0.40            00         Methanol
! HC22      e_ald                     1.00            00         Glyoxal
! HC23      e_ald                     1.00            00         Methylglyoxal
! HC24      e_ald                     1.00            00         Biacetyl
! HC25      e_csl                     1.00            00         Phenols
! HC26      e_csl                     1.00            00         Cresols
! HC27      e_ald                     0.50            00         Methacrolein
! HC27      e_olt                     0.50            00         Methacrolein
! HC28      e_ket                     0.50            00         Methylvinyl ketone
! HC28      e_olt                     0.50            00         Methylvinyl ketone
! HC29      e_ket                     1.00            00         IPRD SAPRAC species (>C4 unsaturated aldehydes)
! HC31      e_ora2                    1.00            00         Acetic Acid
! HC32      e_ora2                    1.00            00         >C2 Acids  (SAPRC PACD species)
! HC33      e_csl                     1.00            00         Xylenols  (SAPRC-11 species)
! HC34      e_csl                     1.00            00         Catechols  (SAPRC-11 species)
! HC36      e_olt                     1.00            00         Propylene
! HC37      e_hc3                     0.40            00         Acetylene
! HC38      e_tol                     0.29            00         Benzene
! HC39      e_hc3                     1.11            00         Butanes
! HC40      e_hc5                     0.97            00         Pentanes
! HC41      e_tol                     1.00            00         Toluene
! HC42      e_xyl                     1.00            00         m-Xylene
! HC43      e_xyl                     1.00            00         p-Xylene
! HC44      e_xyl                     1.00            00         o-Xylene
! HC45      e_hc3                     0.57            00         Propane
! HC46      e_oli                     1.00            00         Dienes
! HC47      e_olt                     1.00            00         Styrenes
! HC47      e_tol                     1.00            00         Styrenes
! HC48      e_hc3                     1.20            00         Ethanol
! HC49      e_hc8                     1.14            00         Ethylene Glycol
! PM01      e_pm25i                   0.20            01         Unspeciated primary PM2.5 - nuclei mode
! PM01      e_pm25j                   0.80            01         Unspeciated primary PM2.5 - accumulation mode
! PM02      e_so4i                    0.20            01         Sulfate PM2.5 - nuclei mode
! PM02      e_so4j                    0.80            01         Sulfate PM2.5 - accumulation mode
! PM03      e_no3i                    0.20            01         Nitrate PM2.5 - nuclei mode
! PM03      e_no3j                    0.80            01         Nitrate PM2.5 - accumulation mode
! PM04      e_orgi                    0.20            01         Organic Carbon PM2.5 - nuclei mode
! PM04      e_orgj                    0.80            01         Organic Carbon PM2.5 - accumulation mode
! PM05      e_eci                     0.20            01         Elemental Carbon PM2.5 - nuclei mode
! PM05      e_ecj                     0.80            01         Elemental Carbon PM2.5 - accumulation mode
! PM10-PRI  e_pm10                    1.00            01         Unspeciated Primary PM10

     DATA NAM /  &
       'CO       ','NOX      ','SO2      ','NH3      ','HC01     ','HC02     ','HC03     ','HC04     ','HC05     ',  &
       'HC06     ','HC07     ','HC08     ','HC09     ','HC10     ','HC12     ','HC13     ','HC14     ',  &
       'HC15     ','HC16     ','HC17     ','HC18     ','HC19     ','HC20     ','HC21     ','HC22     ',  &
       'HC23     ','HC24     ','HC25     ','HC26     ','HC27     ','HC27     ','HC28     ','HC28     ',  &
       'HC29     ','HC31     ','HC32     ','HC33     ','HC34     ','HC36     ','HC37     ','HC38     ',  &
       'HC39     ','HC40     ','HC41     ','HC42     ','HC43     ','HC44     ','HC45     ','HC46     ',  &
       'HC47     ','HC47     ','HC48     ','HC49     ','PM01     ','PM01     ','PM02     ','PM02     ',  &
       'PM03     ','PM03     ','PM04     ','PM04     ','PM05     ','PM05     ','PM10-PRI '/

     DATA NAMRAD /  &
       'e_co     ','e_no     ','e_so2    ','e_nh3    ','e_ch4    ','e_eth    ','e_hc3    ','e_hc3    ','e_hc5    ',  &
       'e_hc8    ','e_ol2    ','e_olt    ','e_oli    ','e_iso    ','e_tol    ','e_xyl    ','e_hcho   ',  &
       'e_ald    ','e_ald    ','e_ald    ','e_ket    ','e_ket    ','e_ket    ','e_hc3    ','e_ald    ',  &
       'e_ald    ','e_ald    ','e_csl    ','e_csl    ','e_ald    ','e_olt    ','e_ket    ','e_olt    ',  &
       'e_ket    ','e_ora2   ','e_ora2   ','e_csl    ','e_csl    ','e_olt    ','e_hc3    ','e_tol    ',  &
       'e_hc3    ','e_hc5    ','e_tol    ','e_xyl    ','e_xyl    ','e_xyl    ','e_hc3    ','e_oli    ',  &
       'e_olt    ','e_tol    ','e_hc3    ','e_hc8    ','e_pm25i  ','e_pm25j  ','e_so4i   ','e_so4j   ',  &
       'e_no3i   ','e_no3j   ','e_orgi   ','e_orgj   ','e_eci    ','e_ecj    ','e_pm10   '/

     DATA FAC /  &
              1.00,       1.00,       1.00,       1.00,       1.00,       1.00,       1.00,       1.11,       0.97,  &
              1.00,       1.00,       1.00,       1.00,       1.00,       1.00,       1.00,       1.00,  &
              1.00,       1.00,       1.00,       0.33,       1.61,       1.61,       0.40,       1.00,  &
              1.00,       1.00,       1.00,       1.00,       0.50,       0.50,       0.50,       0.50,  &
              1.00,       1.00,       1.00,       1.00,       1.00,       1.00,       0.40,       0.29,  &
              1.11,       0.97,       1.00,       1.00,       1.00,       1.00,       0.57,       1.00,  &
              1.00,       1.00,       1.20,       1.14,       0.20,       0.80,       0.20,       0.80,  &
              0.20,       0.80,       0.20,       0.80,       0.20,       0.80,       1.00/

     DATA MWT /  &
                28,         46,         64,         17,         00,         00,         00,         00,         00,  &
                00,         00,         00,         00,         00,         00,         00,         00,  &
                00,         00,         00,         00,         00,         00,         00,         00,  &
                00,         00,         00,         00,         00,         00,         00,         00,  &
                00,         00,         00,         00,         00,         00,         00,         00,  &
                00,         00,         00,         00,         00,         00,         00,         00,  &
                00,         00,         00,         00,         01,         01,         01,         01,  &
                00,         00,         01,         01,         01,         01,         01/


     CONTAINS

     SUBROUTINE  AL2RADM2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     REAL, DIMENSION(IX,JX)         :: XLAT,XLON
     REAL, DIMENSION(NAL2DO,NRADM)  ::  TRANAL2R

     INTEGER                        :: STATUS,SYSTEM

     CHARACTER (len=128)            :: HRDIR
!
! Variables for Point stack emissions
!
     REAL, DIMENSION(IPOINT)        ::  STKHGT, STKDIAM, STKTMP
     REAL, DIMENSION(IPOINT)        ::  STKVEL, STKFLOW, XLNP, XLTP
     INTEGER, DIMENSION(IPOINT)     :: ISTATE, ICOUN, IRTYP

     CHARACTER (len=15), DIMENSION(IPOINT) ::  SITEID
     CHARACTER (len= 8), DIMENSION(IPOINT) ::  ERPTID, UNITID, PROCID

     INTEGER     ::  IAREADO = 1
     INTEGER     ::  IPOINDO = 1
     INTEGER     ::  IFRACGR = 0        ! set to 1 if dx < 10 km

     INTEGER     ::  KK, IHR
     INTEGER     ::  I2, J2
     INTEGER     ::  IEMAX, JEMAX, IMAX
     INTEGER     ::  IDHMX, KSTAK, KSTKDHM, KBOT, KTOP

     REAL        ::  XX, YY
     REAL        ::  EMAX, XLTMX, XLNMX
     REAL        ::  DHM, DHMX, TOP, BOT, ZDIF, FRC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     IF ( dx .LE. 10000. .and. IFRACGR .EQ. 0 ) then
!       WRITE(*,*) ' ERROR: grid dumping no recommended for grid spacing under 10 km '
!       WRITE(*,*) ' Set IFRACGR = 1 '
!       STOP 999
!     ENDIF

! Open output files, 7 = print and log, 19 = unformatted wrf output
     OPEN(19,FILE='wrfem_00to12Z',FORM='UNFORMATTED')
     WRITE(19)NRADM
     WRITE(19)ename
     LDEV=17
     OPEN(LDEV,FILE='al2radm2.outp')
     L1P=LNBLNK(POINDIR)
     WRITE(LDEV,*)'Geometric grid fraction flag = ',IFRACGR
!
! Read point realease info file, print information for maximum stack height in inventory
!
         CALL MAPCF(1.5,1.5,XLTMX,XLNMX)
         WRITE(LDEV,*)' Center of 1,1 grid = ',XLTMX,XLNMX
     IF(IPOINDO.EQ.1)THEN
     WRITE(LDEV,'(/30A1,A,30A1)')('-',I=1,30),' Release point info file ',('-',I=1,30)
     OPEN(8,FILE=POINDIR(1:L1P)//'Relpnt_info.txt')
     SHGTMX=-1.E20
     DO I=1,IPOINT
        READ(8,55)ISTATE(I),ICOUN(I),SITEID(I),ERPTID(I),UNITID(I), &
        PROCID(I),IRTYP(I),STKHGT(I),STKDIAM(I),STKTMP(I),STKVEL(I), &
        STKFLOW(I),XLNP(I),XLTP(I)

        IF(STKHGT(I).GT.SHGTMX)THEN
          IMAX=I
          SHGTMX=STKHGT(I)
        ENDIF
     ENDDO
 55   FORMAT(I2.2,I3.3,A15,3A8,1X,I2,1X,F10.3,F10.3,F10.2,F10.3,F10.4,F11.5,F10.5)
     CLOSE(8)

     WRITE(LDEV,*)'Max stack ht, Record number of max= ',IMAX
     WRITE(LDEV,*)'Max stack ht, FIPS State= ',ISTATE(IMAX)
     WRITE(LDEV,*)'Max stack ht, FIPS County= ',ICOUN(IMAX)
     WRITE(LDEV,*)'Max stack ht, Site ID= ',SITEID(IMAX)
     WRITE(LDEV,*)'Max stack ht, Report ID= ',ERPTID(IMAX)
     WRITE(LDEV,*)'Max stack ht, Unit ID= ',UNITID(IMAX)
     WRITE(LDEV,*)'Max stack ht, Process ID= ',PROCID(IMAX)
     WRITE(LDEV,*)'Max stack ht, Release Type= ',IRTYP(IMAX)
     WRITE(LDEV,*)'Max stack ht, Height= ',STKHGT(IMAX),' meter'
     WRITE(LDEV,*)'Max stack ht, Stack Diameter= ',STKDIAM(IMAX),' meter'
     WRITE(LDEV,*)'Max stack ht, Exit Temperature= ',STKTMP(IMAX),' deg K'
     WRITE(LDEV,*)'Max stack ht, Flow velocity= ',STKVEL(IMAX),' m/s'
     WRITE(LDEV,*)'Max stack ht, Flow rate= ',STKFLOW(IMAX),' m(3)/s'
     WRITE(LDEV,*)'Max stack ht, Longitude= ',XLNP(IMAX)
     WRITE(LDEV,*)'Max stack ht, Latitude= ',XLTP(IMAX)
     ENDIF ! Endif test on IPOINDO=1
!
! Read in table for mapping AL emissions into RADM2 emissions
!
     TRANAL2R=0.
     DO KK=1,NAMFILE2
     DO K=1,NAL2DO
       IF( NAM(KK) .EQ. NAM2EM(K)) THEN
       DO J=1,NRADM
        IF(NAMRAD(KK) .EQ. ENAME(J))THEN
           IF(MWT(KK).EQ.0)THEN
              TRANAL2R( K, J ) = FAC(KK)/(1.E-3*DX)**2  ! Units of VOC in mole/hr --> mole/km(2)/hr
           ELSEIF(MWT(KK).EQ.1)THEN    ! Units of aerosol in ton/hr --> microgram/m(2)/sec
              TRANAL2R( K, J ) = FAC(KK)*9.07184E11 /DX/DX/3600.
           ELSE  ! Units of primary species ton/hr --> mole/km(2)/hr
              TRANAL2R( K, J ) = FAC(KK)*9.07184E5/FLOAT(MWT(KK))/(1.E-3*DX)**2
           ENDIF
        ELSE
!           WRITE(LDEV,*)NAMRAD,' RADM name Not found, stopping'
!           STOP888
        ENDIF
       ENDDO
!       ELSE
!         WRITE(LDEV,*)NAM,' Not found in table emiss. list, stopping'
!         STOP888
       ENDIF
     ENDDO
     ENDDO

     IF(IAREADO .NE. 0) THEN
!
! Read Cross point and Dot point Latitude and Longitudes, print corners
!
        L1A=LNBLNK(AREADIR)
        WRITE(LDEV,'(/30A1,A,30A1)')('-',I=1,30),' Lat, Lon position files ',('-',I=1,30)
        OPEN(8,FILE=AREADIR(1:L1A)//'grid_loc/LAT_dot.txt')
        READ(8,'(12F10.5)')XLATD
        CLOSE(8)
        OPEN(8,FILE=AREADIR(1:L1A)//'grid_loc/LON_dot.txt')
        READ(8,'(12F10.5)')XLOND
        CLOSE(8)
        WRITE(LDEV,'(A,2F11.5)')'(1,1) Dot Point epa4km Lat,Lon= ',XLATD(1,1),XLOND(1,1)
        WRITE(LDEV,'(2(A,I4),A,2F11.5)')'(',IP,',',JP,') Dot Point epa4km Lat,Lon= ',XLATD(IP,JP),XLOND(IP,JP)
        OPEN(8,FILE=AREADIR(1:L1A)//'grid_loc/LAT_xrs.txt')
        READ(8,'(12F10.5)')XLAT
        CLOSE(8)
        OPEN(8,FILE=AREADIR(1:L1A)//'grid_loc/LON_xrs.txt')
        READ(8,'(12F10.5)')XLON
        CLOSE(8)
        WRITE(LDEV,'(A,2F11.5)')'(1,1) Cross Point epa4km Lat,Lon= ',XLAT(1,1),XLON(1,1)
        WRITE(LDEV,'(2(A,I4),A,2F11.5)')'(',IX,',',JX,') Cross Point epa4km Lat,Lon= ',XLAT(IX,JX),XLON(IX,JX)
!
! Read Daily Ozone Season Day Average (short ton/dy) primary emissions,
! print sums, maximums, location of maximums in grid space and lat,lon
!
        WRITE(LDEV,'(/30A1,A,30A1)')('-',I=1,30),' Daily emis. ton/dy  ',('-',I=1,30)
        WRITE(LDEV,'(A,4X,A,7X,A,6X,A,2X,A,2X,A,4X,A)')'Species','Total','Max','I-Max','J-Max','Lat-Max','Lon-Max'

        DO ISP=1,IPRIM
          EM2D=0.
          ETOT = 0.
          EMAX = -1.E20
! Unzip area emission file
          STATUS=SYSTEM('rm -f scratem*')
          CHSPEC=NAMINF(ISP)
          L1D=LNBLNK(CHSPEC)
          STATUS=SYSTEM('cp '//AREADIR(1:L1A)//'area4k/dayav/'//CHSPEC(1:L1D)//'.gz scratem.gz')
          STATUS=SYSTEM('gunzip scratem')
          OPEN(8,FILE='scratem',FORM='FORMATTED')
          READ(8,'(12E9.2)')EMT
          CLOSE(8)
     IF(IFRACGR.EQ.1)THEN   ! IFRACGR=1 means no griddumping of area emissions - call REGRID
          WRITE(LDEV,*)'IFRACGR=1 not supported in emiss_v03 public release'
          STOP'777'
     ELSE  ! IFRACGR not= 1 means griddumping of area emissions
! apportion 4km grid emissions into wrf-27km grids by simple grid dumping
          DO I=1,IX
          DO J=1,JX
            CALL LLXY(XLAT(I,J),XLON(I,J),XX,YY)
!            WRITE(LDEV,*)'i,j,lat,lon,x,y= ',I,J,XLAT(I,J),XLON(I,J),XX,YY
            IF(XX .GE. 1 .AND. XX.LE.IX2 .AND. YY.GE.1 .AND. YY.LE.JX2) THEN
              I2=INT(XX)
              J2=INT(YY)
              EM2D(I2,J2) = EM2D(I2,J2)+EMT(I,J)
            ENDIF
          ENDDO
          ENDDO
      ENDIF   ! end of IFRACGR=1 or not test
!          Get total emission and max emission information
          CALL EMSUM(ETOT,EMAX,IEMAX,JEMAX)

!          Get approximate lat,lon of center of grid point
          CALL MAPCF(FLOAT(IEMAX)+.5,FLOAT(JEMAX)+.5,XLTMX,XLNMX)

!          Output Total and Max info
          WRITE(LDEV,'(A9,1P2E11.3,2I6,0P2F11.5)')NAMINF(ISP),ETOT,EMAX,IEMAX,JEMAX,XLTMX,XLNMX

        ENDDO   ! End of ISP loop (primary species)
     ENDIF

    endhr = min(24,starthr+ maxhr)
    DO IHR=starthr,maxhr                ! DO 172 LOOP
! Loop over all 24 hours
!   DO IHR=1,24                   ! DO 172 LOOP      
       EM3RD=0.

       IF(IAREADO .NE. 0) THEN
          WRITE(HRDIR,'(A,A,I2.2,A)')AREADIR(1:L1A),'area4k/HR',IHR,'/'
          WRITE(LDEV,*)'HRDIR= ',HRDIR
          L1E=LNBLNK(HRDIR)
!
! Read UTC Hour 20, Ozone Season Day Average (mole/hr) speciated VOC emissions,
! print sums, maximums, location of maximums in grid space and lat,lon
!
         WRITE(LDEV,'(/30A1,A,I2.2,A,30A1)')('-',I=1,30),' Hr_',IHR, &
                  ' Area em mole/hr or ton/hr (prim+asol) ',('-',I=1,30)
         WRITE(LDEV,'(A,4X,A,7X,A,6X,A,2X,A,2X,A,4X,A)')'ALemisn','Total', &
                  'Max','I-Max','J-Max','Lat-Max','Lon-Max'

         DO ISP=1,NAL2DO
           ETOT = 0.

           EM2D=0.
! Unzip area emission file
           STATUS=SYSTEM('rm -f scratem*')
           CHSPEC=NAM2EM(ISP)
           L1D=LNBLNK(CHSPEC)

           STATUS=SYSTEM('cp '//HRDIR(1:L1E)//CHSPEC(1:L1D)//'.gz scratem.gz')
           STATUS=SYSTEM('gunzip scratem')
           OPEN(8,FILE='scratem',FORM='FORMATTED')
           READ(8,'(12E9.2)')EMT

           CLOSE(8)
     IF(IFRACGR.EQ.1)THEN   ! IFRACGR=1 means no griddumping of area emissions - call REGRID
          WRITE(LDEV,*)'IFRACGR=1 not supported in emiss_v03 public release'
          STOP'777'
     ELSE  ! IFRACGR not= 1 means griddumping of area emissions
! apportion 4km grid emissions into wrf-27km grids by simple grid dumping
           DO I=1,IX
           DO J=1,JX
             CALL LLXY(XLAT(I,J),XLON(I,J),XX,YY)
             IF(XX .GE. 1 .AND. XX.LE.IX2 .AND. YY.GE.1 .AND. YY.LE.JX2) THEN
               I2=INT(XX)
               J2=INT(YY)
               EM2D(I2,J2)=EM2D(I2,J2)+EMT(I,J)
             ENDIF
           ENDDO
           ENDDO
      ENDIF   ! end of IFRACGR=1 or not test
!       Get total emission and max emission information
           EMAX = -1.E20
           CALL EMSUM(ETOT,EMAX,IEMAX,JEMAX)

!       Get approximate lat,lon of center of grid point
           CALL MAPCF(FLOAT(IEMAX)+.5,FLOAT(JEMAX)+.5,XLTMX,XLNMX)

!       Get real name, depending on primary, VOC or PM2.5
           CHSCRT = NAM2EM(ISP)
           IF(CHSCRT(1:2).EQ.'HC')THEN
             READ(CHSCRT(3:4),'(I2)')IDEX
             CHSPEC=NAMVOC(IDEX)
           ELSEIF(CHSCRT(1:3).EQ.'PM0')THEN
             READ(CHSCRT(3:4),'(I2)')IDEX
             CHSPEC=NAMPM2(IDEX)
           ELSE
             CHSPEC=CHSCRT
           ENDIF
! Output Total and Max info
           WRITE(LDEV,'(A9,1P2E11.3,2I6,0P2F11.5)')CHSPEC,ETOT,EMAX,IEMAX,JEMAX,XLTMX,XLNMX

! Partition Emissions of this species into wrf-RADM2 emission
           DO N=1,NRADM
           DO J=1,JX2
           DO I=1,IX2
             EM3RD(I,1,J,N)=EM3RD(I,1,J,N)+EM2D(I,J)*TRANAL2R(ISP,N)
           ENDDO
           ENDDO
           ENDDO
         ENDDO   ! End of ISP loop

! Get total emission and max emission information for RADM emissions
         WRITE(LDEV,'(/30A1,A,I2.2,A,30A1)')('-',I=1,30),' Hr_',IHR, &
                  ' RADM2 Area em mole/km2/hr or ug/m2/sec (asol)',('-',I=1,30)
         WRITE(LDEV,'(A,4X,A,7X,A,6X,3(A,2X),A,4X,A)')'RADM em','Total', &
                  'Max','I-Max','J-Max','K-Max','Lat-Max','Lon-Max'

         DO N=1,NRADM
         ETOT = 0.
            EMAX = -1.E20
            CALL EMSUM3D(ETOT,EMAX,IEMAX,JEMAX,KEMAX,N)
            CALL MAPCF(FLOAT(IEMAX)+.5,FLOAT(JEMAX)+.5,XLTMX,XLNMX)
            WRITE(LDEV,'(A9,1P2E11.3,3I6,0P2F11.5)')ename(N),ETOT,EMAX,IEMAX,JEMAX,KEMAX,XLTMX,XLNMX
         ENDDO
       ENDIF  ! END IAREADO IF BLOCK

     IF(IPOINDO.EQ.1)THEN
! Begin point stack emissions processing for this hour
       WRITE(HRDIR,'(A,A,I2.2,A)')POINDIR(1:L1P),'HR',IHR,'/'
       WRITE(LDEV,*)'HRDIR= ',HRDIR
       L1E=LNBLNK(HRDIR)
!
! Read UTC Hour 20, Ozone Season Day Average (mole/hr) speciated VOC emissions,
! print sums, maximums, location of maximums in grid space and lat,lon
!
       WRITE(LDEV,'(/30A1,A,I2.2,A,30A1)')('-',I=1,30),' Hr_',IHR,  &
               ' Point em mole/hr or ton/hr (prim+asol) ',('-',I=1,30)
       WRITE(LDEV,'(A,4X,A,7X,A,6X,A,2X,A,2X,A,4X,A)')'ALpoint','Total',   &
               'Max','I-Max','J-Max','Lat-Max','Lon-Max'

       DO ISP=1,NAL2DO
       ETOT = 0.
! Unzip point emission file
         STATUS=SYSTEM('rm -f scratem*')
         CHSPEC=NAM2EM(ISP)
         L1D=LNBLNK(CHSPEC)

         STATUS=SYSTEM('cp '//HRDIR(1:L1E)//CHSPEC(1:L1D)//'.gz scratem.gz')
         STATUS=SYSTEM('gunzip scratem')
         OPEN(8,FILE='scratem',FORM='FORMATTED')
         READ(8,'(12E9.2)')EMP
         CLOSE(8)

!         Get total emission and max emission information
         EMAX = -1.E20
         CALL EMSUMP(ETOT,EMAX,IEMAX)

!         Output Total and Max info
         CHSCRT=NAM2EM(ISP)
         IF(CHSCRT(1:2).EQ.'HC')THEN
           READ(CHSCRT(3:4),'(I2)')IDEX
           CHSPEC=NAMVOC(IDEX)
         ELSEIF(CHSCRT(1:3).EQ.'PM0')THEN
           READ(CHSCRT(3:4),'(I2)')IDEX
           CHSPEC=NAMPM2(IDEX)
         ELSE
           CHSPEC=CHSCRT
         ENDIF
         WRITE(LDEV,'(A9,1P2E11.3,I7,0P2F11.5)')CHSPEC,ETOT,EMAX,IEMAX,XLTP(IEMAX),XLNP(IEMAX)

! apportion 4km grid emissions into wrf-27km grids by simple grid dumping

        DHMX = 0.
        DO I=1,IPOINT
          CALL LLXY(XLTP(I),XLNP(I),XX,YY)
          IF(XX .GE. 1 .AND. XX.LE.IX2 .AND. YY.GE.1 .AND. YY.LE.JX2) THEN
            I2=INT(XX)
            J2=INT(YY)
! Find K level of emissions
            K=1
            ZTOP=REFWZ(K)
            DO WHILE (ZTOP.LE.STKHGT(I))
            K=K+1
             IF(K.GT.KWINP)THEN
             WRITE(LDEV,*)'Stack hght > Wind Climatology Ht, stopping'
             STOP'7777'
             ENDIF
            ZTOP=REFWZ(K)   
            ENDDO
            KSTAK=MAX(1,K-1)
            DHM   = 3.0 * STKDIAM(I) * STKVEL(I) / WSPD(IHR,KSTAK)
            IF(DHM.GT.DHMX)THEN
               DHMX=DHM
               IDHMX=I
               KSTKDHM=KSTAK
            ENDIF
            TOP = STKHGT(I) + 1.5 * DHM
            BOT = STKHGT(I) + 0.5 * DHM

            DO K=KX,1,-1
               IF (ZFA(K+1).GT.BOT) THEN
                 KBOT=K
               ENDIF
            ENDDO

            DO K=KX,1,-1
              IF(ZFA(K+1) .GT. TOP)THEN
                KTOP=K
              ENDIF
            ENDDO

            IF( KBOT .GE. KTOP) THEN
              KTOP = KBOT + 1
            ENDIF

           ZDIF = ZFA(KTOP+1) - ZFA(KBOT)
! Partition Emissions of this species into wrf-RADM2 emission
           DO K=KBOT,KTOP
             FRC=(ZFA(K+1)-ZFA(K))/ZDIF
             DO N=1,NRADM
               EM3RD(I2,K,J2,N)=EM3RD(I2,K,J2,N)+EMP(I)*TRANAL2R(ISP,N)*FRC
             ENDDO
           ENDDO
          ENDIF  ! Endo f 1<=I<=IX2,1<=J<=JX2 test
        ENDDO  ! I=1,IPOINT loop
       ENDDO   ! End of ISP loop
     ENDIF ! Endif test on IPOINDO=1

       WRITE(LDEV,'(/30A1,A,I2.2,A,30A1)')('-',I=1,30),' Hr_',IHR, &
              ' RADM2 Point+Area em mole/km2/hr or ug/m2/sec (asol)',('-',I=1,30)
       WRITE(LDEV,'(A,4X,A,7X,A,6X,3(A,2X),A,4X,A)')'RADM em','Total','Max', &
             'I-Max','J-Max','K-Max','Lat-Max','Lon-Max'

! File size too big - split emissions files in half; 0-12Z , 12-24Z

       IF(IHR.EQ.13)THEN
         CLOSE(19)
         OPEN(19,FILE='wrfem_12to24Z',FORM='UNFORMATTED')
         WRITE(19)NRADM
         WRITE(19)ename
       ENDIF

       WRITE(19)IHR

       DO N=1,NRADM
       ETOT = 0.
         EMAX = -1.E20
         CALL EMSUM3D(ETOT,EMAX,IEMAX,JEMAX,KEMAX,N)
         CALL MAPCF(FLOAT(IEMAX)+.5,FLOAT(JEMAX)+.5,XLTMX,XLNMX)
         WRITE(LDEV,'(A9,1P2E11.3,3I6,0P2F11.5)')ename(N),ETOT,EMAX,IEMAX,JEMAX,KEMAX,XLTMX,XLNMX

!         Write out 3-D emission arrays to unformatted file

         DO I=1,IX2
         DO K=1,KX
         DO J=1,JX2
           EM3RS(I,K,J)=EM3RD(I,K,J,N)
         ENDDO
         ENDDO
         ENDDO
         WRITE(19)EM3RS
       ENDDO  ! end of N=1,NRADM loop for writes

     IF(IPOINDO.EQ.1)THEN
       I=IDHMX
       WRITE(LDEV,*)'DHMX= ',DHMX,' IDHMX= ',I,' KSTK= ',KSTKDHM
       WRITE(LDEV,55)ISTATE(I),ICOUN(I),SITEID(I),ERPTID(I),UNITID(I), PROCID(I), &
                  IRTYP(I),STKHGT(I),STKDIAM(I),STKTMP(I),STKVEL(I), &
                  STKFLOW(I),XLNP(I),XLTP(I)
     ENDIF ! Endif test on IPOINDO=1
     CALL FLUSH(LDEV)

     END DO      ! END IHR=1,24 LOOP, OR DO 172 LOOP

     CLOSE(LDEV)
     CLOSE(19)
     END SUBROUTINE  AL2RADM2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     SUBROUTINE EMSUM(ETOT,EMAX,IEMAX,JEMAX)
!     Calculate total and max
!     EMAX=-1.E25
!     ETOT=0.
     DO I=1,IX2
     DO J=1,JX2
       ETOT=ETOT+EM2D(I,J)
       IF(EM2D(I,J).GT.EMAX)THEN
         EMAX=EM2D(I,J)
         IEMAX=I
         JEMAX=J
       ENDIF
     ENDDO
     ENDDO

     RETURN
     END SUBROUTINE  EMSUM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     SUBROUTINE EMSUM3D(ETOT,EMAX,IEMAX,JEMAX,KEMAX,N)
!     Calculate total and max
!     EMAX=-1.E25
!     ETOT=0.
     DO I=1,IX2
     DO J=1,JX2
     DO K=1,KX
       ETOT=ETOT+EM3RD(I,K,J,N)
       IF(EM3RD(I,K,J,N).GT.EMAX)THEN
         EMAX=EM3RD(I,K,J,N)
         IEMAX=I
         JEMAX=J
         KEMAX=K
       ENDIF
     ENDDO
     ENDDO
     ENDDO

     RETURN
     END SUBROUTINE  EMSUM3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     SUBROUTINE EMSUMP(ETOT,EMAX,IEMAX)
!     Calculate total and max
!     EMAX=-1.E25
!     ETOT=0.
     DO I=1,IPOINT
       ETOT=ETOT+EMP(I)
       IF(EMP(I).GT.EMAX)THEN
         EMAX=EMP(I)
         IEMAX=I
       ENDIF
     ENDDO
     RETURN
     END SUBROUTINE  EMSUMP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     SUBROUTINE MAPCF (XI,YJ,XLAT,XLON)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                     C
!                                                                     C
!     THIS SUBROUTINE COMPUTES THE LATITUDE AND LONGITUDE FROM MODEL  C
!     INDEXES OF A POINT.                                             C
!                                                                     C
!     INPUT :                                                         C
!                                                                     C
!        XI : X COORDINATE OF THE POINT IN MODEL INDEX.               C
!                                                                     C
!        YJ : Y COORDINATE OF THE POINT IN MODEL INDEX.               C
!                                                                     C
!     OUTPUT :                                                        C
!                                                                     C
!        XLAT : LATITUDE OF THE POINT.                                C
!                                                                     C
!        XLON : LONGITUDE OF THE POINT.                               C
!     NOTE ***                                                        C
!                                                                     C
!     A WEST LONGITUDE IS GIVEN BY A NEGATIVE NUMBER; POSITIVE        C
!     ANGLES DENOTE EAST LONGITUDE.                                   C
!                                                                     C
!     A NORTH LATITUDE IS GIVEN BY A POSITIVE NUMBER, AND A NEGATIVE  C
!     NUMBER FOR A SOUTH LATITUDE.                                    C
!                                                                     C
!                                                                     C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     REAL        ::  XN, PSI1, FRTY5D, CONV, RLAT1, RLAT2, POLE
     REAL        ::  PSX, CELL, CELL2, R, C2, XCNTR, PHICTR, YCNTR
     REAL        ::  CNTRJ, CNTRI, XIBD, YJBD, X, Y, FLP, FLPP, RXN
     REAL        ::  CEL1, CEL2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-----------------------------------------------------------------------
!
     XN = 0.
!  9/5/03 Modified for more general Lambert Conformal with 2 standard parallels,CLAT1,CLAT2
     PSI1=CLAT1
     FRTY5D=ATAN(1.)
     CONV=45./FRTY5D
     IF (IPROJ.EQ.1)THEN
!23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
       RLAT1=CLAT1/CONV
       RLAT2=CLAT2/CONV
        IF(CLAT1.NE.CLAT2)THEN
        XN=ALOG(COS(RLAT1)/COS(RLAT2))/ALOG(TAN(FRTY5D+.5*RLAT2)/ &
        TAN(FRTY5D+.5*RLAT1))
        ELSE
        XN=SIN(RLAT1)
        ENDIF
!     WRITE(*,*)'IPROJ,XN,PSI1,CONV,REKM,DX= ',IPROJ,XN,PSI1,CONV,REKM,DX
     ENDIF
     IF (IPROJ.EQ.2) XN = 1.0
!
!-----PSI1 IS LATITUDE WHERE CONE OR PLANE INTERSECTS EARTH
!
     POLE = 90.
     PSI1 = PSI1/CONV
     IF (XLATC.LT.0.)POLE=-90.
!
!-----REKM IS RADIUS OF EARTH IN KM
!     SUBTRACT XROT DEGREES (modify XCNTR,YCNTR) TO ROTATE LAMBERT CONFORMAL PROJECTION
!     CALCULATE R
!
     IF (IPROJ.NE.3) THEN
        PSX = (POLE-XLATC)/CONV
        IF (IPROJ.EQ.1) THEN
           IF(CLAT1.EQ.CLAT2)THEN
           CELL = REKM/tan(psi1)
           CELL2=1.
           ELSE
           CELL =REKM*SIN(PSI1)/XN
           CELL2 = (TAN(PSX/2.))/(TAN(PSI1/2.))
           ENDIF
        ENDIF
        IF (IPROJ.EQ.2) THEN
           CELL = REKM*SIN(PSX)/XN
           CELL2 = (1.+COS(PSI1))/(1.+COS(PSX))
        ENDIF
        R = CELL*(CELL2)**XN
        XCNTR = 0.0
        YCNTR = -R
     ENDIF
!
!-----FOR MERCATOR PROJECTION, THE PROJECTION IS TRUE AT LATITUDE PSI1
!
     IF (IPROJ.EQ.3) THEN
        C2 = REKM*COS(PSI1)
        XCNTR = 0.0
        PHICTR = XLATC/CONV
        CELL = (1.0+SIN(PHICTR))/COS(PHICTR)
        YCNTR = C2*ALOG(CELL)
     ENDIF
!     WRITE (6,1010) XCNTR,YCNTR
!1010 FORMAT (1X,'X COORD GRID CNTR = ',F8.1,' Y COORD = ',F8.1)
!
!-----CALCULATE X AND Y POSITIONS OF GRID
!
     CNTRJ = (JL+1)/2.
     CNTRI = (IL+1)/2.
     IF(INEST1.EQ.1)THEN
     XIBD=XNESSTR+(XI-1)*DX/DXBIGDO
     YJBD=YNESSTR+(YJ-1)*DX/DXBIGDO
     ELSE
     XIBD=XI
     YJBD=YJ
     ENDIF
     X = XCNTR+(XIBD-CNTRI)*DXBIGDO/1000.
     Y = YCNTR+(YJBD-CNTRJ)*DXBIGDO/1000.
!     WRITE (6,1020) X,Y
!1020 FORMAT (1X,'X INDEX = ',F8.1,' Y INDEX = ',F8.1)
!
!-----NOW CALC LAT AND LONG OF THIS POINT
!
     IF (IPROJ.NE.3) THEN
        IF (XLATC.LT.0.0) THEN
           FLP = ATAN2(X,Y)
        ELSE
           FLP = ATAN2(X,-Y)
        ENDIF
        FLPP = (FLP/XN)*CONV+XLONC
        IF (FLPP.LT.-180.) FLPP = FLPP+360.
        IF (FLPP.GT.180.) FLPP = FLPP-360.
        XLON = FLPP
!
!-----NOW SOLVE FOR LATITUDE
!
        R = SQRT(X*X+Y*Y)
        IF (XLATC.LT.0.0) R = -R
        IF (IPROJ.EQ.1) THEN
           IF(CLAT1.EQ.CLAT2)THEN
           XLAT = (1./tan(psi1) + psi1 - r/rekm)*conv
           ELSE
           CELL = (R*XN)/(REKM*SIN(PSI1))
           RXN = 1.0/XN
           CEL1 = TAN(PSI1/2.)*(CELL)**RXN
           CEL2 = ATAN(CEL1)
           PSX = 2.*CEL2*CONV
           XLAT = POLE-PSX
           ENDIF
        ENDIF
        IF (IPROJ.EQ.2) THEN
           CELL = R/REKM
        CEL1 = CELL/(1.0+COS(PSI1))
        CEL2 = ATAN(CEL1)
        PSX = 2.*CEL2*CONV
        XLAT = POLE-PSX
        ENDIF
     ENDIF
!
!     CALCULATIONS FOR MERCATOR LAT, LONG AND MAP SCALES...
!
     IF (IPROJ.EQ.3) THEN
        XLON = XLONC+((X-XCNTR)/C2)*CONV
        CELL = EXP(Y/C2)
        XLAT = 2.0*(CONV*ATAN(CELL))-90.0
     ENDIF
     RETURN
!
     END SUBROUTINE  MAPCF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     SUBROUTINE LLXY (XLAT,XLON,XI,YJ)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     C
!                                                                     C
!     THIS SUBROUTINE COMPUTES THE MODEL INDEXES OF A POINT FROM      C
!     LATITUDE AND LONGITUDE.                                         C
!                                                                     C
!     INPUT :                                                         C
!                                                                     C
!        XLAT : LATITUDE OF THE POINT.                                C
!                                                                     C
!        XLON : LONGITUDE OF THE POINT.                               C
!                                                                     C
!     OUTPUT :                                                        C
!                                                                     C
!        XI : X COORDINATE OF THE POINT IN MODEL INDEX.               C
!                                                                     C
!        YJ : Y COORDINATE OF THE POINT IN MODEL INDEX.               C
!                                                                     C
!     NOTE ***                                                        C
!                                                                     C
!     A WEST LONGITUDE IS GIVEN BY A NEGATIVE NUMBER; POSITIVE        C
!     ANGLES DENOTE EAST LONGITUDE.                                   C
!                                                                     C
!     A NORTH LATITUDE IS GIVEN BY A POSITIVE NUMBER, AND A NEGATIVE  C
!     NUMBER FOR A SOUTH LATITUDE.                                    C
!-----------------------------------------------------------------------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     REAL        ::  XN, PSI1, FRTY5D, CONV, RLAT1, RLAT2, POLE
     REAL        ::  PSX, CELL, DS, R, C1, RR, C2, CELL2, PHI2
     REAL        ::  XXL, XX, YY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  9/5/03 Modified for more general Lambert Conformal with 2 standard parallels,CLAT1,CLAT2
     FRTY5D=ATAN(1.)
     CONV=45./FRTY5D
     IF (IPROJ.EQ.1)THEN
!23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
     RLAT1=CLAT1/CONV
     RLAT2=CLAT2/CONV
      IF(CLAT1.NE.CLAT2)THEN
        XN=ALOG(COS(RLAT1)/COS(RLAT2))/ALOG(TAN(FRTY5D+.5*RLAT2)/ &
        TAN(FRTY5D+.5*RLAT1))
        ELSE
        XN=SIN(RLAT1)
        ENDIF
     ELSEIF(IPROJ.EQ.2)THEN
      XN=1.
     ENDIF
!     WRITE(*,*)'IPROJ,XN,PH1D= ',IPROJ,XN,PH1D
     POLE=90.
     IF(XLATC.LT.0.)POLE=-90.
     DS = DXBIGDO/1000.
!        * COMPUTE LONGITUDE OF X AXIS (C1), ASSUMING ORIGIN AT NORTH POLE.
        C1 = -XLONC-POLE/XN
!        * COMPUTE DISTANCE BETWEEN NORTH POLE AND CENTER OF DOMAIN.
        PSI1 = CLAT1/CONV
        PSX = (POLE-XLATC)/CONV
     IF(IPROJ.EQ.1)THEN
!     * COMPUTE DISTANCE BETWEEN POINT AND NORTH POLE.
           IF(CLAT1.EQ.CLAT2)THEN
           CELL = REKM/tan(psi1)+REKM*(psi1-xlat/conv)
           RR =  CELL
           C2=REKM/TAN(PSI1)
!     WRITE(*,*)' CLAT1=CLAT2,IPROJ=1'
           ELSE
           CELL = REKM*SIN(PSI1)/XN
           CELL2 = CELL/(TAN(PSI1/2.)**XN)
           C2 =  CELL2 * (TAN(PSX/2.)**XN)
           PHI2 = (POLE-XLAT)/CONV
           RR = CELL2 *TAN(PHI2/2.)**XN
           ENDIF
     ELSEIF(IPROJ.EQ.2)THEN
     C2=REKM*SIN(PSX)*(1.+COS(PSI1))/(1.+COS(PSX))
     PHI2=(POLE-CLAT1)/CONV
     RR=REKM*COS(XLAT/CONV)*(1.+SIN(PHI2))/(1.+SIN(XLAT/CONV))
     ELSE
     WRITE(*,*)' Only projections 1 and 2 (Lambert, Polar Stereog.) are supported, Stopping'
     STOP
     ENDIF
     XXL = XN*(XLON+C1)/CONV
     XX = RR*COS(XXL)
     YY = RR*SIN(XXL)+C2
     XI = FLOAT(IL+1)/2.+XX/DS
     YJ = FLOAT(JL+1)/2.+YY/DS
     IF(INEST1.EQ.1)THEN
     XI=1.+(XI-XNESSTR)*DXBIGDO/DX
     YJ=1.+(YJ-YNESSTR)*DXBIGDO/DX
     ENDIF
     RETURN
     END SUBROUTINE  LLXY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     END MODULE EMISSIONS_WPS_V4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     PROGRAM EMISS_TEST

     USE EMISSIONS_WPS_V4
     IMPLICIT NONE

     CALL  AL2RADM2

     END PROGRAM EMISS_TEST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

