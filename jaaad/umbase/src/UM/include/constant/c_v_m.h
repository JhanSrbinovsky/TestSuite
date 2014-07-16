!  Purpose c_v_m.h
!  -------
!
!  Conversion factor from vmr to mmr for each species.
!             vmr*c_species = mmr
!             c_species = m_species/m_air  (m_air=28.97)
!
!  -------------------------------------------------------------------
!
      REAL, PARAMETER :: C_O3P        = 0.5523
      REAL, PARAMETER :: C_O1D        = 0.5523
      REAL, PARAMETER :: C_O3         = 1.657 
      REAL, PARAMETER :: C_NO         = 1.036                           
      REAL, PARAMETER :: C_NO3        = 2.140 
      REAL, PARAMETER :: C_NO2        = 1.588                         
      REAL, PARAMETER :: C_N2O5       = 3.728 
      REAL, PARAMETER :: C_HO2NO2     = 2.727                     
      REAL, PARAMETER :: C_HONO2      = 2.175 
      REAL, PARAMETER :: C_OH         = 0.5868                       
      REAL, PARAMETER :: C_HO2        = 1.139 
      REAL, PARAMETER :: C_H2         = 0.06904                        
      REAL, PARAMETER :: C_H2O2       = 1.174 
      REAL, PARAMETER :: C_CH4        = 0.5523                       
      REAL, PARAMETER :: C_C          = 0.4142 
      REAL, PARAMETER :: C_CO         = 0.9665 
      REAL, PARAMETER :: C_CO2        = 1.5188          
      REAL, PARAMETER :: C_HCHO       = 1.036 
      REAL, PARAMETER :: C_MeOO       = 1.622                       
      REAL, PARAMETER :: C_H2O        = 0.6213 
      REAL, PARAMETER :: C_MeOOH      = 1.657                      
      REAL, PARAMETER :: C_HONO       = 1.622 
      REAL, PARAMETER :: C_O2         = 1.105                         
      REAL, PARAMETER :: C_N2         = 0.9665 
      REAL, PARAMETER :: C_C2H6       = 1.036                        
      REAL, PARAMETER :: C_EtOO       = 2.106 
      REAL, PARAMETER :: C_EtOOH      = 2.140                      
      REAL, PARAMETER :: C_MeCHO      = 1.519
      REAL, PARAMETER :: C_TOTH       = 1.000                        
      REAL, PARAMETER :: C_MeCO3      = 2.589 
      REAL, PARAMETER :: C_PAN        = 4.177                       
      REAL, PARAMETER :: C_C3H8       = 1.519 
      REAL, PARAMETER :: C_PrOO       = 2.589                       
      REAL, PARAMETER :: C_PrOOH      = 2.623 
      REAL, PARAMETER :: C_EtCHO      = 2.002                     
      REAL, PARAMETER :: C_EtCO3      = 3.072 
      REAL, PARAMETER :: C_Me2CO      = 2.002                     
      REAL, PARAMETER :: C_MeCOCH2OO  = 3.072 
      REAL, PARAMETER :: C_MeCOCH2OOH = 3.107            
      REAL, PARAMETER :: C_PPAN       = 4.660 
      REAL, PARAMETER :: C_MeONO2     = 2.658                     
      REAL, PARAMETER :: C_N          = 0.48325 
      REAL, PARAMETER :: C_H          = 0.03452                         
      REAL, PARAMETER :: C_N2O        = 1.5188 
      REAL, PARAMETER :: C_CFCl3      = 4.7480                     
      REAL, PARAMETER :: C_CF2Cl2     = 4.1783 
      REAL, PARAMETER :: C_ClO        = 1.7784                    
      REAL, PARAMETER :: C_HCl        = 1.2604
      REAL, PARAMETER :: C_ClONO2     = 3.3668                    
      REAL, PARAMETER :: C_HOCl       = 1.8129
      REAL, PARAMETER :: C_OClO       = 2.3309                     
      REAL, PARAMETER :: C_BrO        = 3.315
      REAL, PARAMETER :: C_BrONO2     = 4.9034                     
      REAL, PARAMETER :: C_HBr        = 2.7970
      REAL, PARAMETER :: C_HOBr       = 3.3495                      
      REAL, PARAMETER :: C_BrCl       = 3.9884
      REAL, PARAMETER :: C_MeBr       = 3.2805                     
      REAL, PARAMETER :: C_SO2        = 2.2112 
      REAL, PARAMETER :: C_SO3        = 2.7615                      
      REAL, PARAMETER :: C_Me2S       = 2.145
      REAL, PARAMETER :: C_DMS        = 2.145
      REAL, PARAMETER :: C_OCS        = 2.0711                      
      REAL, PARAMETER :: C_SAD        = 4.1255
      REAL, PARAMETER :: C_MSA        = 3.317
      REAL, PARAMETER :: C_S          = 1.1046                           
      REAL, PARAMETER :: C_H2SO4      = 3.385
      REAL, PARAMETER :: C_CF2ClCFCl2 = 6.4722              
      REAL, PARAMETER :: C_CHF2Cl     = 2.9858
      REAL, PARAMETER :: C_MeCCl3     = 4.6082                 
      REAL, PARAMETER :: C_CCl4       = 5.3158
      REAL, PARAMETER :: C_MeCl       = 1.7432                     
      REAL, PARAMETER :: C_CF2ClBr    = 5.7128
      REAL, PARAMETER :: C_CF3Br      = 5.1432                 
      REAL, PARAMETER :: C_Cl         = 1.2261 
      REAL, PARAMETER :: C_Cl2O2      = 3.5568
      REAL, PARAMETER :: C_Br         = 2.7627
      REAL, PARAMETER :: C_C5H8       = 2.3473 
      REAL, PARAMETER :: C_ISO2       = 4.0387       
      REAL, PARAMETER :: C_ISOOH      = 4.0732 
      REAL, PARAMETER :: C_ISON       = 5.3504      
      REAL, PARAMETER :: C_MACR       = 2.4163 
      REAL, PARAMETER :: C_MACRO2     = 4.1077       
      REAL, PARAMETER :: C_MACROOH    = 4.1422 
      REAL, PARAMETER :: C_MPAN       = 5.0742       
      REAL, PARAMETER :: C_HACET      = 2.5544 
      REAL, PARAMETER :: C_MGLY       = 2.4853       
      REAL, PARAMETER :: C_NALD       = 3.6244 
      REAL, PARAMETER :: C_HCOOH      = 1.5878       
      REAL, PARAMETER :: C_MECO3H     = 2.6234 
      REAL, PARAMETER :: C_MECO2H     = 2.0711      
      REAL, PARAMETER :: C_NH3        = 0.5879

!     molecular masses in g/mol of emitted species, 
!     for budget calculations

      REAL, PARAMETER :: m_ch4     =  16. 
      REAL, PARAMETER :: m_co      =  28.
      REAL, PARAMETER :: m_hcho    =  30. 
      REAL, PARAMETER :: m_c2h6    =  30.
      REAL, PARAMETER :: m_c3h8    =  44. 
      REAL, PARAMETER :: m_mecho   =  44.
      REAL, PARAMETER :: m_no2     =  46. 
      REAL, PARAMETER :: m_me2co   =  58.
      REAL, PARAMETER :: m_isop    =  68. 
      REAL, PARAMETER :: m_no      =  30.
      REAL, PARAMETER :: m_c       =  12.

!     molecular masses of stratospheric species, for which surface 
!     mmrs are prescribed

      REAL, PARAMETER :: m_hcl        =  36.5
      REAL, PARAMETER :: m_n2o        =  44.
      REAL, PARAMETER :: m_clo        =  51.5
      REAL, PARAMETER :: m_hocl       =  52.5
      REAL, PARAMETER :: m_oclo       =  67.5
      REAL, PARAMETER :: m_clono2     =  97.5
      REAL, PARAMETER :: m_cf2cl2     = 121. 
      REAL, PARAMETER :: m_cfcl3      = 137.5
      REAL, PARAMETER :: m_hbr        =  81. 
      REAL, PARAMETER :: m_mebr       =  95.
      REAL, PARAMETER :: m_bro        =  96. 
      REAL, PARAMETER :: m_hobr       =  97.
      REAL, PARAMETER :: m_brcl       = 115.5
      REAL, PARAMETER :: m_brono2     = 142.
      REAL, PARAMETER :: m_cf2clcfcl2 = 187.5
      REAL, PARAMETER :: m_chf2cl     =  86.5
      REAL, PARAMETER :: m_meccl3     = 133.5
      REAL, PARAMETER :: m_ccl4       = 154.
      REAL, PARAMETER :: m_mecl       =  50.5
      REAL, PARAMETER :: m_cf2clbr    = 165.5
      REAL, PARAMETER :: m_cf3br      = 149.
      REAL, PARAMETER :: m_ocs        =  60. 
      REAL, PARAMETER :: m_dms        =  62.1 
      REAL, PARAMETER :: m_me2s       =  62.1
      REAL, PARAMETER :: m_msa        =  96.1
      REAL, PARAMETER :: m_s          =  32.07
      REAL, PARAMETER :: m_so2        =  64.06 
      REAL, PARAMETER :: m_so4        =  96.06
      REAL, PARAMETER :: m_nh3        =  17.03

      REAL, PARAMETER :: m_cl         = 35.5 
      REAL, PARAMETER :: m_cl2o2      = 103.
      REAL, PARAMETER :: m_br         = 80.
      REAL, PARAMETER :: m_h2         = 2.016
      REAL, PARAMETER :: m_h2o        = 18.0
      REAL, PARAMETER :: m_mecoch2ooh = 90.0
      REAL, PARAMETER :: m_isooh      = 118.0
      REAL, PARAMETER :: m_mpan       = 147.0
      REAL, PARAMETER :: m_ppan       = 135.0         
!  -------------------------------------------------------------------
