#!/bin/csh 
# -i in place BUT save bu file with extension .sed
#generically 's/old/new/'    \ = delimitting char

#so need to adjust search and replace strings here
#sed -i .sed 's/include \"include/include \"..\/..\/include/' cable_*
#sed -i .sed 's/old              /new                      /' cable_*


sed -i.sed 's/jxs599/$USER/' *cfg
sed -i.sed 's/p66/$PROJECT/' *cfg
#sed -i.sed 's/(r_1)//' cable_*
#sed -i.sed 's/,r_1//' cable_*
#sed -i.sed 's/ r_1,//' cable_*
#
#sed -i.sed 's/(i_d)//' cable_*
#sed -i.sed 's/,i_d//' cable_*
#sed -i.sed 's/ i_d,//' cable_*
#
#sed -i.sed '/define_dim/d' cable_*
#sed -i.sed '/define_dim/d' casa*

#sed -i.sed 's/ssoil%/ssnow%/' cable_*
#sed -i.sed 's/ssoil/ssnow/' cable_*

#sed -i.sed 's/cable_um_init_subrs/cable_um_init_subrs_mod/' cable_*

#sed -i.sed 's/cbm_module/cable_cbm_module/' cable_*
#sed -i.sed 's/cable_define_types/cable_def_types_mod/' cable_*
#sed -i.sed 's/radiation_module/cable_radiation_module/' cable_*
#sed -i.sed 's/roughness_module/cable_roughness_module/' cable_*
#sed -i.sed 's/offline_driver/cable_offline_driver/' cable_*
#sed -i.sed 's/write_module/cable_write_module/' cable_*
#sed -i.sed 's/cable_cable_/cable_/' cable_*
