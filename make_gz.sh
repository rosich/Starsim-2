#!/usr/bin/bash

tar -cf starsim.tar  bin/ data/ filters/ output/ starsim.conf make_gz.sh compile_fmodule.sh starsim.py science/ src/  precompiled_ccf precompiled_spectra precompiled_ccf_VIS_IR  *.conf  output_W52/
gzip -9 starsim_2.tar
