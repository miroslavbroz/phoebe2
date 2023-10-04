#!/usr/bin/gnuplot

set xl "wave [A]"
set yl "fluxes [1]"

set key bottom

p \
  "<awk '($1==0.00)' test_spectroscopy.out" u ($2*1.0e10):3 w lp,\
  "<awk '($1==0.25)' test_spectroscopy.out" u ($2*1.0e10):3 w lp,\

pa -1

set term png small
set out "test_spectroscopy.png"
rep


