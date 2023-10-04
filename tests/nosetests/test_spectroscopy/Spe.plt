#!/usr/bin/gnuplot

set xl "wave [A]"
set yl "flux [1]"

set zeroaxis

p \
  "Spe.dat" u ($2/1.e-10):3 w p,\

pa -1

set term png small
set out "Spe.png"
rep

