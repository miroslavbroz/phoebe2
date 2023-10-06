#!/usr/bin/gnuplot

km = 1.e3		# m
cm = 1.e-2		# m
ang = 1.e-10		# m
erg = 1.e-7		# J
c = 299792458.		# m s^-1
h = 6.67e-34		# J s
k = 1.38e-23		# J K^-1
sigma = 5.68e-8		# W m^-2 K^-4
R_S = 6.957e8		# m
au = 1.496e11		# m
T = 6000.0		# K
R = 1.*R_S
d = 1.*au

B_lambda(lambda, T) = 2.0*h*c**2/lambda**5 / (exp(h*c/(lambda*k*T))-1.0)
B(T) = sigma*T**4

set xl "wave [A]"
set yl "fluxes [W m^{-2} m^{-1}], at Earth"

#set xr [3000:10000]
set yr [0:]
set zeroaxis
set key bottom

tmp = 1.0
call "line.plt" "Halpha" 6562.81

p \
  "<awk '($1==0.00)' test_spectroscopy2.out" u ($2*1.0e10):3 w lp,\
  "<awk '($1==0.25)' test_spectroscopy2.out" u ($2*1.0e10):3 w lp,\
  pi*B_lambda(x*ang, T)*(R/d)**2 w l lw 3 dt 2 lc 'orange' t "pi B_{lambda}(T) (R/d)^2"

pa -1

set term png small
set out "test_spectroscopy2.png"
rep


