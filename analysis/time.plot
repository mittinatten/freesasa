set term post eps enhanced solid size 15cm,7cm 10 color

set encoding utf8
set output "plots/time.eps"

set xlabel "number of atoms"
set ylabel "time (s)"


set multiplot layout 1,2
set size square
set key at graph 0.3,0.95

set style line 1 ps 0.3 pt 7 lc rgb "#CC3333"
set style line 2 ps 0.3 pt 7 lc rgb "#33CC33"
set style line 3 ps 0.3 pt 7 lc rgb "#3333CC"
set style line 4 ps 0.3 pt 7 lc rgb "#CC33CC"
set style line 5 lc rgb "black" lt 1
set style line 6 lc rgb "red" lt 1
set style line 7 lc rgb "blue" lt 1
set logscale

a = 2.73988e-09
b = 1.23509e-05

set title "Lee and Richards"
p "data/time_L_d0.1.dat" t "d = 0.1 Å" ls 1, \
  "data/time_L_d0.3.dat" t "d = 0.2 Å" ls 2, \
  "data/time_L_d0.5.dat" t "d = 0.5 Å" ls 3, \
  "data/time_L_d2.5.dat" t "d = 2.5 Å" ls 4, \
  a*x**2 + b*x t "ax^2 + bx" ls 5, \
  x <=5000 ? b*x : 1/0 t "bx" ls 6, \
  x >= 5000 ? a*x**2 : 1/0 t "ax^2" ls 7

a = 7.63584e-09
b = 2.88918e-05

set title "Shrake and Rupley"
p "data/time_S_n2000.dat" t "n = 2000" ls 1, \
  "data/time_S_n1000.dat" t "n = 1000" ls 2, \
  "data/time_S_n500.dat" t "n = 500" ls 3, \
  "data/time_S_n100.dat" t "n = 100" ls 4, \
  a*x**2 + b*x t "ax^2 + bx" ls 5, \
  x <=5000 ? b*x : 1/0 t "bx" ls 6, \
  x >= 5000 ? a*x**2 : 1/0 t "ax^2" ls 7
