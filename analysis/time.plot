set term post eps enhanced solid size 10cm,10cm 10 color

set encoding utf8
set output "plots/time.eps"

set xlabel "number of atoms"
set ylabel "time (s)"


set multiplot
set size 0.5,0.5 
set size square
set key at graph 0.45,0.95

set style line 1 ps 0.3 pt 7 lc rgb "#CC3333"
set style line 2 ps 0.3 pt 7 lc rgb "#33CC33"
set style line 3 ps 0.3 pt 7 lc rgb "#3333CC"
set style line 4 ps 0.3 pt 7 lc rgb "#CC33CC"
set style line 5 lc rgb "black" lt 1
set style line 6 lc rgb "red" lt 1
set style line 7 lc rgb "blue" lt 1
set logscale

a = 1.41966e-09
b = 1.55958e-05

set origin 0,0.5
set title "Lee and Richards"
p "data/time_LR1_d0.1.dat" t "d = 0.1 Å" ls 1, \
  "data/time_LR1_d0.3.dat" t "d = 0.3 Å" ls 2, \
  "data/time_LR1_d0.5.dat" t "d = 0.5 Å" ls 3, \
  "data/time_LR1_d2.5.dat" t "d = 2.5 Å" ls 4, \
  a*x**2 + b*x t "ax^2 + bx" ls 5, \
  x <=5000 ? b*x : 1/0 t "bx" ls 6, \
  x >= 5000 ? a*x**2 : 1/0 t "ax^2" ls 7

a = 4.10423e-09
b = 2.07026e-05

set origin 0.5,0.5  
set title "Shrake and Rupley"
p "data/time_SR1_2000.dat" t "n = 2000" ls 1, \
  "data/time_SR1_1000.dat" t "n = 1000" ls 2, \
  "data/time_SR1_500.dat" t "n = 500" ls 3, \
  "data/time_SR1_100.dat" t "n = 100" ls 4, \
  a*x**2 + b*x t "ax^2 + bx" ls 5, \
  x <=5000 ? b*x : 1/0 t "bx" ls 6, \
  x >= 5000 ? a*x**2 : 1/0 t "ax^2" ls 7

# Histograms comparing 1- and 2-threaded computations
unset logscale

binwidth=0.05
bin(x,width)=width*floor(x/width)

set xlabel "t(2 threads)/t(1 thread)"
set ylabel "frequency"
unset ytics

set origin 0.02,0
set title "Lee and Richards"
p "< paste ./data/time_LR1_d0.1.dat ./data/time_LR2_d0.1.dat" u (bin($2/$5,binwidth)):(1.0) smooth freq with lines ls 1 t "d = 0.1 Å", \
  "< paste ./data/time_LR1_d0.3.dat ./data/time_LR2_d0.3.dat" u (bin($2/$5,binwidth)):(1.0) smooth freq with lines ls 2 t "d = 0.3 Å", \
  "< paste ./data/time_LR1_d0.5.dat ./data/time_LR2_d0.5.dat" u (bin($2/$5,binwidth)):(1.0) smooth freq with lines ls 3 t "d = 0.5 Å", \
  "< paste ./data/time_LR1_d2.5.dat ./data/time_LR2_d2.5.dat" u (bin($2/$5,binwidth)):(1.0) smooth freq with lines ls 1 t "d = 2.5 Å"

set title "Shrake and Rupley"
set origin 0.52,0
p "< paste ./data/time_SR1_2000.dat ./data/time_SR2_2000.dat" u (bin($2/$5,binwidth)):(1.0) smooth freq with lines t "n = 2000", \
  "< paste ./data/time_SR1_1000.dat ./data/time_SR2_1000.dat" u (bin($2/$5,binwidth)):(1.0) smooth freq with lines t "n = 1000", \
  "< paste ./data/time_SR1_500.dat ./data/time_SR2_500.dat" u (bin($2/$5,binwidth)):(1.0) smooth freq with lines t "n = 500", \
  "< paste ./data/time_SR1_100.dat ./data/time_SR2_100.dat" u (bin($2/$5,binwidth)):(1.0) smooth freq with lines t "n = 100"
