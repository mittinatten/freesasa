set term post eps enhanced solid size 8cm,6cm 10 color

set output "plots/precision.eps"
set encoding utf8

set xlabel "time (ms/atom)"
set ylabel "{/Symbol d} (Å²/atom)"

set style line 1 ps 0.3 pt 7 lc rgb "#CC3333"
set style line 2 ps 0.3 pt 7 lc rgb "#33CC33"

set logscale

p "< grep ^\"STAT SR\" data/precision.dat" u 3:8:4:5:9:10 w boxxyerrorbars ls 1 t "", \
  "< grep ^\"STAT LR\" data/precision.dat" u 3:8:4:5:9:10 w boxxyerrorbars ls 2 t "", \
  "< grep ^\"STAT SR\" data/precision.dat" u 3:8:6:7:11:12 w xyerrorbars ls 1 t "SnR", \
  "< grep ^\"STAT LR\" data/precision.dat" u 3:8:6:7:11:12 w xyerrorbars ls 2 t "LnR"

#"< grep ^LR data/precision.dat" u ($7*1000/$3):8 w d ls 2 t "LnR", \
#"< grep ^SR data/precision.dat" u ($7*1000/$3):8 w d ls 1 t "SnR", \
