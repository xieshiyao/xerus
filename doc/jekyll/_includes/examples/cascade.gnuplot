#!/usr/bin/gnuplot

set term pngcairo size 600, 400
unset key
set grid
set output "cascade.png"
set xlabel "timestep"
set ylabel "mean concentration"

plot 'mean.dat' u 1 w l, \
     'mean.dat' u 2 w l, \
     'mean.dat' u 3 w l, \
     'mean.dat' u 4 w l, \
     'mean.dat' u 5 w l, \
     'mean.dat' u 6 w l, \
     'mean.dat' u 7 w l, \
     'mean.dat' u 8 w l, \
     'mean.dat' u 9 w l, \
     'mean.dat' u 10 w l
