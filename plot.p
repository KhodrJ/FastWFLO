set terminal png size 750,750
set output 'plot.png'

set style arrow 1 linecolor rgb "white"
set view map

# unset colorbox
set xlabel 'x'
set ylabel 'y'
set title 'Velocity Field over Terrain'
splot 'field.txt' u 1:2:5 with pm3d t '', 'field.txt' u 1:2:(1):3:4:(1) w vec arrowstyle 1  t ''
