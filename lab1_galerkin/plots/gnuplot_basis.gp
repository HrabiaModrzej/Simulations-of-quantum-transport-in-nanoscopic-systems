set terminal svg size 700,700 font "Latin Modern Math,18" background rgb 'white'
set output "basis_0.svg"
set size square
set view map
set title "Gaussian k=0"
set xrange [-4:4]
set xlabel "x (nm)"
set ylabel "y (nm)"
splot './../data/gaussian_0.dat' u ($1*0.08 - 4):($2*0.08 -4):3 w pm3d notitle

set out "basis_8.svg"
set title "Gaussian k=8"
splot './../data/gaussian_8.dat' u ($1*0.08 - 4):($2*0.08 -4):3 w pm3d notitle

set out "basis_9.svg"
set title "Gaussian k=9"
splot './../data/gaussian_9.dat' u ($1*0.08 - 4):($2*0.08 -4):3 w pm3d notitle
