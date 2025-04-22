set terminal svg size 700,700 font "Latin Modern Math,18" background rgb 'white'
set size ratio -1
set view map
#unset xtics
#unset ytics
set title "Ground state"
set out "state_0.svg"
set xlabel "x (nm)"
set ylabel "y (nm)"
set xrange [-4:4]
set yrange [-4:4]
delta=0.08
splot './state_0.dat' u (column(1)*delta-4):(column(2)*delta-4):3 w pm3d notitle

set title "1st excited state"
set out "state_1.svg"
splot './state_1.dat' u (column(1)*delta-4):(column(2)*delta-4):3 w pm3d notitle

set title "2nd excited state"
set out "state_2.svg"
splot './state_2.dat' u (column(1)*delta-4):(column(2)*delta-4):3 w pm3d notitle

set title "3rd excited state"
set out "state_3.svg"
splot './state_3.dat' u (column(1)*delta-4):(column(2)*delta-4):3 w pm3d notitle

set title "4th excited state"
set out "state_4.svg"
splot './state_4.dat' u (column(1)*delta-4):(column(2)*delta-4):3 w pm3d notitle

set title "5th excited state"
set out "state_5.svg"
splot './state_5.dat' u (column(1)*delta-4):(column(2)*delta-4):3 w pm3d notitle
