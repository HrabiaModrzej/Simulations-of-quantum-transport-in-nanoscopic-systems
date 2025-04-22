set terminal svg size 700,700 font "Latin Modern Math,18" background rgb 'white'
set output "E_omega.svg"
set size square
set xrange [50:500]
set xlabel "ℏω_x (meV)"
set ylabel "Eigen energy (meV)"
set key top left
plot for [i=2:11] './../data/E_omega.dat' u 1:(column(i)) w l lw 2 title "state ".(i-2)."" ,\
for [j=2:5] './../data/exact.dat' u 1:(column(j)) w l lw 2 lc rgb "black" dashtype 3 notitle