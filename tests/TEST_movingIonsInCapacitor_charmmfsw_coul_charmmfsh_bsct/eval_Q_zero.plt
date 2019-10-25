# gnuplot input
set term png
set output "eval_Q_zero.png"
set title "Induced charge at ion approach"
set xlabel "ion - electrode distance (Ang)"
set ylabel "electrode charge (e)"
set yrange [-3:3]
set xrange [0:18]
plot  'eval_Q_zero.out' using 1:5 every ::1 title 'cation' with lines lt 1, \
      'eval_Q_zero.out' using 2:6 every ::1 title 'anion' with lines lt 2
