# gnuplot input
set term png
set output "eval.png"
set title "Induced charge at ion approach"
set xlabel "ion - electrode distance (Ang)"
set ylabel "electrode charge (e)"
set xrange [0:18]
set yrange [-3:3]
plot  'eval.out' using 1:5 every ::1 title 'cation' with lines lt 1, \
      'eval.out' using 2:6 every ::1 title 'anion' with lines lt 2
