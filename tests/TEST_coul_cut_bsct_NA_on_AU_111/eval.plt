# gnuplot input
set term png
set output "eval.png"
set title "Induced charge at ion approach"
set xlabel "ion - electrode distance (Ang)"
set ylabel "electrode charge (e)"
plot  'eval.out' using 1:5 every ::1 title 'cation' with lines lt 1, \
      'eval.out' using 2:6 every ::1 title 'anion' with lines lt 2