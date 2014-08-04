set terminal svg
set output 'match-size-histogram.svg'
set logscale
set style fill solid
set xlabel "Instance size (bp)"
set ylabel "Number of repeat instances in dm3"
plot 'dm3.match-size' u 1:2 ti '' smooth freq with boxes
