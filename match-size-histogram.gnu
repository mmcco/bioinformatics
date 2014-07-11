set terminal svg
set output 'match-size-histogram.svg'
set title "dm3 Repeat Match Length Distribution" 
set logscale
set style fill solid
set xlabel "Repeat size (bp)"
set ylabel "Number of repeat instances"
plot 'dm3.matchSizeDistr' u 1:2 ti '' smooth freq with boxes
