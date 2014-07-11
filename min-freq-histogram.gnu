set terminal svg
set output 'min-freq-histogram.svg'
set title "31-mer Distribution among 15-mins in dm3 Repeats"
set logscale
set style fill solid
set xlabel "# of associated 31-mers"
set ylabel "# 15-mins"
plot 'dm3.minfreq' u 1:2 ti '' smooth freq with boxes
