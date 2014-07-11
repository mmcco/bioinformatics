set terminal svg
set output 'lca-depth-bar.svg'
set title "Distr. of LCA Node Depth by Kmer in dm3" 
set style fill solid
set logscale y
set boxwidth 0.7
set xtics 1
set xlabel "LCA depth (root has depth 0)"
set ylabel "Number of 31-mers with given depth"
plot 'dm3.depthfreq' u 1:2 ti '' with boxes
