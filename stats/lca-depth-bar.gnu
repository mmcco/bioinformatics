set terminal svg
set output 'lca-depth-bar.svg'
set style fill solid
# set logscale y
set boxwidth 0.7
set xtics 1
set xlabel "Node depth (root has depth 0)"
set ylabel "# Nodes with given depth in hg38"
plot 'hg38.cn-depth' u 1:2 ti '' with boxes
