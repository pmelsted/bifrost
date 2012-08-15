set terminal png
set output 'threadtime_filter1.png'
set grid nopolar
set xlabel 'threads'
set xtics (1,2,3,4)
set ylabel 'time [s]'
set title "Running time of FilterReads as a function of threads" 
plot 'threadtime_filter1.dat' using 1:2 notitle with lines
