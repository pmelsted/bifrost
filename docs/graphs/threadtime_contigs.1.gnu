set terminal png
set output 'threadtime_contigs.1.png'
set grid nopolar
set xlabel 'threads'
set ylabel 'time [s]'
set title "Running time of BuildContigs as a function of threads" 
plot 'threadtime_contigs.1.dat' using 1:2 notitle with lines
