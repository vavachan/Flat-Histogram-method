#!/bin/bash
g++ large_working.cpp
a=1
for T in 0.6 0.8 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7.0 
do
echo $a
./a.out $2 $T $1 > rate.dat
gnuplot << EOF
set terminal png
set output "P_m_$a.png" 
set yrange [0:1]
set title "T=$T"
plot 'rate.dat' u 1:2 w lp 

EOF
((a++))
done
avconv -f image2 -i P_m_%d.png -r 76 -s 800x600 h_e_$2_$1.avi 
rm *.png
