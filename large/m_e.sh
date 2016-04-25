#!/bin/bash
g++ large_working.cpp
./a.out 4 $1 $2 > output.dat
gnuplot<< EOF
set term png
set output "picture.png"
plot 'output.dat' u 1:2 w lp
EOF
 picture.png
