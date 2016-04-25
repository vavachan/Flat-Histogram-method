
for T in 2.0 1.9 1.8 1.7 1.6 1.5 1.4 1.3 1.2 1.1 1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 
#0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0  
do 
echo $a
python3 trimodal.py $1 $2 $T > rate.dat 
gnuplot << EOF
set yrange [-2:2]
set title "T=$T"
set terminal png
set output "self_$a.png" 
plot 'rate.dat' u 1:2 w lp
set output "self_$a.png" 
replot
EOF
((a++))
done
avconv -f image2 -i self_%d.png -r 76 -s 800x600 rec_$2_$1.avi
rm *.png

