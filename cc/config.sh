#!/bin/bash
N=$1
z=$((N**3))
echo $z
for (( c=0; c<=$z; c++ )) 
do 
python cc.py $c $N 
done 
