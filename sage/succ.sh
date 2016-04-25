
for  ((T=1 ; T<=300 ; T++))
do 
	t=$(echo "scale=2;$T/100.0" | bc)
	python3 trimodal.py $1 $2 $t  
done
