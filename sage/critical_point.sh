for ((H=1 ; H<=200 ; H++))
do
	h=$(echo "scale=2;$H/100.0" | bc)
	for ((P=1 ; P<=100 ; P++))
	do
		p=$(echo "scale=2;$P/100.0" | bc)
#		for  ((T=1 ; T<=300 ; T++))
#		do 
#			t=$(echo "scale=2;$T/100.0" | bc)
			python self.py $p $h   
#		done
	done
done
