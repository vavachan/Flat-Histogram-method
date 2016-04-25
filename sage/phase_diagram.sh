for ((H=1 ; H<=100 ; H++))
do
	h=$(echo "scale=3;$H/100.0" | bc)
########	for ((P=1 ; P<=100 ; P++))
########	do
########		p=$(echo "scale=2;$P/100.0" | bc)
		for  ((T=1000 ; T > 1 ; T--))
		do 
			t=$(echo "scale=3;$T/1000.0" | bc)
			a=$(python3 trimodal.py $1 $h $t) 
			flag=$(bc <<< "$a != 0.0")
			beta=$(bc <<< "scale=2;1.0/$t")
			if  (($flag)) ; then
				echo $1 $h $a $t
				break
			fi
			
				
		done
#	done
done
