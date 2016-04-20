for i in {-64..64}
do 
rem=$(( $i % 2 ))
 
if [ $rem -eq 0 ]
then
echo $i
grep -P "\t$i\t" entropy_2para_4.dat > m_$i.dat 
python fitting_scipy.py $i > 
fi
########if [$((i%2)) -eq 0];
########   then 
########	echo $i 
########fi
done 

