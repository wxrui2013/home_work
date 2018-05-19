#!/bin/sh
#find-ico.sh
#find all <0,0,12,0> cluster and record atom's id
set -x
line_number=0
config_number=0
bn=joinnumber
filenumber=51            #!!must definitly right  when you run this script
rm ico.id.dat
for i in `seq $filenumber`		#extract foundmental informations of sample
do 
echo $i				
	cat $bn.$i   | awk 'NF>=10 {print $1,$2,$3,$7,$8,$9,$10} NF==9 {print $1,$2,$3,$7,$8,$9,0} NF==8 {print $1,$2,$3,$7,$8,0,0}' > $bn.$i.n
	#select element data 
	cat $bn.$i.n | awk 'BEGIN { } $4==0 && $5==0 && $6==12 && $7==0 {print $1}' >ico.id.$i.dat
	line_number=`wc -l ico.id.$i.dat |awk '{print $1}'`
	echo "$i $line_number">>ico.id.dat
	cat ico.id.$i.dat  >>ico.id.dat	
	rm $bn.$i.n ico.id.$i.dat

done
echo "ALL DONE"
