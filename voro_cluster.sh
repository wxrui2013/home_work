#!/bin/sh 
#PBS -N cluster 					
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o $PBS_O_WORKDIR
#program 
#	Purpose of this program is to extract CN and cluster's information in output file of voronoi
#Histary
#writer : wei xuerui ;Time :2016/7/11 ; school : beijing normal university.
#
export OMP_NUM_THREADS=1                             
set -x
cd .
bn=joinnumber
filenumber=51            #!!must definitly right  when you run this script
rm $bn.zr.* $bn.cu.* $bn.n.*
for i in `seq $filenumber`		#extract foundmental informations of sample
do 
echo $i				
	cat $bn.$i   | awk 'NF>=10 {print $1,$2,$3,$7,$8,$9,$10} NF==9 {print $1,$2,$3,$7,$8,$9,0} NF==8 {print $1,$2,$3,$7,$8,0,0}' > $bn.$i.n
	cat $bn.$i.n | awk '$2==1 {print} ' > $bn.$i.zr
	cat $bn.$i.n | awk '$2==2 {print} ' > $bn.$i.cu
	cat $bn.$i.zr | awk '{print $3}' | sort -n  >>$bn.zr.cn
	cat $bn.$i.cu | awk '{print $3}' | sort -n  >>$bn.cu.cn
	cat $bn.$i.n  | awk '{print $3}' | sort -n  >>$bn.n.cn
	cat $bn.$i.n  | awk '{print $4,$5,$6,$7}' |sort -n >> $bn.n.cluster
	cat $bn.$i.zr | awk '{print $4,$5,$6,$7}' |sort -n >> $bn.zr.cluster
	cat $bn.$i.cu | awk '{print $4,$5,$6,$7}' |sort -n >> $bn.cu.cluster  
done
echo $i
cat $bn.zr.cn | sort -n | uniq -c | awk 'BEGIN{print "CN count"} {print $2,$1/51}' >ave.zr.cn
cat $bn.cu.cn | sort -n | uniq -c | awk 'BEGIN{print "CN count"} {print $2,$1/51}' >ave.cu.cn
cat $bn.n.cn  | sort -n | uniq -c | awk 'BEGIN{print "CN count"} {print $2,$1/51}' >ave.n.cn
cat $bn.n.cluster  | sort -n | uniq -c | awk '$1>=5100 {print}' | awk 'BEGIN{print "<i3,i4,i5,i6> number"} 
{print "<"$2","$3","$4","$5">",$1/51}' > ave.51.n.cluster
cat $bn.cu.cluster | sort -n | uniq -c | awk '$1>=3315 {print}' | awk 'BEGIN{print "<i3,i4,i5,i6> number"} 
{print "<"$2","$3","$4","$5">",$1/51}' > ave.51.cu.cluster
cat $bn.cu.cluster | sort -n | uniq -c | awk '$1>=3315 {print}' | awk 'BEGIN{print "<i3,i4,i5,i6> number"} 
$2+$3+$4+$5==11 {print "<"$2","$3","$4","$5">",$1/51}' >cluster.cu.cn11
cat $bn.cu.cluster | sort -n | uniq -c | awk '$1>=3315 {print}' | awk 'BEGIN{print "<i3,i4,i5,i6> number"} 
$2+$3+$4+$5==12 {print "<"$2","$3","$4","$5">",$1/51}' >cluster.cu.cn12
cat $bn.cu.cluster | sort -n | uniq -c | awk '$1>=3315 {print}' | awk 'BEGIN{print "<i3,i4,i5,i6> number"} 
$2+$3+$4+$5==13 {print "<"$2","$3","$4","$5">",$1/51}' >cluster.cu.cn13
cat $bn.zr.cluster | sort -n | uniq -c | awk '$1>=1785 {print}' | awk 'BEGIN{print "<i3,i4,i5,i6> number"}
 {print "<"$2","$3","$4","$5">",$1/51}' > ave.51.zr.cluster
cat $bn.zr.cluster | sort -n | uniq -c | awk '$1>=1785 {print}' | awk 'BEGIN{print "<i3,i4,i5,i6> number"}
$2+$3+$4+$5==15 {print "<"$2","$3","$4","$5">",$1/51}' > cluster.zr.cn15
cat $bn.zr.cluster | sort -n | uniq -c | awk '$1>=1785 {print}' | awk 'BEGIN{print "<i3,i4,i5,i6> number"}
$2+$3+$4+$5==16 {print "<"$2","$3","$4","$5">",$1/51}' > cluster.zr.cn16
cat $bn.zr.cluster | sort -n | uniq -c | awk '$1>=1785 {print}' | awk 'BEGIN{print "<i3,i4,i5,i6> number"}
$2+$3+$4+$5==17 {print "<"$2","$3","$4","$5">",$1/51}' > cluster.zr.cn17

mkdir -p ./details
mv $bn.*.n ./details
mv $bn.*.zr ./details
mv $bn.*.cu ./details
echo $i 
echo all done


