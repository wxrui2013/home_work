#!/bin/sh 
#PBS -N voronoi					
#PBS -l nodes=1:ppn=1
#PBS -o $PBS_O_WORKDIR
export OMP_NUM_THREADS=1                             
#cd $PBS_O_WORKDIR
cd .
bn=dump.s3.300.51.lammpstrj
ac=$bn"_weight_voro";
atomnumber=`sed -n '4p' $bn`
filenumber=`grep -c 'ITEM: TIMESTEP' $bn`
echo $atomnumber
echo $filenumber
mkdir -p $ac
mkdir -p ./$ac/affine
mkdir -p ./$ac/number-join
((atq=$atomnumber+10))
((ath=$atomnumber+9))
##################################################################################
for i in `seq  $filenumber`
do
echo $i
	((c=1+($i-1)*( $atomnumber+9)))
	((b=$c+($atomnumber+8)))
	((d=$b+1))
	sed -n "$d q;$c,$b p" $bn >./$ac/dfile.$i
	cd $ac
echo $ac
  sed -n '6,8p' dfile.$i >./index.$i ##### using in the adffine
  a1=`sed -n '1p' index.$i | cut -d' ' -f 1`
	b1=`sed -n '1p' index.$i | cut -d' ' -f 2` 
	a2=`sed -n '2p' index.$i | cut -d' ' -f 1`
	b2=`sed -n '2p' index.$i | cut -d' ' -f 2`
	a3=`sed -n '3p' index.$i | cut -d' ' -f 1`
	b3=`sed -n '3p' index.$i | cut -d' ' -f 2`
	sed "1,9d" dfile.$i | sort -g >data.$i
	awk '{print $1,$2}' data.$i | sort -g >../filetotal
	awk '$2==1 {print $1,$2,$3,$4,$5,1.59} $2 == 2 {print $1,$2,$3,$4,$5,1.28}' data.$i >datar.$i 
  	cut -d " " -f 1,3-6  ./datar.$i > datav.$i
#	cut -d " " -f 1,3-5  ./data.$i >datav.$i 
	mv datav.$i data1.$i
###################±àÒë³öÀ´µÄvoro++##################################################
    voro++  -p -r -c "%i %s %A  @%i %v @%i %s %n"  $a1 $b1 $a2 $b2 $a3 $b3  data1.$i 
##########################################################################################
    cat data1.$i.vol | cut -d '@' -f 2 | sort -n > voro_vol.$i
    cat data1.$i.vol | cut -d '@' -f 1 | sort -n > data_faces.$i
    cat data1.$i.vol | cut -d '@' -f 3 | sort -n >> data.cn
##################################################################################
sort -g ./data_faces.$i > voro_faces_list.$i
join ../filetotal ./voro_faces_list.$i >./number-join/joinnumber.$i
mv  voro_vol.$i data_faces.$i  ./affine 
cd ..
rm ./$ac/dfile.$i  ./$ac/data1.$i.vol \
 ./$ac/data1.$i ./$ac/data.$i  ./$ac/index.$i ./$ac/voro_faces_list.$i  ./$ac/datar.$i
echo $i
done
##################################################################################
cd ./$ac
cat data.cn | awk '
$2==8 {print $0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}\
$2==9 {print $0,0,0,0,0,0,0,0,0,0,0,0,0,0}\
$2==10 {print $0,0,0,0,0,0,0,0,0,0,0,0,0}\
$2==11 {print $0,0,0,0,0,0,0,0,0,0,0,0}\
$2==12 {print $0,0,0,0,0,0,0,0,0,0,0}\
$2==13 {print $0,0,0,0,0,0,0,0,0,0}\
$2==14 {print $0,0,0,0,0,0,0,0,0}\
$2==15 {print $0,0,0,0,0,0,0,0}\
$2==16 {print $0,0,0,0,0,0,0}\
$2==17 {print $0,0,0,0,0,0}\
$2==18 {print $0,0,0,0,0}\
$2==19 {print $0,0,0,0}\
$2==20 {print $0,0,0}\
$2==21 {print $0,0}\
$2==22 {print $0}
 ' > data.cn.formal

echo all done
