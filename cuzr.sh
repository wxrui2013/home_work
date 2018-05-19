#perform a sample of Cu64.5Zr35.5 quench rate is 10^12k/s 
read_restart    restart.cool.1350000	
#units 		metal
#boundary 	p p p
#atom_style      atomic

#lattice 	fcc 3.614
#region		mdbox block 0 20 0 20 0 20    #32000atoms
#create_box	2 mdbox   
#create_atoms	1 region mdbox
#mass            1   63.54
#mass		2   91.22

#set             type  1  type/fraction  2  0.355  1234567   #p463

pair_style      eam/alloy
pair_coeff      * * ZrCu.lammps.eam Zr Cu

group 		cu type 2
group 		zr type 1
timestep        0.002
#velocity        all create 600 4928459 dist gaussian 
run 		0
variable	n equal count(all)
variable 	p equal pe/$n
variable 	v equal vol/$n
variable 	t equal temp 
variable 	s equal step
variable 	pres equal press

neighbor	 3  bin
neigh_modify	every 3 delay 0 check yes

fix	1 all npt temp 300 300 0.25 iso 0 0 0.5
fix     2 all print 1  "$s  $t $v $p ${pres}" file pre-equil.dat title " step temprature volume pe_atom press"

run 	2353	
unfix   1
unfix 	2
#compute MSD of system 
reset_timestep	0
#compute         myRDF all rdf 200 1 1 1 2 2 1 2 2
#compute         rdf all rdf 200
#compute         mymsd all msd com yes average yes
#compute         msdcu cu msd com yes average yes
#compute         msdzr zr msd com yes average yes

#fix             9 all vector 10 c_mymsd[4]
#variable        d equal slope(f_9)/6/(10*dt)
#fix             91 cu vector 10 c_msdcu[4]
#variable        d_cu equal slope(f_91)/6/(10*dt)
#fix             92 zr vector 10 c_msdzr[4]
#variable        d_zr equal slope(f_92)/6/(10*dt)

fix	10 all nvt temp 300 300 0.25 
#fix	21 all ave/time  1000 10 100000 c_myRDF file rdf.300k mode vector
#fix     22 all ave/time  1000 10 100000 c_rdf file rdf1.300k mode vector
#thermo  10
#thermo_style    custom step  temp c_mymsd[4]  v_d c_msdcu[4]  v_d_cu c_msdzr[4] v_d_zr
dump    23 all custom 10000 dump2.s1.300.50.lammpstrj id type x y z
run	500000
write_data 	sample.data
write_restart sample.restart
###################-------heating-------############################################

#################--------equilibrium------###########################################

################--------cooling----------########################################


#fix     1 all npt temp  $t 300 0.25  iso 0 0 0.5 
#fix 	tp all print 1000 "$s $t $v $p ${pres} " file stp.dat title "step temprature volume pe_atom press"
#fix     st2 all print 1000 "$s $t" file s_t2.dat title "step tempprature"
#dump 	1 all atom 5000 dump.atom.$x
#dump    2 all custom 10000 dump.s1.lammpstrj id type xs ys zs
#run     850000
print 	"ALL DONE"

