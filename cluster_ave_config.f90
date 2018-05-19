  module types
	type ::atoms
		integer::id0,types		!atoms's type
		real   ::x,y,z		!coordinates of atoms
		integer::s1 		!place numbers of atoms's status
	end type 
	type :: atomcn
		integer :: id1,cn
		integer :: near_atom(22)	! 22 must check before run the program that match with (data.cn.formal)
	end type
	
end module 
program mine 		!caculate number CN12 and CN16  according to neighbor information about atoms in sample 
					!(from voronoi method)and output CN12 ,CN16,and coupling CN12 whith CN16.
use types
implicit none 
character (len=79):: filename ="dump.s1.com.lammpstrj",string1,string2,string3,string4  !change file name in quote""
real :: xlo,xhi,ylo,yhi,zlo,zhi,boxlx,boxly,boxlz,V,rij,tmp_real !boundary of system
integer:: timestep,natoms,i,j,k,k2,l,l2,m,x,filev,atoms_number1,atoms_number2,atoms_number3,other,tmp,tmp2,tmp3,tmp4,tmp_x,N,triang_number
integer:: config,line_number,bonds_number,counter,min_cluster,min_c_id,min_c,max_c_id,max_c,c_id
integer:: c_number(2500),c_number_counter(2500)
character (len=79):: tmp_config
integer,parameter ::fileid = 1				!id of file
integer,parameter ::pi=3.1415926535,Nline=40009    !***Nline need change!
integer,parameter :: n_con = 0				!select configuration in file (secect config equal n_con+1)
integer,parameter :: n_col =34   !column number for record information about ico_cluster :id,c1..c30,cluster_id,state.bonds_number 



type(atoms),allocatable::atom(:)
type(atomcn),allocatable::cn(:)
integer ,allocatable ::atom_status(:)
integer,allocatable ::ico_id(:)
integer,allocatable ::index_id(:)
integer,allocatable ::index2_ico_id(:)		!
integer,allocatable ::ico_full_id(:,:)
integer,allocatable ::cluster(:,:)
real ,allocatable :: cluster_coeff(:)

open (fileid,file = filename)               !**!open file that you want to read,(dump file)
open (unit=2,file="log.txt")				!**!output important information 
open (unit=4,file="log2.dat")				!
open (unit=3,file="data.cn.formal")			!**!input informatin of CN of atoms in sample (the CN from voronoi method)
open (unit=5,file="ico.id.dat")				!**!input id of <0,0,12,0> cluster
open (unit=6,file="dis_cluster.dat")
open (unit=18,file="cluster_cn12.xyz")      !cn=12        
open (unit=20,file="cluster_cn16.xyz")		!cn=16
open (unit=22,file="cluster_cnx.xyz")		!cn=x
open (unit=24,file="cluster_ico.xyz")		!icosehedra
open (unit=26,file="cluster_n.xyz")			!one type cluster that have n center icosehedra
open (unit=28,file="cluster_coeff_ave.txt") !
do i=1,n_con*Nline								!select configuration in file
	read (fileid,*)	
enddo
do i=1,n_con*40000
	read (3,*)
enddo
do i=1,n_con								!ignore n_con number configurations
	read (5,*) config,line_number
	do	j=1,line_number
		read(5,*)
	enddo
enddo
write(28,*)	"cluster_number cluster_coefficient"

do filev=1,21				!----average of all configurations-------
read (fileid,100) string1
read (fileid,101) timestep
read (fileid,100) string2
read (fileid,101) natoms
read (fileid,100) string3
read (fileid,*) xlo,xhi
read (fileid,*) ylo,yhi
read (fileid,*) zlo,zhi
read (fileid,100) string4  
allocate (atom(natoms))
allocate (cn(natoms))
allocate (atom_status(natoms))
allocate (index_id(natoms))
do i=1,natoms
	read (fileid,*) atom(i)%id0 ,atom(i)%types,atom(i)%x,atom(i)%y,atom(i)%z !***make sure %types (1) represent & 

	k=atom(i)%id0
	index_id(k)=i				!build index
!pause														 !cooper atoms 
	read (3,*) cn(i)%id1 ,cn(i)%cn ,cn(i)%near_atom(:)
!	write (*,*) i,cn(i)%id1 ,cn(i)%cn ,cn(i)%near_atom(:)
!	pause
end do							
100 	format (A79)
101		format (I7)
read (5,*) config,line_number
allocate(index2_ico_id(natoms))
allocate(ico_id(line_number))				!read icosehedral id 	
allocate(ico_full_id(line_number,n_col))
allocate(cluster_coeff(line_number))
do	j=1,line_number
		read(5,*) ico_id(j)
		k=ico_id(j)							!build index array
		index2_ico_id(k)=j
enddo
write (*,*) "reading all done!","configuration=" ,filev
!-------------- finished reading file------------------------
!===============output icosehedra atoms=========================
write(24,*) line_number
write (24,*) "icosehedra atoms <0,0,12,0>"
do i=1,line_number
	tmp=ico_id(i)
	tmp2=index_id(tmp)			!index array of dump file 
	if  (atom(tmp2)%types==1) then 
		write (24,*) "Zr" ,atom(tmp2)%x,atom(tmp2)%y,atom(tmp2)%z	
	else
		write (24,*) "Cu" ,atom(tmp2)%x,atom(tmp2)%y,atom(tmp2)%z	
	endif
enddo 
!stop
!----------------calculate relationship of icosehedra cluster-------
!-----------------find icosehedra near one icosehedra----------------------
ico_full_id(:,:)=0
do i=1,line_number
	tmp=ico_id(i)		!id of atoms
	ico_full_id(i,1)=tmp
	N = cn(tmp)%cn			!cordination number of atoms
	l=2				!counter of row match with c1,c2,...cl...c30:record id of ico_atoms near "tmp" icosehedra
	bonds_number=0	
	do j=1,N	
		do k=1,line_number
			if (cn(tmp)%near_atom(j)==ico_id(k)) then
			ico_full_id(i,l)=cn(tmp)%near_atom(j)
			bonds_number=bonds_number+1
			exit
			endif
		enddo
		l=l+1
	enddo
!	ico_full_id(i,32)=1				!just represent state of one atoms(1 means had checked,0 means no)
	ico_full_id(i,33)= bonds_number
!	write(*,*)  i, ico_full_id(i,23) 
enddo
!--------finish selection of icosehedra in near_atoms ---------------------------
!------------------------array method to achieve cluster-------------------------(don't change )=====
!start point : every atoms in the structure has same cluster_id with their neighbor_atoms  
l=0
do i= 1,line_number						!if checked before (0 means didn't checked)
	if (ico_full_id(i,31) == 0) then	!didn't creat cluster id for the atom
		min_c_id = 10000				!start value ensure min_c_id < c_id worked
		min_c = 10000
		max_c_id = 0
		if (ico_full_id(i,33) /= 0) then
			do j = 2 , 30
				if (ico_full_id(i,j) /= 0) then
					tmp = ico_full_id(i,j)
					tmp2 = index2_ico_id(tmp)
					c_id = ico_full_id(tmp2,31)
					if(min_c_id > c_id ) then		!find minimum of cluster id in near_atom
						min_c_id = c_id
					endif
					if (max_c_id < c_id) then		!find maximum of cluster id in near_atom
						max_c_id = c_id
					endif
					if (c_id /=0 ) then				! for !@12 part
						if (min_c > c_id) then
							min_c = c_id
						endif
					endif
				endif
			enddo
			if (min_c_id == max_c_id) then
				!@11
				if (min_c_id == 0) then
					!@13
					l=l+1
					ico_full_id(i,31) = l
					ico_full_id(i,32) = 1
					do j = 2,30
						if (ico_full_id(i,j) /= 0) then		!icosehedra cluster atoms
							tmp = ico_full_id(i,j)
							tmp2 = index2_ico_id(tmp)
							ico_full_id(tmp2,31) = l
						endif
					enddo
				else
					!@x
					ico_full_id(i,31) = min_c_id
					ico_full_id(i,32) = 1
				endif
			else
				!@12
				ico_full_id(i,31) = min_c
				!@16
				do j = 2 , 30				
					if (ico_full_id(i,j) /= 0) then		!icosehedra cluster atoms
						tmp = ico_full_id(i,j)
						tmp2 = index2_ico_id(tmp)
						if (ico_full_id(tmp2,31) == 0 ) then			! global change of cluster's id
							ico_full_id(tmp2,31) = min_c
						elseif(ico_full_id(tmp2,31) /= min_c) then
							do k = 1,line_number
								if (ico_full_id(k,31) == ico_full_id(tmp2,31)) then		!same c_id means one cluster
									ico_full_id(k,31) = min_c
								endif
							enddo
						endif	
					endif
				enddo
			endif
		else 
			l=l+1
			ico_full_id(i,31) = l
			ico_full_id(i,32) = 1
		endif		
	else								!================
		!@2
		!@4		!@8
		min_c_id=10000
		min_c = 10000
		max_c_id=0
		max_c = 0
		do j = 2 , 30
			if (ico_full_id(i,j) /= 0) then
				tmp = ico_full_id(i,j)
				tmp2 = index2_ico_id(tmp)
				c_id = ico_full_id(tmp2,31)
				if (c_id /=0 ) then
					if(min_c > c_id ) then		!find minimum of cluster id in near_atom
						min_c = c_id
					endif
					if (max_c < c_id) then		!find maximum of cluster id in near_atom
						max_c = c_id
					endif
				endif
	!				if(c_id /= 0) then
	!					if (min_c > c_id) then
	!						min_c = c_id
	!					endif
	!				endif
			endif
		enddo
		if (min_c == max_c) then
			!@5
			if (min_c == ico_full_id(i,31)) then
				ico_full_id(i,32) =1
			endif
		else
			!@6 !@7 
			!@9 !@16
			ico_full_id(i,31) = min_c
			ico_full_id(i,32) = 1 
			do j = 2 , 30				
				if (ico_full_id(i,j) /= 0) then		!icosehedra cluster atoms
					tmp = ico_full_id(i,j)
					tmp2 = index2_ico_id(tmp)
					if (ico_full_id(tmp2,31) == 0 ) then			! global change of cluster's id
						ico_full_id(tmp2,31) = min_c
					elseif(ico_full_id(tmp2,31) /= min_c) then
						do k = 1,line_number
							if (ico_full_id(k,31) == ico_full_id(tmp2,31)) then		!same c_id means one cluster
								ico_full_id(k,31) = min_c
							endif
						enddo
					endif	
				endif
			enddo
		endif
	endif
enddo	
!---check atoms-------------
do i = 1,line_number
	if(ico_full_id(i,32) == 0) then
		ico_full_id(i,32) = 1
		do j= 2,30
			if(ico_full_id(i,j) /= 0) then		!ico cluster id in sample 
				tmp=ico_full_id(i,j)
				tmp2=index2_ico_id(tmp)
				c_id=ico_full_id(tmp2,31)
				if (ico_full_id(i,31) /= c_id) then
				!	write(*,*) i,tmp2,c_id, "ok"
					x=ico_full_id(i,31)
					do k= 1,line_number
						if(ico_full_id(k,31) == x) then
							ico_full_id(k,31) = c_id
						endif
					enddo		
				else
					!write(*,*) i ,tmp2 ,c_id ,ico_full_id(i,31), "mistake"
					!exit do,,,
					!pause
				endif
			endif
		enddo
	else
		! write(*,*) i,"ok"
	endif
enddo


!pause
!-----------------------------------
do i=1,line_number
	do j=2,30
		if(ico_full_id(i,j) /=0 ) then
			tmp=ico_full_id(i,j)
			tmp2=index2_ico_id(tmp)
			c_id=ico_full_id(tmp2,31)
			if (ico_full_id(i,31) /= c_id ) then
			!	write (*,*) i,ico_full_id(i,31),tmp2,c_id," /= mistake"
				do k= 1,line_number
					if(ico_full_id(k,31) == c_id) then
					ico_full_id(k,31) = ico_full_id(i,31)
					endif
				enddo
			endif
		endif
	enddo
enddo

!k=0						!check single icosedral
!do i=1,line_number
!	if (ico_full_id(i,33) == 0) then
!		k=k+1
!		write (*,*) ico_full_id(i,31),ico_full_id(i,1),k
!	endif
!enddo
!pause
!==============================================================(don't change except you know it complately)
!----------------------------------compute triangular number of icosahedral atom-----------------
do i=1,line_number
	triang_number = 0
	do j= 2,29					
		if (ico_full_id(i,j) /= 0 ) then
			tmp=ico_full_id(i,j)
			tmp2=index2_ico_id(tmp)			!index number of icosehaderal
			do l = j+1,30
				if (ico_full_id(i,l) /= 0) then
					tmp3 = ico_full_id(i,l)
					do l2 = 2,30			!check relationship between tmp and tmp3 atom
						if(ico_full_id(tmp2,l2) == tmp3) then
							triang_number=triang_number + 1
							exit
						endif
					enddo
				endif
			enddo
		endif
	enddo
	ico_full_id(i,34) = triang_number
!	write(*,*) i,triang_number
!	pause
enddo
!compute "Xi" the clustering coefficient of atom
tmp_real = 0
cluster_coeff(:)=0
do i=1,line_number
	if (ico_full_id(i,33) >= 2) then
	cluster_coeff(i)=real(2*ico_full_id(i,34))/real(ico_full_id(i,33)*(ico_full_id(i,33)-1))
!	write(*,*) ico_full_id(i,1),cluster_coeff(i)
	!................X=2L/N_bonds*(N_bonds-1)
!	pause
	endif
	tmp_real=tmp_real+cluster_coeff(i)
enddo

write(28,*) line_number,tmp_real/real(line_number)
!out put detail information about clustering coefficient
!do i= 1,line_number
!write(28,*) ico_full_id(i,1),ico_full_id(i,31:34),cluster_coeff(i)
!enddo
!------------------------------------------------------------------------------------------------
!stop
			!!!@continue
!--------------------------------------------------------------
					 
write(*,*) l,line_number		!creat k code for cluster but it doesn't means there are k type cluster(
!pause							!some cluster connected and iliminate their code number
allocate (cluster(line_number,3))			!cluster array record :c_id ,atom_number,bonds
cluster(:,:)=0
write(4,*) "atom_id cluster_id atom_state bond_number"
do i=1,line_number
	write(4,*) ico_full_id(i,1),ico_full_id(i,31:33)
!	write(*,*) ico_full_id(i,1),ico_full_id(i,31:33)
	tmp=ico_full_id(i,31)					!type of cluster id code
	cluster(tmp,1)=tmp
	cluster(tmp,2)=cluster(tmp,2)+1        !number of "tmp" id cluster
	cluster(tmp,3)=cluster(tmp,3)+ico_full_id(i,33)
enddo
!pause
c_number(:)=0
do i=1,line_number
!	write(*,*) cluster(i,:)
!	if (cluster(i,2)>=50) then			!select big cluster
!		write(*,*) cluster(i,:)
!	pause 
!	endif
	if (cluster(i,2) /= 0) then
		c_number(cluster(i,2))=c_number(cluster(i,2))+1 
	endif
enddo	
!pause
!write(6,*) "atoms_number cluster_number all_cluster_number percent"
do i=1,1500
	if (c_number(i) /= 0) then
!		write(6,*) i,c_number(i),i*c_number(i),real(i*c_number(i))/real(line_number)
!		write(*,*) i,c_number(i),i*c_number(i),real(i*c_number(i))/real(line_number)
		c_number_counter(i)=c_number_counter(i)+c_number(i)
	endif
enddo 

!pause------------------------------finish selection----------------------------



!-------------------------------------------------------------------------------
!---------output inter-ico-cluster-----------------------------------------------
atom_status(:)=0
	k=203		!test point!!!!!		!atoms number of cluster (k<=max_cluster_id)  cluster_id
	counter=0						!
do i=1,line_number
	if (ico_full_id(i,31) == k ) then
		tmp=ico_full_id(i,1)
!		tmp_x=index_id(tmp)
		if (atom_status(tmp) == 0) then
			counter=counter+1
			atom_status(tmp) = 1			!atom_status is prepared for output cordination of atoms
		endif
		do j=1,cn(tmp)%cn
			if (atom_status(cn(tmp)%near_atom(j)) == 0) then
		!		atom_status(cn(tmp)%near_atom(j)) = 1			!atom_status is prepared for output cordination of atoms
		!		counter=counter+1
			endif
		enddo
!		write(*,*) i,tmp,cn(tmp)%id1,counter ,cluster(k,2)
	endif
enddo
!pause
write(26,*) counter
write(26,*) "center atom cordination of one type cluster"
do i=1,natoms
	if(atom_status(i) ==1) then
		tmp_x=index_id(i)			!index_id(:) connect with dump.sx.300.51.lammpstrj
!		write(*,*) "id",i,tmp_x
		if  (atom(tmp_x)%types==1) then 
			write (26,*) "Zr" ,atom(tmp_x)%x,atom(tmp_x)%y,atom(tmp_x)%z	
		else
			write (26,*) "Cu" ,atom(tmp_x)%x,atom(tmp_x)%y,atom(tmp_x)%z	
		endif
	endif
enddo
!relase some array for next analyses of structure
deallocate (atom)
deallocate (cn)
deallocate (atom_status)
deallocate (index_id)
deallocate(index2_ico_id)
deallocate(ico_id)				!read icosehedral id 	
deallocate(ico_full_id)
deallocate(cluster)
deallocate(cluster_coeff)



enddo			!average of all configurations


!write(6,*) "atoms_number cluster_number all_cluster_number"
!do i=1,1500
!	if (c_number_counter(i) /= 0) then
!		write(6,*) i,real(c_number_counter(i))/real(filev),i*real(c_number_counter(i))/real(filev) !,&
	!	& real(i*c_number(i))/real(line_number)
!		write(*,*) i,c_number(i),i*c_number(i),real(i*c_number(i))/real(line_number)

!	endif
!enddo 
end
!---------------------------------------------------------------------------------------------------------
!=============================you can ignore this if you output inter-icosahedra cluster correctly========
!counter=0
!do i=1,line_number
!	if (ico_full_id(i,31) == k ) then
!		write(*,*) i
!		tmp=ico_full_id(i,1)
!		
!		tmp_x=index_id(tmp)			!index_id(:) connect with dump.sx.300.51.lammpstrj
!		write(*,*) "id",i,"tmp",tmp,tmp_x
!		if  (atom(tmp_x)%types==1) then 
!			write (26,*) "Zr" ,atom(tmp_x)%x,atom(tmp_x)%y,atom(tmp_x)%z	
!		else
!			write (26,*) "Cu" ,atom(tmp_x)%x,atom(tmp_x)%y,atom(tmp_x)%z	
!		endif
!		counter=counter+1
!		N=cn(tmp)%cn
!		do j=1,N 
!			counter=counter+1
!			tmp3=cn(tmp)%near_atom(j)
!			tmp4=index_id(tmp3)
!			write(*,*) i,N,j,tmp3,tmp4
!			if  (atom(tmp4)%types==1) then 
!				write (26,*) "Zr" ,atom(tmp4)%x,atom(tmp4)%y,atom(tmp4)%z	
!			else
!				write (26,*) "Cu" ,atom(tmp4)%x,atom(tmp4)%y,atom(tmp4)%z	
!			endif
!		enddo
!	write(*,*) counter
!	write(*,*) tmp,cn(tmp)%id1,cn(tmp)%near_atom(:)
!	endif
!enddo
!================================================================================================
!------------------------------------------------------------------------------------------------ 
!stop



!================================================================================================
!----------------------------------- section of caculation---------------------------------
!***********the section caculate relationship between CN12 and CN16 then output configuration--- 
!select atoms that cordination number is 12 and put them into one file
!atoms_number1=0;atoms_number2=0;atoms_number3=0;other=0;atom_status(:)=0
!do i = 1, natoms					!caculate numbers of atoms which have different type CN 
!	if ( cn(i)%cn == 12 ) then
!		atoms_number1=atoms_number1+1
!		atoms_number3=atoms_number3+1
!		N = cn(i)%cn
!		do j=1,N
!			tmp2=cn(i)%near_atom(j)
!			if (cn(tmp2)%cn == 16) then				!check status of atoms to avoid repeatition
!				if (atom_status(tmp2) == 0) then	!
!					atoms_number3=atoms_number3+1   !
!					atom_status(tmp2) = 1			!--------------------------------------
!				else
!					atom_status(tmp2) = atom_status(tmp2) + 1
!					if (atom_status(tmp2) > 10) then
!					write (*,*) i,tmp2,atom_status(tmp2)
!					endif
!				endif
!			endif
!		enddo
!	elseif (cn(i)%cn == 16) then
!		atoms_number2=atoms_number2+1
!	else
!		other=other+1
!	endif
!	pause
!enddo
!write(2,*) "id_CN16 num_col_Cn12"
!do i =1,natoms
!	if (atom_status(i) /= 0) then
!		write(2,*) i,atom_status(i)
!	endif
!enddo
!write(18,*) atoms_number1
!write(18,*) "the file contain CN12 atoms without their neighboring atoms" 
!write(20,*) atoms_number2
!write(20,*) "the file contain CN16 atoms without their neighboring atoms "
!write(22,*) atoms_number3
!write(22,*) "the file contain all atoms of coupling cn12 with CN16"
!atom_status(:)=0
!tmp2=0
!do i=1 , natoms
!	if ( cn(i)%cn == 12 ) then
!		k = cn(i)%id1     !"tmp" just represent one instantaneous number
!		tmp = index_id(k)
!		if  (atom(tmp)%types==1) then 
!			write (18,*) "Zr" ,atom(tmp)%x,atom(tmp)%y,atom(tmp)%z
!			write (22,*) "Zr" ,atom(tmp)%x,atom(tmp)%y,atom(tmp)%z
!		else
!			write (18,*) "Cu" ,atom(tmp)%x,atom(tmp)%y,atom(tmp)%z
!			write (22,*) "Cu" ,atom(tmp)%x,atom(tmp)%y,atom(tmp)%z
!		endif
!		!------------output coupling cluster---------------
!		N = cn(i)%cn
!		do j=1,N
!			k2=cn(i)%near_atom(j) !"tmp2" just represent one label
!			tmp2=index_id(k2)
!			if (cn(tmp2)%cn == 16) then
!				if (atom_status(tmp2)==0) then
!					if  (atom(tmp2)%types==1) then 
!						write (22,*) "Zr" ,atom(tmp)%x,atom(tmp)%y,atom(tmp)%z
!					else
!						write (22,*) "Cu" ,atom(tmp)%x,atom(tmp)%y,atom(tmp)%z
!					endif
!					atom_status(tmp2)=1
!				endif
!			endif
!		enddo
!
!	elseif (cn(i)%cn == 16) then
!		k = cn(i)%id1
!		tmp = index_id(k)
!		if  (atom(tmp)%types==1) then 
!			write (20,*) "Zr" ,atom(tmp)%x,atom(tmp)%y,atom(tmp)%z
!		else
!			write (20,*) "Cu",atom(tmp)%x,atom(tmp)%y,atom(tmp)%z
!		endif
!	endif
!
!enddo
!close (18)
!close (20)
!close (22)
!-----------------------------------------------------------------------
!======================================================================================
!end
