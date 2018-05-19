 module types
	type ::atoms
		integer::id,types		!atoms's type
		real   ::x,y,z		!coordinates of atoms
		integer::s1 		!place numbers of atoms's status
    end type 
    type ::distance
        integer :: id
        real :: ri
    end type
end module
program main

!The program can analyse strucutre of configuration (*.lammpstrj file : id type x y z / id type xs ys zs)
!First . reading file and build subcell
!   2. find neighbor list of atoms and recording them
!   3. sort every list and record them
!   4. compute <ri> of one configuration and record it in one array that store <ri> of all configuration.
!   5. judge the state of configuration file ,if all of the file has been reading then exit loop ,if not then 
!       continue reading and do not exit loop.
!second . compute average <ri> of all configurations.
!   2.  output information about <<ri>> (ensenable average and multiconfiguration average)
!   
!Author : xuerui wei 
!Time   : 2017/7/2
!Institute: BEIJING CONPUTATIONAL SCIENCE RESEARCH CENTER
!
use types
implicit none
character (len=79):: filename ="T300_20Gpa.dump",string1,string2,string3,string4  !change file name in quote""
character (len=79):: filename2,filename3,string5		!create several files with similar style
real :: xlo,xhi,ylo,yhi,zlo,zhi,boxlx,boxly,boxlz	!boundary of system
real :: rij,rx,ry,rz,tmp,tmp2,tmp1,tmpRni,tmp_Ri
real,parameter :: cutoff = 11.0                     !search range 
integer:: timestep,natoms,i,j,k,l,filev
character (len=79):: tmp_config
integer,parameter ::fileid = 1				!id of file
integer,parameter ::pi=3.1415926535, Nline=8009, ri_n=250   !***Nline need change! ,ri_n is max number to collect in neighbor list
real,allocatable :: neighbor_list(:,:)                !id ,i_th neighbor rij (sorted list) 
type(atoms),allocatable :: atom(:)
real ::ri(ri_n,101)=0 ,ri1(ri_n,101)=0,ri2(ri_n,101)=0,ri_ave_multiconfig(ri_n)=0,max_list(500)      !you need to change max_list bigger enough 
real ::Rni(ri_n,101)=0

open(unit = fileid , file = filename)
open(unit = 12, file = "ri_20Gpa.dat")
open(unit = 13, file = "Rni_20Gpa.dat")
open(unit = 14, file = "PDF_1_6_20Gpa.dat")
open(unit = 15, file = "PDF_7_12_20Gpa.dat")
open(unit = 21, file = "ri_Cu_0Gpa.dat")
open(unit = 22, file = "ri_Zr_0Gpa.dat")

write( 12,*) "ri ave_distance"          ! ith nearest atom ,distance 
write (13,*) "Rni ave_distance"         ! Rni=sum(ri)/n :: 1<= i <=n
100 	format (A79)
101		format (I7)
!reading configration
do filev = 1,101
    read (fileid,100) string1
	read (fileid,101) timestep
	read (fileid,100) string2
	read (fileid,101) natoms
	read (fileid,100) string3
	read (fileid,*) xlo,xhi
	read (fileid,*) ylo,yhi
	read (fileid,*) zlo,zhi
	read (fileid,100) string4 
    boxlx=xhi-xlo;boxly=yhi-ylo;boxlz=zhi-zlo
    allocate (atom(natoms))
    allocate (neighbor_list(natoms,ri_n))     !neighbor_list(:,:)
    
    do i=1,natoms
		read (fileid,*) atom(i)%id ,atom(i)%types,atom(i)%x,atom(i)%y,atom(i)%z
        !---unscale coordinate to real space
        atom(i)%x=atom(i)%x*boxlx + xlo
        atom(i)%y=atom(i)%y*boxly + ylo
        atom(i)%z=atom(i)%z*boxlz + zlo
    enddo
    write (*,*) "reading all done!","configuration=" ,filev
    !search neighbor list
    do i= 1,natoms
        l = 0
        max_list(:)=0.0
        do j= 1,natoms                                  !caculate distance between two atoms
            if (j /= i ) then                            ! periodic boundary condition    
                if (abs(atom(i)%x -atom(j)%x) >= boxlx/2 ) then
                    rx=boxlx - abs(atom(i)%x -atom(j)%x)
                else
                    rx=abs(atom(i)%x -atom(j)%x)
                endif
                if (abs(atom(i)%y -atom(j)%y) >= boxly/2 ) then
                    ry=boxly - abs(atom(i)%y -atom(j)%y)
                else
                    ry=abs(atom(i)%y -atom(j)%y)
                endif
                if (abs(atom(i)%z -atom(j)%z) >= boxlz/2 ) then
                    rz=boxlz - abs(atom(i)%z -atom(j)%z)
                else
                    rz=abs(atom(i)%z -atom(j)%z)
                endif
                
                rij=sqrt(rx**2 + ry**2 + rz**2)
                if (rij<=cutoff) then
                    l=l+1
                    max_list(l)=rij
                    if (l>=499) then
                        write(*,*) "you must make max_list(:) bigger or shorten cutoff !"
                        stop
                    endif
                endif
            endif 
        enddo
        !sort neighbor rij list and store it in array neighbor_list(:,:,:)
        call bubble_sort(max_list(1:l),l)
        neighbor_list(i,:)=max_list(1:ri_n)
    enddo
    tmp_Ri=0
    do i=1,ri_n
        tmp1=0
        tmp2=0
        tmp=0
        do j=1,natoms
            if (atom(j)%types == 1) then
                tmp1=tmp1+neighbor_list(j,i)
            else
                tmp2=tmp2+neighbor_list(j,i)
            endif
            tmp=tmp + neighbor_list(j,i)
        enddo
        ri(i,filev) = tmp/natoms
        ri1(i,filev) = tmp1/natoms/2
        ri2(i,filev) = tmp2/natoms/2
        write(12,*) i,ri(i,1)           ! out put single configuration information
        write(21,*) i,ri1(i,1)
        write(22,*) i,ri2(i,1)
        tmp_Ri=tmp_Ri + ri(i,filev)
        Rni(i,filev) = tmp_Ri/i
        write(13,*) i,Rni(i,1)          !output single configuration information
        
    enddo
    write(*,*) "finish one configuration ",filev
    write(*,*) ri(:,1)
    call partical_distribution_function(natoms,1,6,neighbor_list(:,1:6),14)
    call partical_distribution_function(natoms,7,12,neighbor_list(:,7:12),15)
    !########sub section : calculate Dr(n) and DR(n) for different kinds of atoms#########################
   
    
    stop            !multi_configuration compute 
    !release some array
    deallocate(neighbor_list)
    deallocate(atom)
enddo
do i=1,ri_n
    do j=1,101
        ri_ave_multiconfig(i)=ri_ave_multiconfig(i)+ ri(i,j)
    enddo
    ri_ave_multiconfig(i)=ri_ave_multiconfig(i)/real(101)
    write(12,*) i,ri_ave_multiconfig(i)
enddo

close(12)

    end program
    
    subroutine bubble_sort(a,n)
    implicit none
    integer :: n
    real    ::a(n)
    integer ::i,j,k
    real::temp
        do k=1,n
            if(a(k)==0) then
                write(*,*) "sort problem :exist distance 0 between two atoms"
                stop
            endif
        enddo
        do i=n-1,1,-1
            do j=1,i
                if (a(j)>a(j+1)) then
                    temp = a(j)
                    a(j)=a(j+1)
                    a(j+1)=temp
                endif
            enddo
        enddo
        return
    end subroutine
    subroutine partical_distribution_function( natoms, range1,range2 ,neighbor ,filenumber )
    !这一子程序主要是统计第range1个最近邻到第range2个最近邻原子与中心原子之间距离长度的分布形式函数
    !该函数已经经过归一化
    implicit none
    integer:: natoms,range1 ,range2,filenumber
    real :: neighbor(natoms,range1:range2)
    integer,parameter:: bin=100                  !resolation
    real,parameter:: min=1.9 ,max=4.5    !range of partical distribution (you can change this according to your require)
    integer:: i,j,k                         !if the program report error: "over upboundary" ,then you must change variable 'max' bigger   
    real :: dis(bin),dr
    dr=(max-min)/real(bin)
    dis(:)=0
    do i = 1, natoms
        do j= range1,range2
            k=int((neighbor(i,j)-min)/dr)
            dis(k)=dis(k)+1.0
        enddo
    enddo
    pause
    do i=1,bin
    dis(i)=real(dis(i))/real(natoms*(range2-range1))        !had found distribution function (normality)
    enddo
    write(filenumber,*) "r partical_distribution_function"
    do i= 1 ,bin
        write(filenumber,*) min+i*dr,dis(i)
      !  write(*,*) min+i*dr,dis(i),i,dr,min
      !  pause
    enddo
    end subroutine