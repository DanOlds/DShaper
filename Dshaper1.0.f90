Program Dshaper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Dshaper : by Daniel Olds - dpolds@gmail.com
! A code which is designed to convolve PDF-real 
! space data by a spherical shape function
! in order to generate a characteristic shape function 
! atomistic PDF model data for use
! in data refinement.  Input takes the form:
! ./danshape datain width
! datain = 4-column dataset (Q I dI dQ)
!  width = desired width of gaussian 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Copyright (c) 2015, Los Alamos National Security, LLC
!  All rights reserved.
!  Copyright 2015. Los Alamos National Security, LLC. This software was produced
!  under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory 
!  (LANL), which is operated by Los Alamos National Security, LLC for the U.S. Department 
!  of Energy. The U.S. Government has rights to use, reproduce, and distribute this software.
!  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, 
!  EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  
!  If software is modified to produce derivative works, such modified software should
!  be clearly marked, so as not to confuse it with the version available from LANL.
!
!  Additionally, redistribution and use in source and binary forms, with or without
!  modification, are permitted provided that the following conditions are met:
!  1.     Redistributions of source code must retain the above copyright notice,
!   this list of conditions and the following disclaimer. 
!  2.      Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation and/or
!   other materials provided with the distribution. 
!  3.     Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, 
!   LANL, the U.S. Government, nor the names of its contributors may be used to endorse or
!   promote products derived from this software without specific prior written permission. 
!   
!   THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS"
!   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
!   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
!   IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
!   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
!   TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
!   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
!   THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Implicit None
Double Precision, Allocatable, Dimension(:) :: datain,particle_func,dataout,xin
Double Precision, Allocatable, Dimension(:) :: padded_datain, padded_dataout
Integer :: datain_size,i
Double Precision :: rescaler,width,m,b
Double Precision :: shift
Double Precision,parameter :: pi=3.141592653589793d0

print *, "                 Dshaper v1.0"
print *, "        by Daniel Olds, Hsiu-Wen Wang"
print *, "             and Katharine Page"

Call Read_In_Data()

Call Setup_Sphere_Impulse()

Call Generate_Padded_Convolution()

Call Write_Out_Function()

!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Write_Out_Function()
open(unit=22,status="replace",file="shape_function.dat")
do i=1,datain_size
  if(i-shift+1.ge.1) then !shifted data in range
    write(22,*) xin(Int(i-shift+1)),dataout(i)*rescaler
  end if
end do

do i=Int(datain_size-shift+2),datain_size
  write(22,*) xin(i),padded_dataout(Int(i+shift))*rescaler
end do
close(22)

End Subroutine Write_Out_Function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
function convolver(datain,impulsein)
!datain is the signal array
!impulsein is the noise/impulse array
Double Precision, dimension(:), allocatable :: convolver,tempsum
Double Precision, dimension(:) :: datain,impulsein
integer :: impulse_size, data_size
integer :: i,j,k
!determine size of data-array and impulse-array
data_size=Size(datain)  !i.e. 1000
impulse_size=Size(impulsein) !i.e. 20

!set up new arrays for the convolution
allocate(tempsum(data_size))
allocate(convolver(data_size))

!loop over the data (have to seperate into parts bigger and smaller than impulse-signal

do i=impulse_size,data_size !do i=20,1000
  tempsum(i)=0.d0 !zero this point
  j=i !j=20,19,18,ect (see line 3 down)
  do k=1,impulse_size !loop over entire impulse signal
     tempsum(i)=tempsum(i)+datain(j)*impulsein(k) !actual convolution for i-th point
     j=j-1 !decriment j-th value
  end do !end of k=1,kernalsize
end do !end of i=kernalsize,datasize
        
do i=1,impulse_size !going i=1,20
  tempsum(i)=0.d0 !zero sum
  j=i !j=1 (j=5,4,3...)
  k=1 
  do while (j.gt.0)
     tempsum(i)=tempsum(i)+(datain(j)*impulsein(k)) !actual convoltion for i-th point
     j=j-1 !decriment j
     k=k+1 !increment i
  end do
end do
convolver=tempsum !copy sum into function              
end function convolver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Read_In_Data()
Integer :: num_args,reason
Character(128) :: width_name, datain_name
Logical :: is_there
Real :: junk1, junk2, is_on
Double Precision :: x,y
num_args = command_argument_count()

if(num_args.eq.2) then ! assume second argument is real number, describing width
    Call get_command_argument(2,width_name)
    Read(width_name,998) width
    998 Format(F10.0)
else if(num_args.eq.1) then !assume width is in file called "sphere_width.dat"
    Inquire(file="sphere_width.dat",exist=is_there)
    if(is_there) then
        open(unit=69,file="sphere_width.dat")
        read(69,*) width
        close(69)
    else !file not there
        print *, " "
        print *, "Please provide a width, either in the command line as the second"
        print *, "          argument, or in a file 'sphere_width.dat'."
        print *, " "
        Call exit()
    end if
else !some other number of arguments
    print *, " "
    print *, "     Please provide the name of the pdf data file, followed by the diameter "
    print *, "      of sphere to convolve with.  Alternativly, place this width value in"
    print *, "    a file in the same directory as the excecutable, called 'sphere_width.dat'"
    print *, " "
    Call exit()
end if

Call get_command_argument(1,datain_name) !pullling in data(signal) name

!find length of files
Inquire(file=Trim(datain_name),exist=is_there)
!rand_pick_mode=6
if(is_there) then
open(unit=77,position="rewind",file=Trim(datain_name))
  !next count up how many data points there are
  is_on=1
  datain_size=0
  do while(is_on.eq.1)
    read(77,*,iostat=reason) x,y,junk1,junk2
    if(reason.gt.0) then !some kind of error
      print *, "something went terribly wrong",x
      is_on=0
    else if(reason.lt.0) then !end of file reached
      is_on=0
    else !is good
    datain_size=datain_size+1
    end if !end of reason 
  end do !end of do while is_on
close(77)
else
  print *, "file not found:",Trim(datain_name)
end if 

allocate(datain(1:datain_size))
allocate(padded_datain(1:2*datain_size))
allocate(padded_dataout(1:2*datain_size))
allocate(xin(1:datain_size))
allocate(dataout(1:datain_size))

!now read in datain
datain(:)=0.d0
padded_datain(:)=0.d0

open(unit=21,position="rewind",file=Trim(datain_name))
do i=1,datain_size
  read(21,*) xin(i),datain(i),junk1,junk2
end do
close(21)

End Subroutine Read_In_Data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Generate_Padded_Convolution()
dataout(:)=0.d0
padded_dataout(:)=0.d0

padded_dataout=convolver(padded_datain,particle_func)

do i=1,datain_size
  dataout(i)=padded_dataout(i)
end do
End Subroutine Generate_Padded_Convolution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Setup_Sphere_Impulse()
Double Precision :: delx,x
Integer :: xbounds,xmax, particle_size
!set xbounds (width of function to use for 
xbounds=Floor(datain_size/8.d0)

if(Mod(xbounds,2).ne.0.d0) then
  xbounds=xbounds+1
end if

particle_size=Int(xbounds*2+1) !have to add one for the center (0) values 

rescaler=0.d0
delx=xin(4)-xin(3)

xmax=Ceiling(2.d0*width/delx)
particle_size=xmax+1
allocate(particle_func(1:particle_size))
particle_func(:)=0.d0
shift=Real(particle_size+1.d0)/2.d0

do i=1,particle_size

x=xin(i)

if(x.gt.0.d0 .and. x.lt.2*width) then
  particle_func(i)=4.d0*pi*x**2.d0*(1+(x**3.d0/(16.d0*width**3.d0))-(3.d0/4.d0)*(x/width))
else
  particle_func(i)=0.d0
end if
  rescaler=rescaler+particle_func(i)
end do !end of i=-totx/2,totx/2
rescaler=1.d0/(rescaler)

padded_datain(:)=0.d0
do i=1,datain_size
  padded_datain(i)=datain(i)
end do
do i=datain_size+1,2*datain_size
  padded_datain(i)=(xin(i-datain_size)+xin(datain_size))*(m)+b
end do

End Subroutine Setup_Sphere_Impulse

End Program
