Program Dshaper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Implicit None
Double Precision, Allocatable, Dimension(:) :: datain,particle_func,xin
Double Precision, Allocatable, Dimension(:) :: padded_datain,padded_dataout
Integer :: datain_size,i
Double Precision :: rescaler,width,m,b
Double Precision :: shift
Double Precision,parameter :: pi=3.141592653589793d0
!!!!!!!!!!!! new params
Double Precision :: delx,xis
Integer :: sinc_pts
Integer :: j,k,padding_multiplier,extension
Logical :: file_name_from_param
Character(30) :: kernel_choice
Double Precision, Allocatable, Dimension(:) :: dataout,sinc
Character(120) :: ufk_name, pdf_from_param_name

print *, "                aDshaper v1.2"
print *, "        by Daniel Olds, Hsiu-Wen Wang"
print *, "             and Katharine Page"

!to bring in options
Call Select_Modes()

Call Read_In_Data()

if (padding_multiplier.gt.1) Call Adding_Padding(padding_multiplier)

if(kernel_choice.eq.'sinc') then
  Call Setup_Sinc_Kernel()
  Call Convolve_Data()
else if(kernel_choice.eq.'sphere') then 
  Call Setup_Sphere_Impulse()
  Call Generate_Padded_Convolution()
end if

!to write data out
Call Write_Out_Function()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
contains
!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Write_Out_Function()
Integer :: i,j

rescaler=1.d0/rescaler

open(unit=22,status="replace",file="shape_function.dat")
open(unit=23,status="replace",file="corrected_pdf.dat")
do i=1,datain_size
    j=i-shift
    if(j.gt.datain_size) j=j-datain_size
    if(j.lt.1) j=j+datain_size
    write(22,*) xin(i),dataout(j)*rescaler
    write(23,*) xin(i),datain(i) - (dataout(j)*rescaler)
end do
close(22)
close(23)


End Subroutine Write_Out_Function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Read_In_Data()
Integer :: num_args,reason
Character(128) :: width_name, datain_name
Logical :: is_there
Real :: junk1, junk2, is_on
Double Precision :: x,y
Character(160) :: first_line
Integer :: c, col_count

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
    else if(width.ne.0.d0) then !pulled something from param file
        print *, "width from parameter file was:",width
    else !file not there and param file missing it
        print *, " "
        print *, " Please provide a width, either in the command line as the second"
        print *, "argument, in the dshaper.parm file, or in a file 'sphere_width.dat'."
        print *, " "
        Call exit()
    end if
else if(file_name_from_param) then !file name given in parameter file
        datain_name=pdf_from_param_name
        print *, "using run name from parameter file:",Trim(datain_name)
else !some other number of arguments
    print *, " "
    print *, "     Please provide the name of the pdf data file, followed by the diameter "
    print *, "      of sphere to convolve with.  Alternativly, place this width value in"
    print *, "    a file in the same directory as the excecutable, called 'sphere_width.dat'"
    print *, " "
    Call exit()
end if

if(.not.file_name_from_param) Call get_command_argument(1,datain_name) !pullling in data(signal) name


!find length of files
Inquire(file=Trim(datain_name),exist=is_there)
!rand_pick_mode=6
if(is_there) then
!first, figure out how many columns there are
open(unit=77,position="rewind",file=Trim(datain_name))
  read(77,'(a)') first_line
  col_count = count([len_trim(first_line) > 0, (first_line(c:c)/=" " .and. &
      first_line(c:c)/="," .and. first_line(c+1:c+1) == " " .or. &
      first_line(c+1:c+1) == "," , c=1,len_trim(first_line)-1)])
close(77)
    
open(unit=77,position="rewind",file=Trim(datain_name))
  !next count up how many data points there are
  is_on=1
  datain_size=0
  do while(is_on.eq.1)
    if(col_count.eq.4) then
        read(77,*,iostat=reason) x,y,junk1,junk2
    else if(col_count.eq.3) then
        read(77,*,iostat=reason) x,y,junk1
    else if(col_count.eq.2) then
        read(77,*,iostat=reason) x,y
    else 
        print *, "Problem : I'm not sure how many columns your data has."
        print *, "   Please contact Dan Olds at oldsdp@ornl.gov"
        Call Exit() 
    end if 
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
Allocate(dataout(1:datain_size))
dataout(:)=0.d0

!now read in datain
datain(:)=0.d0

open(unit=21,position="rewind",file=Trim(datain_name))
do i=1,datain_size
    if(col_count.eq.4) then
        read(21,*) xin(i),datain(i),junk1,junk2
    else if(col_count.eq.3) then
        read(21,*) xin(i),datain(i),junk1
    else if(col_count.eq.2) then
        read(21,*) xin(i),datain(i)
    end if
end do
close(21)

End Subroutine Read_In_Data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Setup_Sinc_Kernel()
delx=xin(4)-xin(3)
sinc_pts=2*datain_size+extension
shift=-Int(mod(Real(sinc_pts),Real(datain_size)))

allocate(sinc(-sinc_pts:sinc_pts))
sinc(:)=0.d0

!show sinc_function
rescaler=1.d0
open(unit=28,status="replace",file="impulse_was.dat")
do i=-sinc_pts,-1
  xis=delx*Dble(i)
  sinc(i)=sin(xis*width)/(xis*width)
  rescaler=rescaler+sinc(i)
  write(28,*) xis,sinc(i)
end do !end of i=-sinc_pts,-1
sinc(0)=1.d0
write(28,*) 0.0,1.d0
do i=1,sinc_pts
  xis=delx*Dble(i)
  sinc(i)=sin(xis*width)/(xis*width)
  rescaler=rescaler+sinc(i)
  write(28,*) xis,sinc(i)
  end do !end of i=1,sinc_pts
close(28)

End Subroutine Setup_Sinc_Kernel
!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Convolve_Data()
Integer :: i,j,k
    
!Allocate(dataout(1:datain_size))
!dataout(:)=0.d0

do i=1,datain_size
  j=i
  do k=-sinc_pts,sinc_pts
    xis=Dble(k)*delx
    dataout(i)=dataout(i)+(datain(j)*sinc(k))
    j=j-1
    if(j.lt.1) j = j+datain_size
  end do !end of j=-sinc_pts,sinc_pts
end do !end of i=1,datain_size

End Subroutine Convolve_Data
!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Adding_Padding(expando)
Integer :: length_was, length_is
Integer :: expando,i,j,k
Double Precision, Allocatable, Dimension(:) :: new_datain,new_xin

length_was=datain_size
length_is=length_was*expando

Allocate(new_datain(1:length_is))
Allocate(new_xin(1:length_is))
new_datain(:)=0.d0
new_xin(:)=0.d0

!copy old data into new_datain
do i=1,length_was
  new_datain(i)=datain(i)
end do

delx = xin(4)-xin(3)
!now expand new_xis
do i=1,expando
  do j=1,length_was
    k=(i-1)*length_was + j
    new_xin(k)=delx*Dble(k)
  end do
end do

!should have it all setup, get rid of old variables, and reallocate them
Deallocate(xin)
Deallocate(datain)

datain_size = length_is

Allocate(xin(1:length_is))
Allocate(datain(1:length_is))

do i=1,length_is
xin(i) = new_xin(i)
datain(i) = new_datain(i)
end do

open(unit=33, status="replace",file = "newxin.dat")
do i=1,length_is
write(33,*) new_xin(i)
end do
close(33)


Deallocate(new_xin)
Deallocate(new_datain)

End Subroutine Adding_Padding
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Select_Modes()
Logical :: is_there, go_on
Integer :: para_file_length, reason
Character(30) :: com_name
!set defaults
width=0.d0
kernel_choice = "sphere"
padding_multiplier = 1
extension = 0
file_name_from_param = .false.

!kernel_choice = "sphere"
!width = 3.2

Inquire(file="dshaper.parm",exist=is_there)
if(is_there) then!parameter file exists, open and read
  para_file_length=0
  print *, "opening and reading parameter file"
  go_on = .true. !set flag of reading more parameters
  open(unit=100, file="dshaper.parm")
  do while(go_on)
    !read next command
    read(100,*,iostat=reason) com_name
    !error catching
    if(reason.gt.0) then !error
      print *, "something went wrong reading the file"
    else if(reason.lt.0) then !end of file
      go_on=.false.
    else !all is good
    para_file_length=para_file_length+1
    !interprit the command with a case statement
    Select Case(com_name)

    Case("Qmin","qmin","QMIN","width","WIDTH","Width","r","R","rad","RAD","Rad","radius","Radius","RADIUS")
    read(100,*) width
    print *, "using characteristic length:",width    

    Case("KERNEL","KERNAL","kernal","kernel","Kernal","Kernel")
    read(100,*) kernel_choice
    print *, "setting kernel_choice to:",kernel_choice   

    Case("padding","Padding","PADDING")
    read(100,*) padding_multiplier
    print *, "setting padding_multiplier to:",padding_multiplier

    Case("sinc extension","Extension","EXTENSION","extension")
    read(100,*) extension
    print *, "setting extension to:",extension
    
    Case("file","FILE","File")
    read(100,*) pdf_from_param_name
    file_name_from_param = .true.
    print *, "should be good"
    
    Case Default
    print *, "Sorry, I don't know what you meant when you said: ",Trim(com_name)
  
    End Select
    End if !end of if(reason.gt.0) then
  end do !end of do while(go_on)
end if !otherwise, just go with defaults as set above


End Subroutine Select_Modes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!shift=Real(particle_size+1.d0)/2.d0
!shift=-Int(mod(Real(sinc_pts),Real(datain_size)))
!shift=-Int(mod(-Real(datain_size),Real(datain_size)))
shift = datain_size - Real(particle_size+1.d0)/2.d0 + 1

do i=1,particle_size

x=xin(i)

if(x.gt.0.d0 .and. x.lt.2*width) then
  particle_func(i)=4.d0*pi*x**2.d0*(1+(x**3.d0/(16.d0*width**3.d0))-(3.d0/4.d0)*(x/width))
else
  particle_func(i)=0.d0
end if
  rescaler=rescaler+particle_func(i)
end do !end of i=-totx/2,totx/2
!rescaler=1.d0/(rescaler)

padded_datain(:)=0.d0
do i=1,datain_size
  padded_datain(i)=datain(i)
end do
do i=datain_size+1,2*datain_size
  padded_datain(i)=(xin(i-datain_size)+xin(datain_size))*(m)+b
end do

End Subroutine Setup_Sphere_Impulse
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Generate_Padded_Convolution()
dataout(:)=0.d0
padded_dataout(:)=0.d0

padded_dataout=convolver(padded_datain,particle_func)

do i=1,datain_size
  dataout(i)=padded_dataout(i)
end do

End Subroutine Generate_Padded_Convolution
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
End Program
