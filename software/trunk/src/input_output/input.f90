module input
use strings       

private :: variable_read_double,variable_read_int,variable_read_string
private :: variable_read_vec_double,variable_read_vec_int,variable_read_vec_string
interface variable_read
    module procedure variable_read_double
    module procedure variable_read_int
    module procedure variable_read_string
end interface

interface variable_read_vec
    module procedure variable_read_vec_double
    module procedure variable_read_vec_int
    module procedure variable_read_vec_string
end interface


!---------------------------------------------------------------------------
contains 
!---------------------------------------------------------------------------


subroutine file_read(inp_line,filename,args,inp_line_true)


implicit none

integer i,narg,inp_line_true,ios,inp_line
character(len=80),allocatable :: input_file(:)
character(len=40),intent(inout) :: args(:,:)
character(len=*),intent(in) :: filename

allocate(input_file(inp_line))

open(1,file=filename,status='old',action='read')

do i = 1,inp_line
!    print *, i
    call readline(1,input_file(i),ios)
!    print *, input_file(i)
    call parse(input_file(i),'=',args(:,i),narg)

    if(ios.ne.0) then 
        inp_line_true = i
        print *, 'Reached end of file. Final input index = ',inp_line_true
        exit
    endif
enddo

!inp_line = inp_line_true

close(1)

end subroutine file_read


!--------------------------------------------------------------------------


subroutine variable_read_double(args,string,var,inp_line)

integer, intent(in) :: inp_line
character(len=40),intent(in) :: args(:,:)
character(len=*),intent(in) :: string
double precision,intent(out) :: var

integer i, ios
logical found

found = .false.
do i = 1,inp_line
    if(string.eq.trim(adjustl(args(1,i)))) then
        call value(args(2,i),var,ios)
        found = .true.
        exit
    elseif(i.eq.inp_line.and..not.found) then
        print *, 'Variable ',string,' not found. Check input file '
    endif
enddo

end subroutine variable_read_double

!--------------------------------------------------------------------------


subroutine variable_read_int(args,string,var,inp_line)

integer, intent(in) :: inp_line
character(len=40),intent(in) :: args(:,:)
character(len=*),intent(in) :: string
integer,intent(out) :: var

integer i, ios
logical found

found = .false.
do i = 1,inp_line
    if(string.eq.trim(adjustl(args(1,i)))) then
        call value(args(2,i),var,ios)
        found = .true.
        exit
    elseif(i.eq.inp_line.and..not.found) then
        print *, 'Variable ',string,' not found. Check input file '
    endif
enddo



end subroutine variable_read_int


!--------------------------------------------------------------------------


subroutine variable_read_string(args,string,var,inp_line)

integer, intent(in) :: inp_line
character(len=40),intent(in) :: args(:,:)
character(len=*),intent(in) :: string
character(len=*),intent(out) :: var

integer i, ios
logical found

found = .false.
do i = 1,inp_line
    if(string.eq.trim(adjustl(args(1,i)))) then
        var = args(2,i)
        found = .true.
        exit
    elseif(i.eq.inp_line.and..not.found) then
        print *, 'Variable ',string,' not found. Check input file '
    endif
enddo



end subroutine variable_read_string


!---------------------------------------------------------------------------

subroutine variable_read_vec_double(args,string,var,inp_line,num_species) 

integer, intent(in) :: inp_line,num_species
character(len=*),intent(in) :: args(:,:)
character(len=*),intent(in) :: string
double precision,intent(out) :: var(:)

integer i, ios
character(len=40) temp
logical found

found = .false.

do i = 1,inp_line
    if(string.eq.trim(adjustl(args(1,i)))) then 
        do j = 1,(num_species-1)
            call split(args(2,i),',',temp)
            call value(temp,var(j),ios)
        enddo
        found = .true.
        call value(args(2,i),var(num_species),ios)
    elseif(i.eq.inp_line.and..not.found) then
        print *, 'Variable ',string,' not found. Check input file '
    endif
enddo

end subroutine variable_read_vec_double

!---------------------------------------------------------------------------

subroutine variable_read_vec_int(args,string,var,inp_line,num_species) 

integer, intent(in) :: inp_line,num_species
character(len=*),intent(in) :: args(:,:)
character(len=*),intent(in) :: string
integer,intent(out) :: var(:)

integer i, ios
character(len=40) temp
logical found

found = .false.

do i = 1,inp_line
    if(string.eq.trim(adjustl(args(1,i)))) then 
        do j = 1,(num_species-1)
            call split(args(2,i),',',temp)
            call value(temp,var(j),ios)
        enddo
        found = .true.
        call value(args(2,i),var(num_species),ios)
    elseif(i.eq.inp_line.and..not.found) then
        print *, 'Variable ',string,' not found. Check input file '

    endif
enddo

end subroutine variable_read_vec_int

!---------------------------------------------------------------------------

subroutine variable_read_vec_string(args,string,var,inp_line,num_species) 

integer, intent(in) :: inp_line,num_species
character(len=*),intent(in) :: args(:,:)
character(len=*),intent(in) :: string
character(len=*),intent(out) :: var(:)

integer i, ios
character(len=40) temp
logical found

found = .false.

do i = 1,inp_line
    if(string.eq.trim(adjustl(args(1,i)))) then 
        do j = 1,(num_species-1)
            call split(args(2,i),',',temp)
            var(j) = trim(adjustl(temp))
!            call value(args(2,i),var(j),ios)
        enddo
        found = .true.
        var(num_species) = trim(adjustl(args(2,i)))
    elseif(i.eq.inp_line.and..not.found) then
        print *, 'Variable ',string,' not found. Check input file '

    endif
enddo

end subroutine variable_read_vec_string

end module input
