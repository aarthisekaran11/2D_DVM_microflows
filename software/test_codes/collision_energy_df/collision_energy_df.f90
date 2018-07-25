program collision_energy_df

  implicit none
  
  double precision, parameter :: pi = 3.14159265358979d0
  double precision, parameter :: double_tol = 1.0d-12
  double precision, parameter :: omega = 1.0d0

  integer, parameter :: num_points_x = 35
  integer, parameter :: num_points_y = 35
  integer, parameter :: num_points_z = 35

  double precision, parameter :: xv_min = -3.5d0
  double precision, parameter :: xv_max =  3.5d0
  double precision, parameter :: yv_min = -3.5d0
  double precision, parameter :: yv_max =  3.5d0
  double precision, parameter :: zv_min = -3.5d0
  double precision, parameter :: zv_max =  3.5d0

  double precision :: T
  double precision, parameter :: n = 1.0d0
  double precision, parameter :: m = 1.0d0

  double precision, dimension(1:num_points_x,1:num_points_y,1:num_points_z) :: phi
  double precision, dimension(1:num_points_x) :: x
  double precision, dimension(1:num_points_y) :: y
  double precision, dimension(1:num_points_z) :: z

  double precision, dimension(:), allocatable :: g_full
  double precision, dimension(:), allocatable :: g
  double precision, dimension(:), allocatable :: g_df

  integer :: num_points

  integer :: i_min, i_max
  integer :: j_min, j_max
  integer :: k_min, k_max

  integer :: i, j, k, ii, jj, kk
  integer :: l, ll, p

  double precision :: beta_x, beta_y, beta_z, beta3
  double precision :: x_min, y_min, z_min

  double precision :: c_x, c_y, c_z, c2
  double precision :: coeff

  integer :: count, g_length

  double precision :: g_x, g_y, g_z, gr
  double precision :: m_r, Et
  double precision :: sum1, sum2
  double precision :: avg1, avg2

  double precision :: temp
  integer :: iswap(1), itemp, iswap1

  integer :: file_unit
  character(len=128) :: filename

  integer :: status

  double precision :: cdf1, cdf2

  write(*,*) "Temperature: "
  read(*,*)  T

  write(*,*) "double_tol = ",double_tol

  write(*,*) "started..."

  file_unit = 9
  filename = "collision_energy_df2.dat"
  open( unit=file_unit, file=filename, status='unknown' )

  write( file_unit, fmt="(a)" ) 'title = "Collision Probabilities"'
  write( file_unit, fmt="(a)" ) 'variables = "g", "Et", "P(Et)", "B(Et)", "P(Ec)", "B(Ec)", &
       "cdf1", "cdf2"'

  num_points = num_points_x*num_points_y*num_points_z
!!$  num_velocities = int( ( num_points + 1 )*( num_points )/2 )

  i_min = 1
  i_max = num_points_x
  j_min = 1
  j_max = num_points_y
  k_min = 1
  k_max = num_points_z

  beta_x = (xv_max - xv_min)/dble(num_points_x-1)
  beta_y = (yv_max - yv_min)/dble(num_points_y-1)
  beta_z = (zv_max - zv_min)/dble(num_points_z-1)
  beta3  = beta_x*beta_y*beta_z

  ! velocities
  do i = i_min, i_max
     x(i) = xv_min + dble(i)*beta_x
  end do

  do j = j_min, j_max
     y(j) = yv_min + dble(j)*beta_y
  end do

  do k = k_min, k_max
     z(k) = zv_min + dble(k)*beta_z
  end do

  allocate( g_full(1:num_points) )
  
  ! calculate Maxwellian and make relative velocity list
  coeff = n*sqrt( m*m*m )/sqrt( pi*pi*pi*T*T*T )

  x_min = x(i_min)
  y_min = y(j_min)
  z_min = z(k_min)

  ! Calculate phi and relative velocities
  do l = 1, num_points

     k = k_min + (l-1)/(num_points_x*num_points_y)
     j = j_min + (l-1)/num_points_x - (k - k_min)*num_points_y
     i = i_min + (l-1) - (j - j_min)*num_points_x - (k - k_min)*num_points_x*num_points_y

     c_z = z(k)*z(k)
     c_y = y(j)*y(j)
     c_x = x(i)*x(i)
     c2 = c_x + c_y + c_z

     phi(i,j,k) = coeff*exp(-c2*m/T)*beta3

     g_z = ( z(k) - z_min )*( z(k) - z_min )
     g_y = ( y(j) - y_min )*( y(j) - y_min )
     g_x = ( x(i) - x_min )*( x(i) - x_min )   
     gr = sqrt( g_x + g_y + g_z )

     g_full(l) = gr

  end do

  ! Order g_full vector
  do i = 1, num_points - 1

     iswap = minloc( g_full(i:num_points) )
     iswap1 = iswap(1) + i - 1

     if( iswap1 .ne. i )then
        temp = g_full(i)
        g_full(i) = g_full(iswap1)
        g_full(iswap1) = temp
     end if

  end do

  ! Count number of unique values
  count = 1
  do i = 2, num_points
     if( abs( g_full(i) - g_full(i-1) ) .lt. double_tol )then
        cycle
     else
        count = count + 1
     end if
  end do
  g_length = count

  write(*,*) "found ",g_length," number of unique relative velocities"

  allocate( g(1:g_length), g_df(1:g_length) )

  ! Set g vector
  count = 1
  g(count) = g_full(count)
  do i = 2, num_points
     if( abs( g_full(i) - g_full(i-1) ) .lt. double_tol )then
        cycle
     else
        count = count + 1
        g(count) = g_full(i)
     end if
  end do

  ! Calculate g_df
  g_df = 0.0d0
  sum1 = 0.0d0

  do l = 1, num_points

     k = k_min + (l-1)/(num_points_x*num_points_y)
     j = j_min + (l-1)/num_points_x - (k - k_min)*num_points_y
     i = i_min + (l-1) - (j - j_min)*num_points_x - (k - k_min)*num_points_x*num_points_y

     do ll = l, num_points

        kk = k_min + (ll-1)/(num_points_x*num_points_y)
        jj = j_min + (ll-1)/num_points_x - (kk - k_min)*num_points_y
        ii = i_min + (ll-1) - (jj - j_min)*num_points_x - (kk - k_min)*num_points_x*num_points_y
        
        g_z = ( z(k) - z(kk) )*( z(k) - z(kk) )
        g_y = ( y(j) - y(jj) )*( y(j) - y(jj) )
        g_x = ( x(i) - x(ii) )*( x(i) - x(ii) )

        gr = sqrt( g_x + g_y + g_z )

        do p = 1, g_length
           if( abs( gr*gr - g(p)*g(p) ) .lt. double_tol )then
              g_df(p) = g_df(p) + phi(i,j,k) * phi(ii,jj,kk) * gr**(2.0d0-2.0d0*omega)
              sum1 = sum1 + phi(i,j,k) * phi(ii,jj,kk) * gr**(2.0d0-2.0d0*omega)
              exit
           end if
           if( p .eq. g_length )then
              write(*,*) "relative velocity not found. g = ",gr
              write(*,*) g_x, g_y, g_z
              write(*,*) i, j, k, ii, jj, kk
              stop
           end if
        end do

     end do
  end do

  g_df = g_df/sum1

  write(*,*) "completed calculation of g_df"

  sum2 = 0.0d0
  do i = 1, g_length
     Et = 0.25d0*m*g(i)*g(i)
!!$     if( i .gt. 1 .and. i .lt. g_length )then
!!$        gr = 0.25d0*m*abs( g(i-1)*g(i-1) - g(i+1)*g(i+1) )*0.5d0
!!$     else if( i .lt. 1 )then
!!$        gr = 0.25d0*m*abs( g(i)*g(i) - g(i+1)*g(i+1) )*0.5d0
!!$     else
!!$        gr = 0.25d0*m*abs( g(i)*g(i) - g(i-1)*g(i-1) )*0.5d0
!!$     end if

     sum2 = sum2 + Et**(1.5-omega)*exp(-2.0d0*Et/T)!*gr
  end do

  avg1 = 0.0d0
  avg2 = 0.0d0

  cdf1 = 0.0d0
  cdf2 = 0.0d0

  write( file_unit, fmt="(a7,i6)" ) "zone i=",g_length
  do i = 1, g_length
     Et = 0.25d0*m*g(i)*g(i)
     
     cdf1 = cdf1 + g_df(i)
     cdf2 = cdf2 + Et**(1.5-omega)*exp(-2.0d0*Et/T)/sum2

     write( file_unit, fmt="(8e20.12)" ) &
          g(i), Et, &
          g_df(i), Et**(1.5-omega)*exp(-2.0d0*Et/T)/sum2, &
          g_df(i)/exp(-2.0d0*Et/T), Et**(1.5-omega)/sum2, &
          cdf1, cdf2

     avg1 = avg1 + g_df(i)*Et
     avg2 = avg2 + Et**(1.5-omega)*exp(-2.0d0*Et/T)*Et

  end do

  avg1 = avg1
  avg2 = avg2/sum2

  write(*,*) "dvm formulation average = ",avg1
  write(*,*) "LB formulation average = ",avg2
  write(*,*) "(5/2-w)*T/2 = ",(2.5d0-omega)*T*0.5d0

  deallocate( g_full, g, g_df )

  close( file_unit )

end program collision_energy_df
