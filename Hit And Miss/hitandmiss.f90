program hitandmiss

  !Compilation:
  !$ make -f makehitandmiss

  !Running:
  !$ ./hitandmiss.exe

  !Calculating the volume of a N-sphere with Monte Carlo integration using
  !hit-and-miss method in dimensions 1-15. Number of generated points is
  !50000000. Output is the volume and statistical error of the result.
  
  ! Modules for Mersenne twister
  use mtmod
  
  implicit none
  
  integer, parameter :: rkk = selected_real_kind(10,40)
  integer :: npoints = 50000000, ndim, npointsinside, n, i

  ! v(:) will have all the coordinates values for the points generated
  ! in the hypercube
  real(kind=rkk), allocatable :: v(:) 
  real(kind=rkk) :: sq, V_e, p = 0.0, f_h, Integral, err 

  ! Initializing the random number generator
  call sgrnd(getseed(info=1))
  
  ! We will calculate the volume in dimensions 1...15, this main loop goes over
  ! them
  do ndim = 1,15
     allocate(v(ndim))

     ! Calculating the volume by simulating poins in a cube
     
     do n = 1,npoints
        do i = 1,ndim
           v(i) = 3*grnd()-1.5 ! RNG's generated from [0,1[, that's why big cube
        end do
        ! Evaluate if the point lands inside the sphere
        sq = DOT_PRODUCT(v, v)
        if (sq<1) then
           p = p+1    ! p tracks how many points land inside the ball
        end if
     end do
     
     V_e = 3**ndim    ! Side of cube is 3
     f_h = p/npoints  ! Fraction of points that hit the ball
     Integral = V_e*f_h      ! Result of the integral
     err = V_e*(sqrt((f_h-f_h**2)/npoints)) ! Error of the average

     ! Print results on screen
     print '("V is", x, f10.7, x, "with error", x, f10.8, x, "for dim", x, i2)'&
          &,Integral, err, ndim
     
     ! Reset variables
     deallocate(v)
     p = 0
     
  end do

end program hitandmiss
     

     
        
        
        
