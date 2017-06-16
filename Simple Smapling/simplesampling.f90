program simplesampling

  !Compilation:
  !$ make -f makesimplesampling

  !Running:
  !$ ./simplesampling.exe

  !Evaluating the volume of an N-dimensional hypersphere (radius 1) with
  !Monte Carlo integration using simple sampling. Integration is done for
  !dimension 1 to 15, and as an output the program prints the volume and
  !statistical error. Number of points generated for each integration is
  !50000000. 
  
  ! Modules for Mesenne twister
  use mtmod
  
  implicit none
  
  integer, parameter :: rkk = selected_real_kind(10,40)
  integer :: npoints = 50000000, ndim, npointsinside, n, i, p=0

  ! v(:) will have all the coordinates values for the points generated
  ! w(:) will have the volume terms V(m) of the formula V(m)<f>
  real(kind=rkk), allocatable :: v(:), w(:) 
  real(kind=rkk) :: sq, pi, f, ff, fsum = 0.0, ffsum = 0.0, fm, ffm, Intgrl, err

  ! Initializing the random number generator
  call sgrnd(getseed(info=1))

  ! Defining pi
  pi = acos(-1.0) 

  ! Allocate w(:)
  allocate(w(15))
  w(1) = 2
  w(2) = pi
  w(3) = 4*pi/3
  w(4) = (pi**2)/2
  w(5) = (8*pi**2)/15
  w(6) = (pi**3)/6
  w(7) = (16*pi**3)/105
  w(8) = (pi**4)/24
  w(9) = (32*pi**4)/945
  w(10) = (pi**5)/120
  w(11) = (64*pi**5)/10395
  w(12) = (pi**6)/720
  w(13) = (128*pi**6)/135135
  w(14) = (pi**7)/5040
  w(15) = (256*pi**7)/2027025
  
  ! We will calculate the volume in dimensions 1...15, this main loop goes over
  ! them
  do ndim = 1,15
     allocate(v(ndim))

     ! Calculating the volume by simulating points

     do n = 1,npoints
        do i = 1,ndim
           v(i) = 3*grnd()-1.5 ! RNG's generated from [0,1[
        end do
        ! Evaluate if the point lands inside the sphere
        sq = DOT_PRODUCT(v, v)
        if (sq<1) then 
           p = p+1
           f = sqrt(1-sq)
           fsum = fsum+f 
           ff = 1-sq
           ffsum = ffsum+ff
        end if
     end do

     fm = fsum/p
     ffm = ffsum/p
     Intgrl = 2*w(ndim)*fm
     
     err = w(ndim)*sqrt((ffm-fm**2)/npoints) ! Error of the average
     
     ! Print results on screen
     print '("V is", x, f10.7, x, "with error", x, f10.8, x, "for dim", x, i2)'&
          &,Intgrl, err, ndim
     
     ! Reset variables
     deallocate(v)
     fsum = 0.0
     ffsum = 0.0
     p = 0

  end do

end program simplesampling
