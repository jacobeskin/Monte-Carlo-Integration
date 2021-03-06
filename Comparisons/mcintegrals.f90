program mcintegrals
! Modules for Mesenne twister
  use mtmod
  
  implicit none

  
  integer, parameter :: rkk = selected_real_kind(10,40)
  integer :: npoints, n, j, k, p

  real(kind=rkk) :: x, y, f, fsum, i_ds, i_ss, i_hm, i_pss, i_is, box
  real(kind=rkk) :: i_ex = 0.7034804524336, t_ds, t_ss, t_hm, t_is, start, end
  real(kind=rkk) :: di_ds, di_ss, di_hm, di_pss, di_is, V, t_pss

  ! Initializing the random number generator
  call sgrnd(getseed(info=1))

  print *,

  ! Calculations are done with npoints from 10^2 up to 10^6
  do n = 2,6
     npoints = 10**n
     
!----- Direct sampling -----
     call cpu_time(start)
     fsum = 0.0
     do j = 1,npoints 
        x = grnd()
        if(x==0.0)then
           x = grnd()
        end if
        f = 1.0/(x**(1.0/3.0)+x**(1.0/4.0))
        fsum = fsum+f
     end do
     i_ds = fsum/npoints
     di_ds = i_ds-i_ex
     call cpu_time(end)
     t_ds = end-start

!----- Stratied sampling -----
     call cpu_time(start)
     fsum = 0.0
     do j = 1,npoints
        x = (j+grnd())*(1.0/npoints)
        if(x==0.0)then
           x = (j+grnd())*(1.0/npoints)
        end if
        f = 1.0/(x**(1.0/3.0)+x**(1.0/4.0))
        fsum = fsum+f
     end do
     i_ss = fsum/npoints
     di_ss = i_ss-i_ex
     call cpu_time(end)
     t_ss = end-start

!----- Hit & miss sampling -----
     call cpu_time(start)
     V = huge(x)
     p = 0
     do j = 1,npoints
        x = grnd()
        if(x==0.0)then
           x = grnd()
        end if
        y = grnd()*huge(x)
        f = 1.0/(x**(1.0/3.0)+x**(1.0/4.0))
        if(y<f)then
           p = p+1
        end if
     end do
     i_hm = V*(p/npoints)
     di_hm = i_hm-i_ex
     call cpu_time(end)
     t_hm = end-start

!----- Partially stratied sampling -----
     call cpu_time(start)
     fsum = 0.0
     do k = 1,(npoints/100)
        do j = 1,100
           x = (j+grnd())*(1.0/100.0)
           if(x==0.0)then
              x = (j+grnd())*(1.0/100.0)
           end if
           f = 1.0/(x**(1.0/3.0)+x**(1.0/4.0))
           fsum = fsum+f
        end do
     end do
     i_pss = fsum/npoints
     di_pss = i_pss-i_ex
     call cpu_time(end)
     t_pss = end-start

!----- Importance sampling -----

     call cpu_time(start)
     fsum = 0.0
     do j = 1,npoints
        x = (5.0*log(2.0)*grnd()/46.0)**(3.0/5.0)
        f = (12.0/log(2.0))*(x**(2.0/3.0))
        fsum = fsum+f
     end do
     i_is = fsum/npoints
     di_is = i_is-i_ex
     call cpu_time(end)
     t_is = end-start

!----- Printing results -----
     print '("For N=", i7)', npoints
     print '("DS:", x, "DI=", x, f14.12, x, "Time=", f11.9)', di_ds, t_ds
     print '("SS:", x, "DI=", x, f14.12, x, "Time=", f11.9)', di_ss, t_ss
     print '("HM:", x, "DI=", x, f14.12, x, "Time=", f11.9)', di_hm, t_hm
     print '("PSS:", x, "DI=", x, f14.12, x, "Time=", f11.9)', di_pss, t_pss
     print '("IS:", x, "DI=", x, f15.12, x, "Time=", f11.9)', di_is, t_is
     print *,
     
  end do

end program mcintegrals


           
           
           
     
     
