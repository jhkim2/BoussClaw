
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use qinit_module, only: qinit_type,add_perturbation
    use geoclaw_module, only: sea_level,grav
    
    implicit none
    
    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Locals
    integer :: i,j,m,num_points,status
    real(kind=8) :: x,y,xim
    real(kind=8) :: x1,y1,a0,a1,k,c,eta,rx,g,dist
      ! Parameters for problem
    real(kind=8), parameter :: a =3.d0
    real(kind=8), parameter :: h0 =80.d0

    ! Other storage
    real(kind=8) :: omega,disc,xi,eta1

    ! Set flat state based on sea_level
    q = 0.d0
    forall(i=1:mx, j=1:my)
        q(1,i,j) = max(0.d0, sea_level - aux(1,i,j))
    end forall
    
    ! Add perturbation to initial conditions
    if (qinit_type > 0) then
        call add_perturbation(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    endif

  omega = sqrt((3.0d0/4.0d0)*(a/(h0**3)))
  !omega = 0.5*dsqrt((a/(h0**3))) 
  
  do i=1-mbc,mx+mbc
     x = xlower + (i - 0.5d0)*dx
     do j=1-mbc,my+mbc
        y = ylower + (j - 0.5d0) * dy

        eta=0.0d0
        eta1=0.0d0
        disc=0.0d0
        
        if (x >=6594 .and. x<9891) then
           xi=x-8243
           eta = a*(1.0d0/(cosh(omega*xi)))*(1.0d0/cosh(omega*xi))
           disc=eta*sqrt(grav/h0)

           q(1,i,j) = eta+q(1,i,j)
           q(2,i,j) = disc*q(1,i,j)
           q(3,i,j) = 0.0d0
        else
           q(1,i,j) = max(0.d0, sea_level-aux(1,i,j))
           q(2,i,j) = 0.0d0
           q(3,i,j) = 0.0d0
        endif
        
     enddo
  enddo

end subroutine qinit
