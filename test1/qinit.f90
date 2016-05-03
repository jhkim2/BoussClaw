
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
    
    ! Set flat state based on sea_level
    q = 0.d0
    forall(i=1:mx, j=1:my)
        q(1,i,j) = max(0.d0, sea_level - aux(1,i,j))
    end forall
    
    ! Add perturbation to initial conditions
    if (qinit_type > 0) then
        call add_perturbation(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    endif

    g = grav
    x1 = 0.d0
    y1 = 0.d0
    a0 = 1.d0
    a1 = 0.2d0
    k  = dsqrt(3.d0*a1/(a0+a1))/2.d0/a0
    c = dsqrt(g*(a0+a1))

    do j=1-mbc,my+mbc
        y=ylower + (j-0.5d0)*dy
        do i=1-mbc,mx+mbc
            x=xlower+ (i-0.5d0)*dx
            if (sqrt(x**2+y**2)<50.d0) then
            dist = abs(x+y)/sqrt(2.d0)
            rx= k*dist
            eta= a1*(1.d0/cosh(rx))**2.d0

            if (eta<1d-4) eta=0.d0
            q(1,i,j) = q(1,i,j) + eta
            if (eta>0.d0) then
               q(2,i,j) = q(1,i,j)*c*(1.d0-a0/q(1,i,j))/sqrt(2.d0)
               q(3,i,j) = q(1,i,j)*c*(1.d0-a0/q(1,i,j))/sqrt(2.d0)
            endif
            endif
        enddo
    enddo

    
end subroutine qinit
