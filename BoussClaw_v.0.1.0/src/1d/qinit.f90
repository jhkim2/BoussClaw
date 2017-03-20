
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use qinit_module, only: qinit_type,add_perturbation
    use geoclaw_module, only: sea_level,grav,dry_tolerance
    
    implicit none
    
    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Locals
    integer :: i,j,m
    real(kind=8) :: x1,a0,a1,k,g,c,rx,x,y,eta

    g=grav
    
    ! Set flat state based on sea_level
    q = 0.d0

    x1 = 24.60338d0
    a0 = 1.d0
    a1 = 0.2d0
    k  = dsqrt(3.d0*a1/(a0+a1))/2.d0/a0
    c = dsqrt(g*(a0+a1))
    x1 = x1 + 6.d0/sqrt(g)*c

    do j=1-mbc,my+mbc
        do i=1-mbc,mx+mbc
            q(1,i,j) = max(0.d0, sea_level - aux(1,i,j))
        enddo
    enddo

    do j=1,my
        y=ylower + (j-0.5d0)*dy
        do i=1,mx
            x=xlower+ (i-0.5d0)*dx
            rx= k*(x-x1)
            eta= a1*(1.d0/cosh(rx))**2.d0
            if (eta<1d-4) eta=0.d0

            !q(1,i,j) = q(1,i,j) + eta
            if (eta>0.d0) then
            !   q(2,i,j) = -q(1,i,j)*c*(1.d0-a0/q(1,i,j))
            endif
            q(3,i,j)=0.d0

        enddo
    enddo

    do j=1,my
        y=ylower + (j-0.5d0)*dy
        do i=1,mx
            x=xlower+ (i-0.5d0)*dx
            rx= k*(x-x1)
            eta= a0*a1/(a0+a1)*(1.d0/cosh(rx))**2.d0 &
                +a1**2/(a0+a1)*(1.d0/cosh(rx))**4.d0
            if (eta<1d-8) eta=0.d0

            q(1,i,j) = q(1,i,j) + eta
            if (eta>0.d0) then
               q(2,i,j) = -q(1,i,j)*g*a1/c*(1.d0/cosh(rx))**2.d0
            endif
            q(3,i,j)=0.d0

        enddo
    enddo
  
    if (.false.) then
        open(23, file='fort.aux',status='unknown',form='formatted')
        print *,'Writing out aux arrays'
        print *,' '
        do j=1,my
            do i=1,mx
                write(23,*) i,j,(q(m,i,j),m=1,meqn)
            enddo
        enddo
        close(23)
    endif
    
end subroutine qinit
