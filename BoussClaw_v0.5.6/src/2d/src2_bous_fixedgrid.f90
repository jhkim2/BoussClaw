
  subroutine src2_bous(meqn,mbc,mx,my, &
                  xlower,ylower,dx,dy,q,maux,aux,t,dt)

      use geoclaw_module, only: g => grav, coriolis_forcing, coriolis
      use geoclaw_module, only: friction_forcing, friction_depth
      use geoclaw_module, only: manning_coefficient,sea_level
      use bous_module

      implicit none
            
      integer(kind=4) i,j,mx,my,k,k1,k2,mxy,JOB
      integer(kind=4), parameter :: ndim = 2
      
      real(kind=8)   q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8)  q0(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) psi(1:2*(mx+2)*(my+2))
      
      real(kind=8) slope(1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) hbar_tt(1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) hbar_xtt(1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) hbar_ytt(1-mbc:mx+mbc,1-mbc:my+mbc)
      
      real(kind=8)   xx(2*(mx+2)*(my+2))
      real(kind=8)   xx0(4,2*(mx+2)*(my+2))
     
      integer(kind=4), allocatable, dimension (:) :: IA
      integer(kind=4), allocatable, dimension (:) :: JA
      real(kind = 8 ), allocatable, dimension (:) :: SA
            
      real(kind=8) :: delt,ERR,TOL2

      real(kind=8) tol,w,t,dt,s,x,xc,xi0,y,dt1
      real(kind=8) dx,dy,xlower,ylower,hx,hy
      integer(kind=4) INFO,maux,mbc,meqn
     
      integer(kind=4), parameter :: NSAVE=10
      integer(kind=4) :: LENW
      real(kind=8), allocatable, dimension(:) :: RWORK
      integer(kind=4) :: LENIW
      integer(kind=4), allocatable, dimension(:) :: IWORK            

      integer(kind=4) :: IERR,ISYM,ITER
      integer(kind=4) :: ITMAX,ITOL,IUNIT,NELT,m,ii,nstep
      
      integer ( kind = 4 ), parameter :: ligw = 20
      integer ( kind = 4 ), parameter :: maxl = 20

      integer ( kind = 4 ) lrgw
      
      integer ( kind = 4 ) igwk(ligw)
      real ( kind = 8 ), allocatable ::rgwk(:)
      real ( kind = 8 ), allocatable :: sb(:)
      real ( kind = 8 ), allocatable :: sx(:)
      
      !EXTERNAL MATVEC_TRIAD, MSOLVE_IDENTITY
      
      logical verbose

      verbose=.False.
           
!     ----------------------------------------------------------------
!
!     Boussinesq type-------------------------------------------------

      mxy=(mx+2)*(my+2)
     
      nstep=ceiling(max(dt/dx,dt/dy)*2*sqrt(2.))

      hbar_xtt=0.d0
      hbar_ytt=0.d0
      hbar_tt=0.d0
      slope(:,:) = aux(1,:,:)-sea_level

      delt=dt/nstep

      if (t == 0.d0 ) then

      	allocate (A_inverse(1:2*mxy,1:2*mxy))

     	JOB=0
        NELT=1
         
        call read_A(mx,my,meqn,mbc,dx,dy, &
             q,maux,aux,IA,JA,SA,NELT,slope(:,:),JOB,hbar_xtt,hbar_ytt,ndim)

        allocate ( IA(1:nelt) )
        allocate ( JA(1:nelt) )
        allocate ( SA(1:nelt) )
     
        JOB=1
     
        call read_A(mx,my,meqn,mbc,dx,dy, &
              q,maux,aux,IA,JA,SA,NELT,slope(:,:),JOB,hbar_xtt,hbar_ytt,ndim)

        ISYM = 0
        TOL2 = 1d-6
        TOL  = TOL2
        IUNIT= 0
        ITMAX= 50000

        ITOL = 0
        LRGW = 100 + 2*mxy* ( maxl + 6 ) + maxl * ( maxl + 3 )
        !LRGW = 78732417
        
        allocate(RWORK(1))
        allocate(IWORK(1))
        allocate(RGWK(1:LRGW))
        allocate(SB(1:2*mxy))
        allocate(SX(1:2*mxy))
        
        igwk(1) = maxl
        igwk(2) = maxl
        igwk(3) = 0
        igwk(4) = 0
        igwk(5) = 60

        A_inverse = 0.d0

        do i = 1,2*mxy
        	A_inverse(i,i) = 1.d0
        enddo

        CALL DGMRES (2*mxy, A_inverse, A_inverse, NELT, IA, JA, SA, ISYM, MATVEC_TRIAD, &
                     MSOLVE_IDENTITY, ITOL, TOL, ITMAX, ITER, ERR, IERR, &
                     IUNIT, SB, SX, RGWK, LRGW, &
                     IGWK, LIGW, RWORK, IWORK)

        deallocate(SB)
        deallocate(SX)
        deallocate(RGWK)
      
        deallocate(RWORK)
        deallocate(IWORK)
        
        deallocate ( IA )
        deallocate ( JA )
        deallocate ( SA )

    endif
    
    do ii=1,nstep
        
        ! RK4   
        
        ! First Stage
        
        call read_psi(mx,my,meqn,mbc, &
        dx,dy,q,maux,aux,psi,g,slope(:,:) &
        ,hbar_xtt,hbar_ytt)
        
        XX(1:2*mxy) = MATMUL(A_inverse, psi)

        XX0(1,1:2*mxy) = XX(1:2*mxy)
        
        do j=1,my
            do i=1,mx
                k1=i+j*(mx+2)+1
                k2=mxy+k1
                q0(2,i,j)=q(2,i,j)-dt/2.d0*xx(k1)
                q0(3,i,j)=q(3,i,j)-dt/2.d0*xx(k2)
            end do
        end do
        
        ! Second Stage
        
        call read_psi(mx,my,meqn,mbc,dx,dy,q0,maux,aux,psi, &
        g,slope(:,:),hbar_xtt,hbar_ytt)

        XX(1:2*mxy) = MATMUL(A_inverse, psi)

        XX0(2,1:2*mxy) = XX(1:2*mxy)
        
        do j=1,my
            do i=1,mx
                k1=i+j*(mx+2)+1
                k2=mxy+k1
                q0(2,i,j)=q(2,i,j)-dt/2.d0*xx(k1)
                q0(3,i,j)=q(3,i,j)-dt/2.d0*xx(k2)
            end do
        end do
        
        ! Third Stage
        
        call read_psi(mx,my,meqn,mbc,dx,dy,q0,maux,aux,psi, &
        g,slope(:,:),hbar_xtt,hbar_ytt)
        
        XX(1:2*mxy) = MATMUL(A_inverse, psi)

        XX0(3,1:2*mxy) = XX(1:2*mxy)

        do j=1,my
            do i=1,mx
                k1=i+j*(mx+2)+1
                k2=mxy+k1
                q0(2,i,j)=q(2,i,j)-dt*xx(k1)
                q0(3,i,j)=q(3,i,j)-dt*xx(k2)
            end do
        end do
        
        ! Fourth Stage
        
        call read_psi(mx,my,meqn,mbc,dx,dy,q0,maux,aux,psi, &
          g,slope(:,:),hbar_xtt,hbar_ytt)
                
        XX(1:2*mxy) = MATMUL(A_inverse, psi)

        XX0(4,1:2*mxy) = XX(1:2*mxy)

!=======================================================================

        do j=1,my
          do i=1,mx
            k1=i+j*(mx+2)+1
            k2=mxy+k1
                              
            q(2,i,j)=q(2,i,j)-delt/6.d0*( xx0(1,k1) + 2.d0*xx0(2,k1) &
               +2.d0*xx0(3,k1) + xx0(4,k1))
            q(3,i,j)=q(3,i,j)-delt/6.d0*( xx0(1,k2) + 2.d0*xx0(2,k2) &
               +2.d0*xx0(3,k2) + xx0(4,k2))

          enddo
        enddo
     
      enddo
      
 999  continue

      aux(5,:,:)=aux(4,:,:)
      aux(4,:,:)=aux(1,:,:)
      aux(4,-1,-1)=dt
      
  end subroutine src2_bous
  
! =========================================================
      subroutine read_A(mx,my,meqn,mbc,dx,dy, &
     &      q,maux,aux,IA,JA,SA,NELT,slope,JOB,hbar_xtt,hbar_ytt,ndim)
! =========================================================
      use bous_module

      implicit none
     
      integer(kind=4) mx,my,meqn,mbc,maux,NELT,mxy,JOB
      integer(kind=4) i,j,k,k1,NCHK,ndim
      real(kind=8) dx,dy

      real(kind=8)   q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) hbar_xtt(1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) hbar_ytt(1-mbc:mx+mbc,1-mbc:my+mbc)
      
      integer(kind=4) IA(1:NELT)
      integer(kind=4) JA(1:NELT)
      real(kind = 8 ) SA(1:NELT)
      
      real(kind=8)  HH(1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) HHH(1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) slope(1-mbc:mx+mbc,1-mbc:my+mbc)
      
      mxy=(mx+2)*(my+2)
      
      if (JOB == 0 ) then
         NELT=0
        do i=0,mx+1
           do j=0,my+1
               NELT=NELT+2
           enddo
        enddo
      
         do j=1,my
    
           do i=1,mx

              if (q(1,i,j).lt.1.d-1.or.slope(i,j)>-1.d-1) then

              ! do nothing
        
              elseif (maxval(slope(i-2:i+2,j-2:j+2))  > sw_depth) then

              ! do nothing

              elseif (maxval(abs(hbar_xtt(i-2:i+2,j-2:j+2)))>1.d-14 .or. &
                   maxval(abs(hbar_ytt(i-2:i+2,j-2:j+2)))>1.d-14    ) then

              ! do nothing

              elseif (maxval(q(1,i-2:i+2,j-2:j+2)) < 1.d-1) then

              ! do nothing

              elseif (maxval(slope(i-1:i+1,j-1:j+1)) <= sw_depth ) then              

      ! D1 part     

              NELT=NELT+3
     
     ! D2 part
     
              NELT=NELT+4        
                  
     ! D3 part
        
              NELT=NELT+4

     ! D4 part
         
              NELT=NELT+3  
        
               endif
            
 997    continue
            
           enddo
        enddo
     
     elseif (JOB==1) then
     
         NCHK=0
              
         do j=1-mbc,my+mbc
            do i=1-mbc,mx+mbc
               HH(i,j)=(max(0.,-slope(i,j)))**2
              HHH(i,j)=(max(0.,-slope(i,j)))**3
            enddo
         enddo
    
        do i=0,mx+1
           do j=0,my+1
              k=i+j*(mx+2)+1
              
              NCHK = NCHK+1
              IA(NCHK)=k
              JA(NCHK)=k
              SA(NCHK)=1.d0

              NCHK = NCHK+1
              IA(NCHK)=k+mxy
              JA(NCHK)=k+mxy
              SA(NCHK)=1.d0
           enddo
        enddo
      
         do j=1,my
    
           do i=1,mx
        
              if (q(1,i,j).lt.1.d-1.or.slope(i,j)>-1.d-1) then

              ! do nothing
        
              elseif (maxval(slope(i-2:i+2,j-2:j+2))  > sw_depth) then

              ! do nothing

              elseif (maxval(abs(hbar_xtt(i-2:i+2,j-2:j+2)))>1.d-14 .or. &
                   maxval(abs(hbar_ytt(i-2:i+2,j-2:j+2)))>1.d-14    ) then

              ! do nothing

              elseif (maxval(q(1,i-2:i+2,j-2:j+2)) < 1.d-1) then

              ! do nothing

              elseif (maxval(slope(i-1:i+1,j-1:j+1)) <= sw_depth ) then              

      ! D1 part     

           k= j*(mx+2)+i+1
           
              NCHK = NCHK+1
              IA(NCHK)=k
              JA(NCHK)=k
              SA(NCHK)=2.d0*(B_param+.5d0)*HH(i,j)/dx**2 &
     &        -2.d0/6.d0*HHH(i,j)/(-slope(i,j))/dx**2
     
              NCHK = NCHK+1
              IA(NCHK)=k
              JA(NCHK)=k+1
              SA(NCHK)=-(B_param+.5d0)*HH(i,j)/dx**2 &
     &        +1.d0/6.d0*HHH(i,j)/(-slope(i+1,j))/dx**2
     
              NCHK = NCHK+1
              IA(NCHK)=k
              JA(NCHK)=k-1
              SA(NCHK)=-(B_param+.5d0)*HH(i,j)/dx**2 &
     &        +1.d0/6.d0*HHH(i,j)/(-slope(i-1,j))/dx**2
    
     
     ! D2 part
     
           k1= mxy+k
           
              NCHK = NCHK+1
              IA(NCHK)=k
              JA(NCHK)=k1-(mx+2)+1
              SA(NCHK)=(B_param+.5)*hh(i,j)/4.d0/dx/dy &
     &       -1.d0/6.d0*HHH(i,j)/(-slope(i+1,j-1))/4.d0/dx/dy

              NCHK = NCHK+1
              IA(NCHK)=k
              JA(NCHK)=k1-(mx+2)-1
              SA(NCHK)=-(B_param+.5)*hh(i,j)/4.d0/dx/dy &
     &       +1.d0/6.d0*HHH(i,j)/(-slope(i-1,j-1))/4.d0/dx/dy

              NCHK = NCHK+1
              IA(NCHK)=k
              JA(NCHK)=k1+(mx+2)+1
              SA(NCHK)=-(B_param+.5)*hh(i,j)/4.d0/dx/dy &
     &       +1.d0/6.d0*HHH(i,j)/(-slope(i+1,j+1))/4.d0/dx/dy

              NCHK = NCHK+1
              IA(NCHK)=k
              JA(NCHK)=k1+(mx+2)-1
              SA(NCHK)= (B_param+.5)*hh(i,j)/4.d0/dx/dy &
     &       -1.d0/6.d0*HHH(i,j)/(-slope(i-1,j+1))/4.d0/dx/dy           
                  
     ! D3 part
        
              k= j*(mx+2)+i+mxy+1
              k1=k-mxy
           
              NCHK = NCHK+1
              IA(NCHK)=k
              JA(NCHK)=k1-(mx+2)+1
              SA(NCHK)= (B_param+.5)*hh(i,j)/4.d0/dx/dy &
     &       -1.d0/6.d0*HHH(i,j)/(-slope(i+1,j-1))/4.d0/dx/dy

              NCHK = NCHK+1
              IA(NCHK)=k
              JA(NCHK)=k1-(mx+2)-1
              SA(NCHK)=-(B_param+.5)*hh(i,j)/4.d0/dx/dy &
     &       +1.d0/6.d0*HHH(i,j)/(-slope(i-1,j-1))/4.d0/dx/dy

              NCHK = NCHK+1
              IA(NCHK)=k
              JA(NCHK)=k1+(mx+2)+1
              SA(NCHK)=-(B_param+.5)*hh(i,j)/4.d0/dx/dy &
     &       +1.d0/6.d0*HHH(i,j)/(-slope(i+1,j+1))/4.d0/dx/dy

              NCHK = NCHK+1
              IA(NCHK)=k
              JA(NCHK)=k1+(mx+2)-1
              SA(NCHK)= (B_param+.5)*hh(i,j)/4.d0/dx/dy &
     &       -1.d0/6.d0*HHH(i,j)/(-slope(i-1,j+1))/4.d0/dx/dy

     ! D4 part
         
              k=mxy+j*(mx+2)+i+1
           
              NCHK = NCHK+1
              IA(NCHK)=k
              JA(NCHK)=k
              SA(NCHK)= 2.d0*(B_param+.5d0)*HH(i,j)/dy**2 &
     &        -2.d0*HHH(i,j)/6.d0/(-slope(i,j))/dy**2

              NCHK = NCHK+1
              IA(NCHK)=k
              JA(NCHK)=k+mx+2
              SA(NCHK)=-(B_param+.5d0)*HH(i,j)/dy**2 &
     &        +HHH(i,j)/6.d0/(-slope(i,j+1))/dy**2

              NCHK = NCHK+1
              IA(NCHK)=k
              JA(NCHK)=k-(mx+2)
              SA(NCHK)=-(B_param+.5d0)*HH(i,j)/dy**2 &
     &        +HHH(i,j)/6.d0/(-slope(i,j-1))/dy**2  
        
            endif
         
998     continue
            
           enddo
        enddo
               
        if (NCHK.ne.NELT) then
           print*,'ERROR NELT is not correct'
           print*,'NCHK=',NCHK,'NELT=',NELT
           stop
        endif
     
     else
        print*,'Error with read_A in boud_module'
        stop
      endif
             
      return
      end
   
!======================================================================
      subroutine read_psi(mx,my,meqn,mbc, &
     &      dx,dy,q,maux,aux,psi,g,slope,hbar_xtt,hbar_ytt)

      use geoclaw_module, only: sea_level
      use bous_module
     
      implicit none

      integer(kind=4) mx,my,meqn,mbc,maux,mxy,i,j,k,k2,m,ndim,kk
      real   (kind=8) dx,dy,g
     
      real(kind=8)   q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) slope(1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) hbar_xtt(1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) hbar_ytt(1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) psi(1:2*(mx+2)*(my+2))
      real(kind=8)   hh(1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8)  hhh(1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8)  hetax(-1:mx+2,-1:my+2)
      real(kind=8)  hetay(-1:mx+2,-1:my+2)
      real(kind=8)  eta(1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8)  huv(1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8)  hu2(1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8)  hv2(1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8)  detax(-1:mx+2,-1:my+2)
      real(kind=8)  detay(-1:mx+2,-1:my+2)
      real(kind=8)  s1(0:mx+1,0:my+1)
      real(kind=8)  s2(0:mx+1,0:my+1)
      real(kind=8)  s1_H(0:mx+1,0:my+1)
      real(kind=8)  s2_H(0:mx+1,0:my+1)
      real(kind=8)  tol,topo
     
     logical swe
     
     swe=.false.

      mxy=(mx+2)*(my+2)
      
      do j=1-mbc,my+mbc
         do i=1-mbc,mx+mbc
           if (q(1,i,j).gt.1d-4) then
              huv(i,j)= q(2,i,j)*q(3,i,j)/q(1,i,j)
              hu2(i,j)= q(2,i,j)**2/q(1,i,j)
              hv2(i,j)= q(3,i,j)**2/q(1,i,j)
           else
              huv(i,j)= 0.d0
              hu2(i,j)= 0.d0
              hv2(i,j)= 0.d0
           endif
           eta(i,j)= q(1,i,j) + slope(i,j) 
            HH(i,j)=(max(0.d0,-slope(i,j)))**2
           HHH(i,j)=(max(0.d0,-slope(i,j)))**3
        enddo
     enddo
     
     hetax = 0.d0
     hetay = 0.d0
     detax = 0.d0
     detay = 0.d0

      do j= 1,my
         do i= 1,mx
         if (slope(i+1,j)<0.d0.and.slope(i-1,j)<0.d0.and.slope(i,j)<0.d0) then
            if (i==1) then
                hetax(i,j)= q(1,i,j)*(eta(i+1,j)-eta(i,j))/dx
                detax(i,j)=-slope(i,j)*(eta(i+1,j)-eta(i,j))/dx
            elseif (i==mx) then
                hetax(i,j)= q(1,i,j)*(eta(i,j)-eta(i-1,j))/dx
                detax(i,j)=-slope(i,j)*(eta(i,j)-eta(i-1,j))/dx
            else
                hetax(i,j)= q(1,i,j)*(eta(i+1,j)-eta(i-1,j))/dx/2.
                detax(i,j)=-slope(i,j)*(eta(i+1,j)-eta(i-1,j))/dx/2.
            endif
         endif
         if (slope(i,j+1)<0.d0.and.slope(i,j-1)<0.d0.and.slope(i,j)<0.d0) then
            if (j==1) then
                hetay(i,j)= q(1,i,j)*(eta(i,j+1)-eta(i,j))/dy
                detay(i,j)=-slope(i,j)*(eta(i,j+1)-eta(i,j))/dy
            elseif (j==my) then
                hetay(i,j)= q(1,i,j)*(eta(i,j)-eta(i,j-1))/dy
                detay(i,j)=-slope(i,j)*(eta(i,j)-eta(i,j-1))/dy
            else
                hetay(i,j)= q(1,i,j)*(eta(i,j+1)-eta(i,j-1))/dy/2.
                detay(i,j)=-slope(i,j)*(eta(i,j+1)-eta(i,j-1))/dy/2.
            endif   
         endif
        enddo
     enddo
     
     s1=0.d0
     s2=0.d0
     s1_H=0.d0
     s2_H=0.d0

     do i=1,mx
        do j=1,my
          if (i==1) then
             s1(i,j)= (hu2(i+1,j)-hu2(i,j))/dx
             s2(i,j)= g*hetay(i,j)+(huv(i+1,j)-huv(i,j))/dx
          elseif (i==mx) then
             s1(i,j)= (hu2(i,j)-hu2(i-1,j))/dx
             s2(i,j)= g*hetay(i,j)+(huv(i,j)-huv(i-1,j))/dx
          else
             s1(i,j)= (hu2(i+1,j)-hu2(i-1,j))/2.d0/dx
             s2(i,j)= g*hetay(i,j)+(huv(i+1,j)-huv(i-1,j))/2.d0/dx
          endif
          
          if (j==1) then
             s1(i,j)= s1(i,j)+g*hetax(i,j)+(huv(i,j+1)-huv(i,j))/dy
             s2(i,j)= s2(i,j)+(hv2(i,j+1)-hv2(i,j))/dy
          elseif (j==my) then
             s1(i,j)= s1(i,j)+g*hetax(i,j)+(huv(i,j)-huv(i,j-1))/dy
             s2(i,j)= s2(i,j)+(hv2(i,j)-hv2(i,j-1))/dy
          else
             s1(i,j)= s1(i,j)+g*hetax(i,j)+(huv(i,j+1)-huv(i,j-1))/2.d0/dy
             s2(i,j)= s2(i,j)+(hv2(i,j+1)-hv2(i,j-1))/2.d0/dy
          endif
      
            if (slope(i,j)<-1d-1) then
              s1_H(i,j)=s1(i,j)/(-slope(i,j))
              s2_H(i,j)=s2(i,j)/(-slope(i,j))
           else
              s1_H(i,j)=0.d0
              s2_H(i,j)=0.d0
           endif
        enddo
     enddo
     
      tol = 1d-8

      psi(:)=0.d0

     do j=1,my
        do i=1,mx
           k= i+j*(mx+2)+1

           topo =-slope(i,j)

           if (q(1,i,j).lt.1.d-1) then
              psi(k) = 0.d0
              psi(k+mxy) = 0.d0

           elseif (slope(i,j)> sw_depth) then
              psi(k) = 0.d0
              psi(k+mxy) = 0.d0
         
           elseif (use_bous_sw_thresh.and.abs(q(1,i,j)-topo)>bous_sw_thresh*topo) then
              psi(k) = 0.d0
              psi(k+mxy) = 0.d0

           elseif (abs(hbar_xtt(i,j))>1.d-14.or.abs(hbar_ytt(i,j))>1.d-14) then
              psi(k) = 0.d0
              psi(k+mxy) = 0.d0

           elseif (maxval(abs(hbar_xtt(i-2:i+2,j-2:j+2)))>1.d-14 .or. &
                   maxval(abs(hbar_ytt(i-2:i+2,j-2:j+2)))>1.d-14    ) then
              psi(k) = 0.d0
              psi(k+mxy) = 0.d0

           elseif (maxval(slope(i-2:i+2,j-2:j+2))> sw_depth ) then
            
              psi(k)=0.d0
              psi(k+mxy)=0.d0

           elseif (maxval(q(1,i-2:i+2,j-2:j+2)) < 1.d-1) then
              psi(k) = 0.d0
              psi(k+mxy) = 0.d0
            
           else
            
            if (i==1) then
               
               if (j==1) then
                  
                     psi(k)=(B_param+.5d0)*hh(i,j)*( &
                      (s1(i+2,j)-2.d0*s1(i+1,j)+s1(i,j))/dx**2. &
                     +(s2(i+1,j+1)-s2(i+1,j)-s2(i,j+1)+s2(i,j)) &
                      /dx/dy) &
                     -B_param*g*HH(i,j)*( &
                      (detax(i+2,j)-2.d0*detax(i+1,j)+detax(i,j))/dx**2 &
                     +(detay(i+1,j+1)-detay(i+1,j) &
                      -detay(i,j+1)+detay(i,j))/dx/dy) &
                     -HHH(i,j)/6.d0*(s1_H(i+2,j)-2.d0*s1_H(i+1,j) &
                     +s1_H(i,j))/dx**2 &
                     -HHH(i,j)/6.d0*(s2_H(i+1,j+1)-s2_H(i+1,j) &
                     -s2_H(i,j+1)+s2_H(i,j))/dx/dy &
                     -1.d0/2.d0*hh(i,j)*hbar_xtt(i,j)
         
                    k= k + mxy
           
                    psi(k)=(B_param+.5d0)*hh(i,j)*( &
                      (s2(i,j+2)-2.d0*s2(i,j+1)+s2(i,j))/dy**2 &
                     +(s1(i+1,j+1)-s1(i+1,j)-s1(i,j+1)+s1(i,j) &
                      )/dx/dy) &
                     -B_param*g*HH(i,j)*( &
                      (detay(i,j+2)-2.d0*detay(i,j+1)+detay(i,j))/dy**2 &
                     +(detax(i+1,j+1)-detax(i+1,j) &
                      -detax(i,j+1)+detax(i,j))/dx/dy) &
                     -HHH(i,j)/6.d0*(s2_H(i,j+2)-2.d0*s2_H(i,j+1) &
                     +s2_H(i,j))/dy**2 &
                     -HHH(i,j)/6.d0*(s1_H(i+1,j+1)-s1_H(i+1,j) &
                     -s1_H(i,j+1)+s1_H(i,j))/dx/dy &
                     -1.d0/2.d0*hh(i,j)*hbar_ytt(i,j)
                  
               elseif (j==my) then
                  
                     psi(k)=(B_param+.5d0)*hh(i,j)*( &
                      (s1(i+2,j)-2.d0*s1(i+1,j)+s1(i,j))/dx**2. &
                     +(s2(i+1,j)-s2(i+1,j-1)-s2(i,j)+s2(i,j-1)) &
                      /dx/dy) &
                     -B_param*g*HH(i,j)*( &
                      (detax(i+2,j)-2.d0*detax(i+1,j)+detax(i,j))/dx**2 &
                     +(detay(i+1,j)-detay(i+1,j-1) &
                      -detay(i,j)+detay(i,j-1))/dx/dy) &
                     -HHH(i,j)/6.d0*(s1_H(i+2,j)-2.d0*s1_H(i+1,j) &
                     +s1_H(i,j))/dx**2 &
                     -HHH(i,j)/6.d0*(s2_H(i+1,j)-s2_H(i+1,j-1) &
                     -s2_H(i,j)+s2_H(i,j-1))/dx/dy &
                     -1.d0/2.d0*hh(i,j)*hbar_xtt(i,j)
         
                    k= k + mxy
           
                    psi(k)=(B_param+.5d0)*hh(i,j)*( &
                      (s2(i,j)-2.d0*s2(i,j-1)+s2(i,j-2))/dy**2 &
                     +(s1(i+1,j)-s1(i+1,j-1)-s1(i,j)+s1(i,j-1) &
                      )/dx/dy) &
                     -B_param*g*HH(i,j)*( &
                      (detay(i,j)-2.d0*detay(i,j-1)+detay(i,j-2))/dy**2 &
                     +(detax(i+1,j)-detax(i+1,j-1) &
                      -detax(i,j)+detax(i,j-1))/dx/dy) &
                     -HHH(i,j)/6.d0*(s2_H(i,j)-2.d0*s2_H(i,j-1) &
                     +s2_H(i,j-2))/dy**2 &
                     -HHH(i,j)/6.d0*(s1_H(i+1,j)-s1_H(i+1,j-1) &
                     -s1_H(i,j)+s1_H(i,j-1))/dx/dy &
                     -1.d0/2.d0*hh(i,j)*hbar_ytt(i,j)
                  
               else
                  
                     psi(k)=(B_param+.5d0)*hh(i,j)*( &
                      (s1(i+2,j)-2.d0*s1(i+1,j)+s1(i,j))/dx**2. &
                     +(s2(i+1,j+1)-s2(i+1,j-1)-s2(i,j+1)+s2(i,j-1)) &
                      /2.d0/dx/dy) &
                     -B_param*g*HH(i,j)*( &
                      (detax(i+2,j)-2.d0*detax(i+1,j)+detax(i,j))/dx**2 &
                     +(detay(i+1,j+1)-detay(i+1,j-1) &
                      -detay(i,j+1)+detay(i,j-1))/2.d0/dx/dy) &
                     -HHH(i,j)/6.d0*(s1_H(i+2,j)-2.d0*s1_H(i+1,j) &
                     +s1_H(i,j))/dx**2 &
                     -HHH(i,j)/6.d0*(s2_H(i+1,j+1)-s2_H(i+1,j-1) &
                     -s2_H(i,j+1)+s2_H(i,j-1))/2.d0/dx/dy &
                     -1.d0/2.d0*hh(i,j)*hbar_xtt(i,j)
         
                    k= k + mxy
           
                    psi(k)=(B_param+.5d0)*hh(i,j)*( &
                      (s2(i,j+1)-2.d0*s2(i,j)+s2(i,j-1))/dy**2 &
                     +(s1(i+1,j+1)-s1(i+1,j-1)-s1(i,j+1)+s1(i,j-1) &
                      )/2.d0/dx/dy) &
                     -B_param*g*HH(i,j)*( &
                      (detay(i,j+1)-2.d0*detay(i,j)+detay(i,j-1))/dy**2 &
                     +(detax(i+1,j+1)-detax(i+1,j-1) &
                      -detax(i,j+1)+detax(i,j-1))/2.d0/dx/dy) &
                     -HHH(i,j)/6.d0*(s2_H(i,j+1)-2.d0*s2_H(i,j) &
                     +s2_H(i,j-1))/dy**2 &
                     -HHH(i,j)/6.d0*(s1_H(i+1,j+1)-s1_H(i+1,j-1) &
                     -s1_H(i,j+1)+s1_H(i,j-1))/2.d0/dx/dy &
                     -1.d0/2.d0*hh(i,j)*hbar_ytt(i,j)
                  
               endif
               
            elseif (i==mx) then
               
               if (j==1) then
                  
                     psi(k)=(B_param+.5d0)*hh(i,j)*( &
                      (s1(i,j)-2.d0*s1(i-1,j)+s1(i-2,j))/dx**2. &
                     +(s2(i,j+1)-s2(i,j)-s2(i-1,j+1)+s2(i-1,j)) &
                      /dx/dy) &
                     -B_param*g*HH(i,j)*( &
                      (detax(i,j)-2.d0*detax(i-1,j)+detax(i-2,j))/dx**2 &
                     +(detay(i,j+1)-detay(i,j) &
                      -detay(i-1,j+1)+detay(i-1,j))/dx/dy) &
                     -HHH(i,j)/6.d0*(s1_H(i,j)-2.d0*s1_H(i-1,j) &
                     +s1_H(i-2,j))/dx**2 &
                     -HHH(i,j)/6.d0*(s2_H(i,j+1)-s2_H(i,j) &
                     -s2_H(i-1,j+1)+s2_H(i-1,j))/dx/dy &
                     -1.d0/2.d0*hh(i,j)*hbar_xtt(i,j)
         
                    k= k + mxy
           
                    psi(k)=(B_param+.5d0)*hh(i,j)*( &
                      (s2(i,j+2)-2.d0*s2(i,j+1)+s2(i,j))/dy**2 &
                     +(s1(i,j+1)-s1(i,j)-s1(i-1,j+1)+s1(i-1,j) &
                      )/dx/dy) &
                     -B_param*g*HH(i,j)*( &
                      (detay(i,j+2)-2.d0*detay(i,j+1)+detay(i,j))/dy**2 &
                     +(detax(i,j+1)-detax(i,j) &
                      -detax(i-1,j+1)+detax(i-1,j))/dx/dy) &
                     -HHH(i,j)/6.d0*(s2_H(i,j+2)-2.d0*s2_H(i,j+1) &
                     +s2_H(i,j))/dy**2 &
                     -HHH(i,j)/6.d0*(s1_H(i,j+1)-s1_H(i,j) &
                     -s1_H(i-1,j+1)+s1_H(i-1,j))/dx/dy &
                     -1.d0/2.d0*hh(i,j)*hbar_ytt(i,j)
                  
               elseif (j==my) then
                  
                     psi(k)=(B_param+.5d0)*hh(i,j)*( &
                      (s1(i,j)-2.d0*s1(i-1,j)+s1(i-2,j))/dx**2. &
                     +(s2(i,j)-s2(i,j-1)-s2(i-1,j)+s2(i-1,j-1)) &
                      /dx/dy) &
                     -B_param*g*HH(i,j)*( &
                      (detax(i,j)-2.d0*detax(i-1,j)+detax(i-2,j))/dx**2 &
                     +(detay(i,j)-detay(i,j-1) &
                      -detay(i-1,j)+detay(i-1,j-1))/dx/dy) &
                     -HHH(i,j)/6.d0*(s1_H(i,j)-2.d0*s1_H(i-1,j) &
                     +s1_H(i-2,j))/dx**2 &
                     -HHH(i,j)/6.d0*(s2_H(i,j)-s2_H(i,j-1) &
                     -s2_H(i-1,j)+s2_H(i-1,j-1))/dx/dy &
                     -1.d0/2.d0*hh(i,j)*hbar_xtt(i,j)
         
                    k= k + mxy
           
                    psi(k)=(B_param+.5d0)*hh(i,j)*( &
                      (s2(i,j)-2.d0*s2(i,j-1)+s2(i,j-2))/dy**2 &
                     +(s1(i,j)-s1(i,j-1)-s1(i-1,j)+s1(i-1,j-1) &
                      )/dx/dy) &
                     -B_param*g*HH(i,j)*( &
                      (detay(i,j)-2.d0*detay(i,j-1)+detay(i,j-2))/dy**2 &
                     +(detax(i,j)-detax(i,j-1) &
                      -detax(i-1,j)+detax(i-1,j-1))/dx/dy) &
                     -HHH(i,j)/6.d0*(s2_H(i,j)-2.d0*s2_H(i,j-1) &
                     +s2_H(i,j-2))/dy**2 &
                     -HHH(i,j)/6.d0*(s1_H(i,j)-s1_H(i,j-1) &
                     -s1_H(i-1,j)+s1_H(i-1,j-1))/dx/dy &
                     -1.d0/2.d0*hh(i,j)*hbar_ytt(i,j)
                  
               else
                     psi(k)=(B_param+.5d0)*hh(i,j)*( &
                      (s1(i,j)-2.d0*s1(i-1,j)+s1(i-2,j))/dx**2. &
                     +(s2(i,j+1)-s2(i,j-1)-s2(i-1,j+1)+s2(i-1,j-1)) &
                      /2.d0/dx/dy) &
                     -B_param*g*HH(i,j)*( &
                      (detax(i,j)-2.d0*detax(i-1,j)+detax(i-2,j))/dx**2 &
                     +(detay(i,j+1)-detay(i,j-1) &
                      -detay(i-1,j+1)+detay(i-1,j-1))/2.d0/dx/dy) &
                     -HHH(i,j)/6.d0*(s1_H(i,j)-2.d0*s1_H(i-1,j) &
                     +s1_H(i-2,j))/dx**2 &
                     -HHH(i,j)/6.d0*(s2_H(i,j+1)-s2_H(i,j-1) &
                     -s2_H(i-1,j+1)+s2_H(i-1,j-1))/2.d0/dx/dy &
                     -1.d0/2.d0*hh(i,j)*hbar_xtt(i,j)
         
                    k= k + mxy
           
                    psi(k)=(B_param+.5d0)*hh(i,j)*( &
                      (s2(i,j+1)-2.d0*s2(i,j)+s2(i,j-1))/dy**2 &
                     +(s1(i,j+1)-s1(i,j-1)-s1(i-1,j+1)+s1(i-1,j-1) &
                      )/2.d0/dx/dy) &
                     -B_param*g*HH(i,j)*( &
                      (detay(i,j+1)-2.d0*detay(i,j)+detay(i,j-1))/dy**2 &
                     +(detax(i,j+1)-detax(i,j-1) &
                      -detax(i-1,j+1)+detax(i-1,j-1))/2.d0/dx/dy) &
                     -HHH(i,j)/6.d0*(s2_H(i,j+1)-2.d0*s2_H(i,j) &
                     +s2_H(i,j-1))/dy**2 &
                     -HHH(i,j)/6.d0*(s1_H(i,j+1)-s1_H(i,j-1) &
                     -s1_H(i-1,j+1)+s1_H(i-1,j-1))/2.d0/dx/dy &
                     -1.d0/2.d0*hh(i,j)*hbar_ytt(i,j)
                  
               endif
               
            elseif (j==1) then
               
                  psi(k)=(B_param+.5d0)*hh(i,j)*( &
                   (s1(i+1,j)-2.d0*s1(i,j)+s1(i-1,j))/dx**2. &
                  +(s2(i+1,j+1)-s2(i+1,j)-s2(i-1,j+1)+s2(i-1,j)) &
                   /2.d0/dx/dy) &
                  -B_param*g*HH(i,j)*( &
                   (detax(i+1,j)-2.d0*detax(i,j)+detax(i-1,j))/dx**2 &
                  +(detay(i+1,j+1)-detay(i+1,j) &
                   -detay(i-1,j+1)+detay(i-1,j))/2.d0/dx/dy) &
                  -HHH(i,j)/6.d0*(s1_H(i+1,j)-2.d0*s1_H(i,j) &
                  +s1_H(i-1,j))/dx**2 &
                  -HHH(i,j)/6.d0*(s2_H(i+1,j+1)-s2_H(i+1,j) &
                  -s2_H(i-1,j+1)+s2_H(i-1,j))/2.d0/dx/dy &
                  -1.d0/2.d0*hh(i,j)*hbar_xtt(i,j)
         
                 k= k + mxy
           
                 psi(k)=(B_param+.5d0)*hh(i,j)*( &
                   (s2(i,j+2)-2.d0*s2(i,j+1)+s2(i,j))/dy**2 &
                  +(s1(i+1,j+1)-s1(i+1,j)-s1(i-1,j+1)+s1(i-1,j) &
                   )/2.d0/dx/dy) &
                  -B_param*g*HH(i,j)*( &
                   (detay(i,j+2)-2.d0*detay(i,j+1)+detay(i,j))/dy**2 &
                  +(detax(i+1,j+1)-detax(i+1,j) &
                   -detax(i-1,j+1)+detax(i-1,j))/2.d0/dx/dy) &
                  -HHH(i,j)/6.d0*(s2_H(i,j+2)-2.d0*s2_H(i,j+1) &
                  +s2_H(i,j))/dy**2 &
                  -HHH(i,j)/6.d0*(s1_H(i+1,j+1)-s1_H(i+1,j) &
                  -s1_H(i-1,j+1)+s1_H(i-1,j))/2.d0/dx/dy &
                  -1.d0/2.d0*hh(i,j)*hbar_ytt(i,j)
               
            elseif (j==my) then
               
                  psi(k)=(B_param+.5d0)*hh(i,j)*( &
                   (s1(i+1,j)-2.d0*s1(i,j)+s1(i-1,j))/dx**2. &
                  +(s2(i+1,j)-s2(i+1,j-1)-s2(i-1,j)+s2(i-1,j-1)) &
                   /2.d0/dx/dy) &
                  -B_param*g*HH(i,j)*( &
                   (detax(i+1,j)-2.d0*detax(i,j)+detax(i-1,j))/dx**2 &
                  +(detay(i+1,j)-detay(i+1,j-1) &
                   -detay(i-1,j)+detay(i-1,j-1))/2.d0/dx/dy) &
                  -HHH(i,j)/6.d0*(s1_H(i+1,j)-2.d0*s1_H(i,j) &
                  +s1_H(i-1,j))/dx**2 &
                  -HHH(i,j)/6.d0*(s2_H(i+1,j)-s2_H(i+1,j-1) &
                  -s2_H(i-1,j)+s2_H(i-1,j-1))/2.d0/dx/dy &
                  -1.d0/2.d0*hh(i,j)*hbar_xtt(i,j)
         
                 k= k + mxy
           
                 psi(k)=(B_param+.5d0)*hh(i,j)*( &
                   (s2(i,j)-2.d0*s2(i,j-1)+s2(i,j-2))/dy**2 &
                  +(s1(i+1,j)-s1(i+1,j-1)-s1(i-1,j)+s1(i-1,j-1) &
                   )/2.d0/dx/dy) &
                  -B_param*g*HH(i,j)*( &
                   (detay(i,j)-2.d0*detay(i,j-1)+detay(i,j-2))/dy**2 &
                  +(detax(i+1,j)-detax(i+1,j-1) &
                   -detax(i-1,j)+detax(i-1,j-1))/2.d0/dx/dy) &
                  -HHH(i,j)/6.d0*(s2_H(i,j)-2.d0*s2_H(i,j-1) &
                  +s2_H(i,j-2))/dy**2 &
                  -HHH(i,j)/6.d0*(s1_H(i+1,j)-s1_H(i+1,j-1) &
                  -s1_H(i-1,j)+s1_H(i-1,j-1))/2.d0/dx/dy &
                  -1.d0/2.d0*hh(i,j)*hbar_ytt(i,j)
               
            else
               
                  psi(k)=(B_param+.5d0)*hh(i,j)*( &
                   (s1(i+1,j)-2.d0*s1(i,j)+s1(i-1,j))/dx**2. &
                  +(s2(i+1,j+1)-s2(i+1,j-1)-s2(i-1,j+1)+s2(i-1,j-1)) &
                   /4.d0/dx/dy) &
                  -B_param*g*HH(i,j)*( &
                   (detax(i+1,j)-2.d0*detax(i,j)+detax(i-1,j))/dx**2 &
                  +(detay(i+1,j+1)-detay(i+1,j-1) &
                   -detay(i-1,j+1)+detay(i-1,j-1))/4.d0/dx/dy) &
                  -HHH(i,j)/6.d0*(s1_H(i+1,j)-2.d0*s1_H(i,j) &
                  +s1_H(i-1,j))/dx**2 &
                  -HHH(i,j)/6.d0*(s2_H(i+1,j+1)-s2_H(i+1,j-1) &
                  -s2_H(i-1,j+1)+s2_H(i-1,j-1))/4.d0/dx/dy &
                  -1.d0/2.d0*hh(i,j)*hbar_xtt(i,j)
  
                 k= k + mxy
           
                 psi(k)=(B_param+.5d0)*hh(i,j)*( &
                   (s2(i,j+1)-2.d0*s2(i,j)+s2(i,j-1))/dy**2 &
                  +(s1(i+1,j+1)-s1(i+1,j-1)-s1(i-1,j+1)+s1(i-1,j-1) &
                   )/4.d0/dx/dy) &
                  -B_param*g*HH(i,j)*( &
                   (detay(i,j+1)-2.d0*detay(i,j)+detay(i,j-1))/dy**2 &
                  +(detax(i+1,j+1)-detax(i+1,j-1) &
                   -detax(i-1,j+1)+detax(i-1,j-1))/4.d0/dx/dy) &
                  -HHH(i,j)/6.d0*(s2_H(i,j+1)-2.d0*s2_H(i,j) &
                  +s2_H(i,j-1))/dy**2 &
                  -HHH(i,j)/6.d0*(s1_H(i+1,j+1)-s1_H(i+1,j-1) &
                  -s1_H(i-1,j+1)+s1_H(i-1,j-1))/4.d0/dx/dy &
                  -1.d0/2.d0*hh(i,j)*hbar_ytt(i,j)
               
            endif
         
            endif
         
 993  continue
            
        enddo
     enddo

     if (.not.use_bous_sw_thresh) go to 994

      do j=2,my-1
         do i=2,mx-1
            k= i+j*(mx+2)+1

            topo = -slope(i,j)

            if (q(1,i,j)>1.d-1.and.topo>1.d-1.and. &
                q(1,i,j)-topo>bous_sw_thresh*topo) then
 
               psi(k) = 0.d0
               psi(k+mxy) = 0.d0
               !psi(:)=0.d0
               !go to 994

               xaxis1: do kk=i,1,-1
                  if (eta(kk,j)-eta(kk-1,j)> 0.d0 ) then
                  !if ((q(1,kk,j)+aux(1,kk,j))> sea_level ) then
                     k= kk+j*(mx+2)+1
                     psi(k) = 0.d0
                     psi(k+mxy) = 0.d0
                  else
                     exit xaxis1
                  endif
               end do xaxis1

               xaxis2: do kk=i,mx
                  if (eta(kk+1,j)-eta(kk,j)< 0.d0 ) then
                  !if ((q(1,kk,j)+aux(1,kk,j))> sea_level ) then
                     k= kk+j*(mx+2)+1
                     psi(k) = 0.d0
                     psi(k+mxy) = 0.d0
                  else
                     exit xaxis2
                  endif
               end do xaxis2

               yaxis1: do kk=j,1,-1
                  if (eta(i,kk)-eta(i,kk-1)> 0.d0 ) then
                  !if ((q(1,i,kk)+aux(1,i,kk))> sea_level ) then
                     k= i+kk*(mx+2)+1
                     psi(k) = 0.d0
                     psi(k+mxy) = 0.d0
                  else
                     exit yaxis1
                  endif
               end do yaxis1

               yaxis2: do kk=j,my
                  if (eta(i,kk+1)-eta(i,kk)< 0.d0 ) then
                  !if ((q(1,i,kk)+aux(1,i,kk))> sea_level ) then
                     k= i+kk*(mx+2)+1
                     psi(k) = 0.d0
                     psi(k+mxy) = 0.d0
                  else
                     exit yaxis2
                  endif
               end do yaxis2

            elseif (q(1,i,j)>1.d-1.and.topo>1.d-1.and. &
                q(1,i,j)-topo<-bous_sw_thresh*topo) then
 
               psi(k) = 0.d0
               psi(k+mxy) = 0.d0

               xaxis11: do kk=i,1,-1
                  if (eta(kk,j)-eta(kk-1,j)< 0.d0 ) then
                     k= kk+j*(mx+2)+1
                     psi(k) = 0.d0
                     psi(k+mxy) = 0.d0
                  else
                     exit xaxis11
                  endif
               end do xaxis11

               xaxis12: do kk=i,mx
                  if (eta(kk+1,j)-eta(kk,j)> 0.d0 ) then
                     k= kk+j*(mx+2)+1
                     psi(k) = 0.d0
                     psi(k+mxy) = 0.d0
                  else
                     exit xaxis12
                  endif
               end do xaxis12

               yaxis11: do kk=j,1,-1
                  if (eta(i,kk)-eta(i,kk-1)< 0.d0 ) then
                     k= i+kk*(mx+2)+1
                     psi(k) = 0.d0
                     psi(k+mxy) = 0.d0
                  else
                     exit yaxis11
                  endif
               end do yaxis11

               yaxis12: do kk=j,my
                  if (eta(i,kk+1)-eta(i,kk)> 0.d0 ) then
                     k= i+kk*(mx+2)+1
                     psi(k) = 0.d0
                     psi(k+mxy) = 0.d0
                  else
                     exit yaxis12
                  endif
               end do yaxis12


            endif

 995     continue

         end do
      end do

 994  continue

      return
      end
