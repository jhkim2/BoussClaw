
  subroutine src2_bous(meqn,mbc,mx,my, &
                  xlower,ylower,dx,dy,q,maux,aux,t,dt)

      use geoclaw_module, only: g => grav, coriolis_forcing, coriolis
      use geoclaw_module, only: friction_forcing, friction_depth
      use geoclaw_module, only: manning_coefficient,sea_level
      use bous_module

      implicit none
            
      integer(kind=4) i,j,mx,my,k,ii,nstep
      integer(kind=4), parameter :: ndim = 1
      
      real(kind=8)   q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      real(kind=8), allocatable, dimension(:,:) ::  q0(:,:)
      real(kind=8), allocatable, dimension(:,:) :: q1d(:,:)
      real(kind=8), allocatable, dimension(:,:) :: aux1d(:,:)
      real(kind=8), allocatable, dimension(:,:) :: xx0(:,:)
      real(kind=8), allocatable, dimension(:) :: psi(:)
            
      real(kind=8) :: delt

      real(kind=8) tol,t,dt,x,y
      real(kind=8) dx,dy,xlower,ylower
      integer(kind=4) INFO,maux,mbc,meqn
     
      INTEGER            LDB
      INTEGER, allocatable, dimension(:) :: IPIV( : )
      real(kind=8), allocatable, dimension(:) :: D(:), DL(:), DU(:), DU2(:)             
      
      logical verbose

      verbose=.False.
           
!      ----------------------------------------------------------------
!
!      Boussinesq type-------------------------------------------------
    
     nstep=ceiling(max(dt/dx,dt/dy)*2*sqrt(2.))

     delt=dt/nstep
  
     do ii=1,nstep

        !----------------------------------------------------------
        ! y-sweep

        LDB = my+2

        allocate( q1d(meqn,1-mbc:my+mbc) )
        allocate( q0(meqn,1-mbc:my+mbc) )
        allocate( aux1d(maux,1-mbc:my+mbc) )
        allocate( xx0(1:my,4) )
        allocate( psi( LDB ) )
        allocate( IPIV( LDB ))
        allocate( D(LDB) )
        allocate( DL(LDB-1) )
        allocate( DU(LDB-1) )
        allocate( DU2(LDB-2) )

        do i=1,mx

          q1d = q(:,i,:)
          aux1d= aux(:,i,:)

          call read_diag(my,meqn,mbc,dy,q1d,maux,aux1d, DL, D, DU)

          call DGTTRF( LDB , DL, D, DU, DU2, IPIV, INFO )
    
          psi = 0.d0
          xx0 = 0.d0
          q0  = 0.d0
          q0(1,:)  = q(1,i,:)
          q0(2,:)  = q(3,i,:)

          ! RK4   

          ! First Stage
      
          call read_psi(my,meqn,mbc,dy,q0,maux,aux1d,psi,g)

          call DGTTRS( 'N', LDB , 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )

          XX0(1:my,1) = psi(2:my+1)
        
          q0(2,1:my)=q(3,i,1:my)-delt/2.d0*xx0(1:my,1)
        
          ! Second Stage

          call read_psi(my,meqn,mbc,dy,q0,maux,aux1d,psi,g)

          call DGTTRS( 'N', my+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )

          XX0(1:my,2)=psi(2:my+1)
        
          q0(2,1:my)=q(3,i,1:my)-delt/2.d0*xx0(1:my,2)
        
          ! Third Stage
        
          call read_psi(my,meqn,mbc,dy,q0,maux,aux1d,psi,g)

          call DGTTRS( 'N', my+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )
                
          XX0(1:my,3)=psi(2:my+1)
        
          q0(2,1:my)=q(3,i,1:my)-delt*xx0(1:my,3)
        
          ! Fourth Stage
        
          call read_psi(my,meqn,mbc,dy,q0,maux,aux1d,psi,g)

          call DGTTRS( 'N', my+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )

          XX0(1:my,4)=psi(2:my+1)

!=======================================================================
                              
            q(3,i,1:my) = q(3,i,1:my)- delt/6.d0*( xx0(1:my,1) + 2.d0*xx0(1:my,2) &
               +2.d0*xx0(1:my,3) + xx0(1:my,4))

        enddo

        deallocate( q1d )
        deallocate( q0 )
        deallocate( aux1d )
        deallocate( xx0,psi )
        deallocate( IPIV )
        deallocate( D, DL, DU, DU2 )

        !--------------------------------------------------------------------
        ! x-sweep

        LDB = mx+2

        allocate( q1d(meqn,1-mbc:mx+mbc) )
        allocate( q0(meqn,1-mbc:mx+mbc) )
        allocate( aux1d(maux,1-mbc:mx+mbc) )
        allocate( xx0(1:mx,4) )
        allocate( psi(mx+2) )
        allocate( IPIV( LDB ))
        allocate( D(LDB) )
        allocate( DL(LDB-1) )
        allocate( DU(LDB-1) )
        allocate( DU2(LDB-2) )

        do j=1,my

          q1d = q(:,:,j)
          aux1d= aux(:,:,j)

          call read_diag(mx,meqn,mbc,dx,q1d,maux,aux1d, DL, D, DU)

          call DGTTRF( mx+2, DL, D, DU, DU2, IPIV, INFO )
    
          psi = 0.d0
          xx0 = 0.d0
          q0  = q(:,:,j)

          ! RK4   

          ! First Stage
      
          call read_psi(mx,meqn,mbc,dx,q0,maux,aux1d,psi,g)

          call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )

          XX0(1:mx,1) = psi(2:mx+1)
        
          q0(2,1:mx)=q(2,1:mx,j)-delt/2.d0*xx0(1:mx,1)
        
          ! Second Stage

          call read_psi(mx,meqn,mbc,dx,q0,maux,aux1d,psi,g)

          call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )

          XX0(1:mx,2)=psi(2:mx+1)
        
          q0(2,1:mx)=q(2,1:mx,j)-delt/2.d0*xx0(1:mx,2)
        
          ! Third Stage
        
          call read_psi(mx,meqn,mbc,dx,q0,maux,aux1d,psi,g)

          call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )
                
          XX0(1:mx,3)=psi(2:mx+1)
        
          q0(2,1:mx)=q(2,1:mx,j)-delt*xx0(1:mx,3)
        
          ! Fourth Stage
        
          call read_psi(mx,meqn,mbc,dx,q0,maux,aux1d,psi,g)

          call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )

          XX0(1:mx,4)=psi(2:mx+1)

!=======================================================================
                              
            q(2,1:mx,j) = q(2,1:mx,j)- delt/6.d0*( xx0(1:mx,1) + 2.d0*xx0(1:mx,2) &
               +2.d0*xx0(1:mx,3) + xx0(1:mx,4))

        enddo

        deallocate( q1d )
        deallocate( q0 )
        deallocate( aux1d )
        deallocate( xx0,psi )
        deallocate( IPIV )
        deallocate( D, DL, DU, DU2 )
     
      enddo
      
 999  continue
      
  end subroutine src2_bous
  
! =========================================================
      subroutine read_diag(mx,meqn,mbc,dx,q,maux,aux,DL,D,DU)
! =========================================================
      use bous_module
      use geoclaw_module, only: sea_level

      implicit none
     
      integer(kind=4) mx,meqn,mbc,maux
      integer(kind=4) i,k
      real(kind=8) dx

      real(kind=8)   q(meqn,1-mbc:mx+mbc)
      real(kind=8) aux(maux,1-mbc:mx+mbc)

      DOUBLE PRECISION   D(mx+2), DL(mx+1), DU(mx+1)                
      
      real(kind=8)  HH(1-mbc:mx+mbc)
      real(kind=8) HHH(1-mbc:mx+mbc)
      real(kind=8) slope(1-mbc:mx+mbc)
      
      slope = aux(1,:)-sea_level  
     
            do i=1-mbc,mx+mbc
               HH(i)=(max(0.,-slope(i)))**2
              HHH(i)=(max(0.,-slope(i)))**3
            enddo
      
        D = 1.d0
        DU= 0.d0
        DL= 0.d0

  
        do i=1,mx
        
            if (q(1,i).lt.1.d-1.or.slope(i)> sw_depth) then

              ! do nothing
        
            elseif (maxval(slope(i-2:i+2))  > sw_depth) then

              ! do nothing

            elseif (maxval(q(1,i-2:i+2)) < 1.d-1) then

              ! do nothing

            elseif (maxval(slope(i-1:i+1)) <= sw_depth) then              

      ! D1 part     

              D(i+1) = 1.d0 + 2.d0*(B_param+.5d0)*HH(i)/dx**2 &
              -2.d0/6.d0*HHH(i)/(-slope(i))/dx**2
     
              DU(i+1)=-(B_param+.5d0)*HH(i)/dx**2 &
              +1.d0/6.d0*HHH(i)/(-slope(i+1))/dx**2
     
              DL(i)=-(B_param+.5d0)*HH(i)/dx**2 &
              +1.d0/6.d0*HHH(i)/(-slope(i-1))/dx**2
     
            endif
            
        enddo
             
      return
      end
   
!======================================================================
      subroutine read_psi(mx,meqn,mbc,dx,q1d,maux,aux1d,psi,g)

      use geoclaw_module, only: sea_level
      use bous_module
     
      implicit none

      integer(kind=4) mx,meqn,mbc,maux,i,j,k,kk
      real   (kind=8) dx,g

      real(kind=8)   q1d(meqn,1-mbc:mx+mbc)
      real(kind=8) aux1d(maux,1-mbc:mx+mbc)
      real(kind=8) slope(1-mbc:mx+mbc)
      real(kind=8) psi(mx+2)

      real(kind=8)   hh(1-mbc:mx+mbc)
      real(kind=8)  hhh(1-mbc:mx+mbc)
      real(kind=8)  hetax(-1:mx+2)
      real(kind=8)  eta(1-mbc:mx+mbc)
      real(kind=8)  hu2(1-mbc:mx+mbc)
      real(kind=8)  detax(-1:mx+2)
      real(kind=8)  s1(0:mx+1)
      real(kind=8)  s1_H(0:mx+1)
      real(kind=8)  tol,topo
  
      slope = aux1d(1,:)-sea_level
   
      do i=1-mbc,mx+mbc
           if (q1d(1,i).gt.1d-4) then
              hu2(i)= q1d(2,i)**2/q1d(1,i)
           else
              hu2(i)= 0.d0
           endif
           eta(i)= q1d(1,i) + slope(i) 
            HH(i)=(max(0.d0,-slope(i)))**2
           HHH(i)=(max(0.d0,-slope(i)))**3
      enddo
     
      hetax = 0.d0
      detax = 0.d0

      do i= 1,mx
         if (slope(i+1)<0.d0.and.slope(i-1)<0.d0.and.slope(i)<0.d0) then
            if (i==1) then
                hetax(i)= q1d(1,i)*(eta(i+1)-eta(i))/dx
                detax(i)=-slope(i)*(eta(i+1)-eta(i))/dx
            elseif (i==mx) then
                hetax(i)= q1d(1,i)*(eta(i)-eta(i-1))/dx
                detax(i)=-slope(i)*(eta(i)-eta(i-1))/dx
            else
                hetax(i)= q1d(1,i)*(eta(i+1)-eta(i-1))/dx/2.
                detax(i)=-slope(i)*(eta(i+1)-eta(i-1))/dx/2.
            endif
         endif
      enddo

      s1=0.d0
      s1_H=0.d0

      do i=1,mx
          if (i==1) then
             s1(i)= (hu2(i+1)-hu2(i))/dx
          elseif (i==mx) then
             s1(i)= (hu2(i)-hu2(i-1))/dx
          else
             s1(i)= (hu2(i+1)-hu2(i-1))/2.d0/dx
          endif

          s1(i)= s1(i)+g*hetax(i)
      
          if (slope(i)<-1d-1) then
             s1_H(i)=s1(i)/(-slope(i))
          else
             s1_H(i)=0.d0
          endif
      enddo
     
      tol = 1d-8

      psi=0.d0

        do i=1,mx

           k = i+1

           topo =-slope(i)

           if (q1d(1,i).lt.1.d-1) then
              psi(k) = 0.d0

           elseif (slope(i)> sw_depth) then
              psi(k) = 0.d0
         
           elseif (use_bous_sw_thresh.and.abs(q1d(1,i)-topo)>bous_sw_thresh*topo) then
              psi(k) = 0.d0

           elseif (maxval(slope(i-2:i+2))> sw_depth ) then
            
              psi(k)=0.d0

           elseif (maxval(q1d(1,i-2:i+2)) < 1.d-1) then

              psi(k) = 0.d0

           else
            
            if (i==1) then
                 
                     psi(k)=(B_param+.5d0)*hh(i)*( &
                      (s1(i+2)-2.d0*s1(i+1)+s1(i))/dx**2.) &
                     -B_param*g*hh(i)*( &
                      (detax(i)-2.d0*detax(i+1)+detax(i))/dx**2) &
                     -Hhh(i)/6.d0*(s1_H(i+2)-2.d0*s1_H(i+1) &
                     +s1_H(i))/dx**2
              
            elseif (i==mx) then

                     psi(k)=(B_param+.5d0)*hh(i)*( &
                      (s1(i)-2.d0*s1(i-1)+s1(i-2))/dx**2. ) &
                     -B_param*g*hh(i)*( &
                      (detax(i)-2.d0*detax(i-1)+detax(i-2))/dx**2 ) &
                     -Hhh(i)/6.d0*(s1_H(i)-2.d0*s1_H(i-1) &
                     +s1_H(i-2))/dx**2 
                          
            else
               
                  psi(k)=(B_param+.5d0)*hh(i)*( &
                   (s1(i+1)-2.d0*s1(i)+s1(i-1))/dx**2.) &
                  -B_param*g*hh(i)*( &
                   (detax(i+1)-2.d0*detax(i)+detax(i-1))/dx**2 ) &
                  -Hhh(i)/6.d0*(s1_H(i+1)-2.d0*s1_H(i) &
                  +s1_H(i-1))/dx**2 
               
            endif
           endif
                
        enddo

     if (.not.use_bous_sw_thresh) go to 994

         do i=2,mx-1
            k= i+1

            topo = -slope(i)

            if (q1d(1,i)>1.d-1.and.topo>1.d-1.and. &
                q1d(1,i)-topo>bous_sw_thresh*topo) then
 
               psi(k) = 0.d0

               xaxis1: do kk=i,1,-1
                  if (eta(kk)-eta(kk-1)> 0.d0 ) then
                  !if ((q(1,kk,j)+aux(1,kk,j))> sea_level ) then
                     psi(kk) = 0.d0
                  else
                     exit xaxis1
                  endif
               end do xaxis1

               xaxis2: do kk=i,mx
                  if (eta(kk+1)-eta(kk)< 0.d0 ) then
                  !if ((q(1,kk,j)+aux(1,kk,j))> sea_level ) then
                     psi(kk) = 0.d0
                  else
                     exit xaxis2
                  endif
               end do xaxis2

            elseif (q1d(1,i)>1.d-1.and.topo>1.d-1.and. &
                q1d(1,i)-topo<-bous_sw_thresh*topo) then
 
               psi(k) = 0.d0

               xaxis11: do kk=i,1,-1
                  if (eta(kk)-eta(kk-1)< 0.d0 ) then
                     psi(kk) = 0.d0
                  else
                     exit xaxis11
                  endif
               end do xaxis11

               xaxis12: do kk=i,mx
                  if (eta(kk+1)-eta(kk)> 0.d0 ) then
                     psi(kk) = 0.d0
                  else
                     exit xaxis12
                  endif
               end do xaxis12

            endif

         end do

 994  continue

      return
      end
