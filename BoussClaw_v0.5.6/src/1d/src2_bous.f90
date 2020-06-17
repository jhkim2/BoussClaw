
  subroutine src2_bous(meqn,mbc,mx,my, &
                  xlower,ylower,dx,dy,q,maux,aux,t,dt)

      use geoclaw_module, only: g => grav, coriolis_forcing, coriolis
      use geoclaw_module, only: friction_forcing, friction_depth
      use geoclaw_module, only: manning_coefficient,sea_level
      use bous_module

      implicit none
            
      integer(kind=4) i,j,mx,my,k,ii,nstep,kk
      integer(kind=4), parameter :: ndim = 1
      
      real(kind=8)   q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8)  q0(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) psi(mx+2)
     
      real(kind=8)   xx0(1:mx,4)
            
      real(kind=8) :: delt

      real(kind=8) tol,t,dt,x,y
      real(kind=8) dx,dy,xlower,ylower
      integer(kind=4) INFO,maux,mbc,meqn
     
      INTEGER            LDB
      INTEGER            IPIV( mx+2 )
      DOUBLE PRECISION   D(mx+2), DL(mx+1), DU(mx+1), DU2(1:mx)             
      
      logical verbose

      verbose=.False.
           
!      ----------------------------------------------------------------
!
!      Boussinesq type-------------------------------------------------
    
     nstep=ceiling(max(dt/dx,dt/dy)*2*sqrt(2.))

     k=int(my/2)

     q(3,:,:)=0.d0

     delt=dt/nstep

     LDB = mx+2
     D  =0.d0
     DU =0.d0
     DL =0.d0
     DU2=0.d0
    
     do ii=1,nstep

          call read_diag(mx,my,meqn,mbc,dx,dy,q,maux,aux, DL, D, DU)

          call DGTTRF( mx+2, DL, D, DU, DU2, IPIV, INFO )
    
          psi = 0.d0
          xx0 = 0.d0
          q0  = q

          ! RK4   

          ! First Stage
      
          call read_psi(mx,my,meqn,mbc,dx,dy,q0,maux,aux,psi,g)

          call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )

          XX0(1:mx,1) = psi(2:mx+1)
        
          do j=1-mbc,my+mbc
             q0(2,1:mx,j)=q(2,1:mx,j)-delt/2.d0*xx0(1:mx,1)
          enddo
        
          ! Second Stage

          call read_psi(mx,my,meqn,mbc,dx,dy,q0,maux,aux,psi,g )

          call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )

          XX0(1:mx,2)=psi(2:mx+1)
        
          do j=1-mbc,my+mbc
             q0(2,1:mx,j)=q(2,1:mx,j)-delt/2.d0*xx0(1:mx,2)
          enddo
        
          ! Third Stage
        
          call read_psi(mx,my,meqn,mbc,dx,dy,q0,maux,aux,psi,g )

          call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )
                
          XX0(1:mx,3)=psi(2:mx+1)
        
          do j=1-mbc,my+mbc
             q0(2,1:mx,j)=q(2,1:mx,j)-delt*xx0(1:mx,3)
          enddo
        
          ! Fourth Stage
        
          call read_psi(mx,my,meqn,mbc,dx,dy,q0,maux,aux,psi,g )

          call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )

          XX0(1:mx,4)=psi(2:mx+1)

!=======================================================================

        do j=1-mbc,my+mbc
                              
            q(2,1:mx,j) = q(2,1:mx,j)- delt/6.d0*( xx0(1:mx,1) + 2.d0*xx0(1:mx,2) &
               +2.d0*xx0(1:mx,3) + xx0(1:mx,4))

        enddo
     
      enddo
      
 999  continue
      
  end subroutine src2_bous
  
! =========================================================
      subroutine read_diag(mx,my,meqn,mbc,dx,dy,q,maux,aux,DL,D,DU)
! =========================================================
      use bous_module
      use geoclaw_module, only: sea_level

      implicit none
     
      integer(kind=4) mx,my,meqn,mbc,maux,mxy
      integer(kind=4) i,j,k
      real(kind=8) dx,dy

      real(kind=8)   q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      DOUBLE PRECISION   D(mx+2), DL(mx+1), DU(mx+1)                
      
      real(kind=8)  HH(1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) HHH(1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) slope(1-mbc:mx+mbc,1-mbc:my+mbc)
      
      slope = aux(1,:,:)-sea_level  
     
         do j=1-mbc,my+mbc
            do i=1-mbc,mx+mbc
               HH(i,j)=(max(0.,-slope(i,j)))**2
              HHH(i,j)=(max(0.,-slope(i,j)))**3
            enddo
         enddo
      
        D = 1.d0
        DU= 0.d0
        DL= 0.d0

        j = int(my/2)
   
        do i=1,mx
        
            if (q(1,i,j).lt.1.d-1.or.slope(i,j)>-1.d-1) then

              ! do nothing
        
            elseif (maxval(slope(i-2:i+2,j-2:j+2))  >-1d-1) then

              ! do nothing

            elseif (maxval(q(1,i-2:i+2,j-2:j+2)) < 1.d-1) then

              ! do nothing

            elseif (maxval(slope(i-1:i+1,j-1:j+1)) <-1d-1 ) then              

      ! D1 part     

              D(i+1) = 1.d0 + 2.d0*(B_param+.5d0)*HH(i,j)/dx**2 &
              -2.d0/6.d0*HHH(i,j)/(-slope(i,j))/dx**2
     
              DU(i+1)=-(B_param+.5d0)*HH(i,j)/dx**2 &
              +1.d0/6.d0*HHH(i,j)/(-slope(i+1,j))/dx**2
     
              DL(i)=-(B_param+.5d0)*HH(i,j)/dx**2 &
              +1.d0/6.d0*HHH(i,j)/(-slope(i-1,j))/dx**2
     
            endif
            
        enddo
             
      return
      end
   
!======================================================================
      subroutine read_psi(mx,my,meqn,mbc,dx,dy,q,maux,aux,psi,g )

      use geoclaw_module, only: sea_level
      use bous_module
     
      implicit none

      integer(kind=4) mx,my,meqn,mbc,maux,i,j,k,kk,iL,iR
      real   (kind=8) dx,dy,g
     
      real(kind=8)   q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) slope(1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8) psi(mx+2)

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
    
      slope = aux(1,:,:)-sea_level
     
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

      psi=0.d0

      j = int(my/2)

        do i=1,mx

           k = i+1

           topo =-slope(i,j)

           if (q(1,i,j).lt.1.d-1) then
              psi(k) = 0.d0

           elseif (slope(i,j)>-1.d-1) then
              psi(k) = 0.d0
         
           elseif (use_bous_sw_thresh.and.abs(q(1,i,j)-topo)>bous_sw_thresh*topo) then
              psi(k) = 0.d0
              !go to 994

           !elseif (-eta(i,j)>.1*q(1,i,j)) then
           !   psi(k) = 0.d0

           elseif (maxval(slope(i-2:i+2,j-2:j+2))> sw_depth ) then
            
              psi(k)=0.d0

           elseif (maxval(q(1,i-2:i+2,j-2:j+2)) < 1.d-1) then

              psi(k) = 0.d0

           else
            
              if (i==1) then
                 
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
                     -s2_H(i,j+1)+s2_H(i,j-1))/2.d0/dx/dy
              
              elseif (i==mx) then

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
                     -s2_H(i-1,j+1)+s2_H(i-1,j-1))/2.d0/dx/dy 
                          
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
                  -s2_H(i-1,j+1)+s2_H(i-1,j-1))/4.d0/dx/dy 
               
              endif
           endif
                
        enddo

        ! if a wave packet intercept the slope, then use the SWE

        do i=2,mx-1

           if (q(1,i,j)>0.01d0 .and. &
               (eta(i,j)-eta(i-1,j))*(eta(i+1,j)-eta(i,j))< 0.d0 ) then
               ! local min or max

               iL = i - abs(eta(i,j))*1./dx
               iL = max(1,iL)
               iR = i + abs(eta(i,j))*1./dx
               iR = min(mx,iR)

               if (maxval(slope(iL:iR,j))>0.d0) then
                  psi(iL:iR) = 0.d0
               endif

           endif

        end do

        ! Limit the dispersion for large troughs
        do i=2,mx-1

           if (-eta(i,j)>.5*q(1,i,j)) then

               iL = i - abs(eta(i,j))*2/dx
               iL = max(1,iL)
               iR = i + abs(eta(i,j))*2/dx
               iR = min(mx,iR)

               !psi(iL:iR) = 0.d0

           endif

        end do

        if (.not.use_bous_sw_thresh) go to 994

        do i=2,mx-1
            k= i+1

            topo = -slope(i,j)

            if (q(1,i,j)>1.d-1.and.topo>1.d-1.and. &
                q(1,i,j)-topo>bous_sw_thresh*topo) then
 
               psi(k) = 0.d0

               xaxis1: do kk=i,1,-1
                  if (eta(kk,j)-eta(kk-1,j)> 0.d0 ) then
                  !if ((q(1,kk,j)+aux(1,kk,j))> sea_level ) then
                     psi(kk) = 0.d0
                  else
                     exit xaxis1
                  endif
               end do xaxis1

               xaxis2: do kk=i,mx
                  if (eta(kk+1,j)-eta(kk,j)< 0.d0 ) then
                  !if ((q(1,kk,j)+aux(1,kk,j))> sea_level ) then
                     psi(kk) = 0.d0
                  else
                     exit xaxis2
                  endif
               end do xaxis2

            elseif (q(1,i,j)>1.d-1.and.topo>1.d-1.and. &
                q(1,i,j)-topo<-bous_sw_thresh*topo) then
 
               psi(k) = 0.d0

               xaxis11: do kk=i,1,-1
                  if (eta(kk,j)-eta(kk-1,j)< 0.d0 ) then
                     psi(kk) = 0.d0
                  else
                     exit xaxis11
                  endif
               end do xaxis11

               xaxis12: do kk=i,mx
                  if (eta(kk+1,j)-eta(kk,j)> 0.d0 ) then
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
