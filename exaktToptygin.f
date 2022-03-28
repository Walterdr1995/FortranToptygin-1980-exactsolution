      program exakt
      double precision pi,dx,dy,ymin,ymax,xmin,xmax,t
      double precision, dimension(1001) :: x,r
      double precision, dimension(301) :: y,s
      double precision, dimension(1001,301) :: exasol 
c-----------------------------------------------------------------------------
c    Declaration of basic grid parameters
      pi = 4.0*datan(1.d0)
      ymax = 7.5
      ymin = -7.5
      xmax = 75.0
      xmin = -75.0
      dx = (xmax-xmin)/dble(1000)
      dy = (ymax-ymin)/dble(300)
c.............................................................................
c.............................................................................
c     Implementation of equidistant the x-grid and a possible
c     redifinementm, in case a non-equidistant grid is needed
      do i=1,1001
      
      x(i) = xmin + dx*dble(i-1)
      
      r(i) = x(i)
c      r(i) = 1.00001*x(i) - 25.0*datan((x(i))/25.d0)
      
      enddo
c     Implementation of equidistant the y-grid and a possible
c     redifinementm, in case a non-equidistant grid is needed  
      do j=1,301
      
      y(j) = ymin + dy*dble(j-1)
      
      s(j) = y(j)
c      s(j) = 1.00001*y(j)   - 3.0*datan(y(j)/3.0)
      
      enddo
c.............................................................................      
c.............................................................................
c     Calculation of the exact Toptygin solution      
c      t = pi**2.0/6.0  ! time variable
      t= 10.0
      do j=1,301
      do i=1,1001
      
     

c     Returns the exact Toptygin solution of Toptyghin(1980) Eq. 5.9 & 5.10
      call TopTyg(r(i),s(j),t,exasol(i,j)) 
c     Return the time integral over the Toptygin solution 
c      call TimeIntgrl(r,s,t,exasol(i,j)) 

      
      enddo
      enddo
c.......................................................................
c.......................................................................
c     Here the solution we will written into a .dat file
      open(22,file='solutionx.dat')
      do j=1,301
      do i=1,1001
      if(j.eq.201) then
        write(22,'(4e18.5E4)') r(i),s(j),t,exasol(i,j)
      endif  
      enddo
      enddo
      close(22)
c.............................................................................
c-----------------------------------------------------------------------------
      end

        
      SUBROUTINE TimeIntgrl(r,s,t,loes)   
      double precision r,s,t,loes
      double precision tau, dtau, inte,add
c     Integration Parameters are defined     
      dtau = 0.001
      inte = 0.0
      tau = 0.0
c     Integration takes place 
      do while(tau.lt.t)
      call TopTyg(r,s,t-tau,add)
      inte = inte + add*dtau
      tau = tau + dtau
      enddo
c     Solution gets returned 
      loes = inte
      return
      end
      
      
      SUBROUTINE TopTyg(r,s,t,add)
      double precision r,s,t,add,q,p_0,pi
      double precision V1,V2,kappa1,kappa2,deltap
      double precision exp1,erf11,erf12,lambda1, amp
      double precision exp2,erf21,erf22,lambda2
      q = 4.0  !shock ratio
      V1    = 1.0  !upstream velocity
      V2    = V1/q !downstream velocity
      kappa1 = 1.0 !upstream diffusion coefficient
      kappa2 = kappa1/(q)**2.0 !downstream diffusion coefficient
      p_0 = 1.0 ! reference momentum
      pi = 4.0*datan(1.d0) ! pi
      x0 = 0.0 !injection spot
       
c     Some substitutions        
c      exp1  = V1*r/2.0/kappa1
c      erf11 = dabs(r)/dsqrt(4.0*kappa1*t)
c      erf12 = dsqrt(V1**2.0*t/4.0/kappa1)
c      exp2  = V2*r/2.0/kappa2
c      erf21 = dabs(r)/dsqrt(4.0*kappa2*t)
c      erf22 = dsqrt(V2**2.0*t/4.0/kappa2)
c      amp = dsqrt(4.0*kappa1/pi/V1**2.0/t)
      lambda1 =dabs(r+x0)
     .          +3.0/(V1-V2)*kappa1*(1.0+dsqrt(kappa2/kappa1))*s
      lambda2 = r+3.0/(V1-V2)*kappa2*(1.0+dsqrt(kappa1/kappa2))*s
     .          -dsqrt(kappa2/kappa1)*x0 
      deltap = 1.0/dsqrt(pi * 0.01) *dexp(- (s/0.01)**2.0)
      
      
      if(r.lt.0.0) then
            if(s.lt.0.0) then
             add = 0.0
            else
             add = deltap/2.0/p_0**2.0/dsqrt(pi*kappa1*t)
     .          *dexp(V1/2.0/kappa1*(r-x0)-V1**2.0*t/4.0/kappa1)
     .          *(dexp(-(r-x0)**2.0/4.0/kappa1/t)
     .            - dexp(-(r+x0)**2.0/4.0/kappa1/t))        
     .         + p_0**(-3.0)* 
     .         3.0*lambda1/2.0/(V1-V2)/sqrt(pi*kappa1)*dexp(-1.5*s)
     .         *t**(-1.5) *dexp(-lambda1**2.0/4.0/kappa1/t+V1*r/2.0
     .          /kappa1-V1**2.0*t/4.0/kappa1-V1*x0/2.0/kappa1) 
               !upstream solution
            endif
        
          else
            if(s.lt.0.0) then
              add = 0.0
            
            else
            
              add = p_0**(-3.0)* 
     .         3.0*lambda2/2.0/(V1-V2)/sqrt(pi*kappa2)*dexp(-1.5*s)
     .         *t**(-1.5) *dexp(-lambda2**2.0/4.0/kappa2/t+V2*r/2.0
     .          /kappa2-V2**2.0*t/4.0/kappa2-V1*x0/2.0/kappa1) 
                !downstream solution
          
            endif
        
          endif
      return
      end  
        
