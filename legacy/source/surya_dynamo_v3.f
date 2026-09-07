C     This is a 2-dimensional code for solving the kinematic dynamo
C     equation in the solar convection zone with a meridional
C     circulation. 

       implicit real*8(a-h,o-z)
       parameter (nmax=257,lmax=96)
       common u(nmax,nmax),t,it
       common/lmx/l
       common /plsincos/pl(nmax,lmax),sn(nmax)
       double precision a(nmax,nmax),b(nmax,nmax),c(nmax,nmax),
     &  d(nmax,nmax),fun(nmax),uint(nmax),
     &  e(nmax,nmax),f(nmax,nmax),vp(2*nmax,2*nmax),vq(2*nmax,2*nmax)
     &  ,a1(nmax),b1(nmax),c1(nmax),r(nmax),phi(nmax,nmax),
     &  phib(nmax,nmax),oldu(nmax,nmax),ss1(nmax),uub(nmax),
     &  uu(nmax),al(nmax,nmax),dom(nmax,nmax),ub(nmax,nmax),ra(nmax)
        double precision ab(nmax,nmax),bb(nmax,nmax),cb(nmax,nmax),
     &  db(nmax,nmax),eb(nmax,nmax),fb(nmax,nmax),eta(nmax,nmax),
     &  dror(nmax,nmax),drot(nmax,nmax),vp1(2*nmax,2*nmax),
     &  vq1(2*nmax,2*nmax),psi(2*nmax,2*nmax),etab(nmax,nmax),
     &  ss1_p(nmax),etab2(2*nmax,2*nmax),dvp(2*nmax,2*nmax),
     &  vpb(2*nmax,2*nmax),deta(2*nmax,2*nmax),uin(nmax,nmax)
        double precision aw(nmax,nmax),bw(nmax,nmax),cw(nmax,nmax),
     &  dw(nmax,nmax),ew(nmax,nmax),fw(nmax,nmax),rho(nmax,nmax),
     &  wr(nmax,nmax),wt(nmax,nmax),ft1(2*nmax,2*nmax),
     &  gt1(2*nmax,2*nmax),ft2(2*nmax,2*nmax),gt2(2*nmax,2*nmax),
     &  fr1(2*nmax,2*nmax),fr2(2*nmax,2*nmax),gr1(2*nmax,2*nmax),
     &  gr2(2*nmax,2*nmax),nuT(nmax,nmax),nuT2(2*nmax,2*nmax),
     &  rho2(2*nmax,2*nmax), nuD(nmax,nmax), lorz(nmax,nmax)
        double precision V00(2*nmax,2*nmax),H1(2*nmax,2*nmax),
     &  smoth(2*nmax), dror0(nmax,nmax),drot0(nmax,nmax)   
c     &  ,cor(2*nmax,2*nmax)
       double precision phiw(nmax,nmax),omg(nmax,nmax),
     &  aw1(nmax),bw1(nmax),cw1(nmax),rw(nmax),ww(nmax)
       external ss,erf
c        character*10 file_name, file_id
        
       n=nmax-1
       pi=4.0d0*atan(1.0d0)
C
C----------------------------------------------------------------------
C	
C	PART I. This is the ONLY part of the code in which Level I users
C	may want to make some changes.  Read Sect. 3 of the 'Guide' to
C	learn about the variables specified in this part.  
C
        irelax=1
        tmax = 10.0d0

        v0= -20.0d0
        et0=3.0d0
        et1=0.04d0
        al0=30.0d0
        omg0=2*43.28*pi 
        dt = 0.00005d0
c        dt=0.0002d0             old value
C
C----------------------------------------------------------------------
C
C	PART II.  This is the part where the alpha coefficient, 
C	diffusivity, differential rotation and meridional circulation 
C	are specified. Level II Users wishing to make changes in 
C	this part should read Sect. 4 of the 'Guide'.
C       
        pm=6.96d0
        pb=0.55d0*pm
        pk=0.73d0*pm
        pw=2.5d0*pm
        qm=pi
        dp=(pm-pb)/float(n)
        dq=-qm/float(n)
        ita=int((0.7d0-0.55d0)*pm/dp)+1
        fac=1.0d8/(3600.*24*365) 
C           
27    format(2(f13.7,2x)) 
C      Here we define a smoothing function at the 
C      boundary of convection layer    
c       open(25, file='smooth.dat',status='unknown')   
       do i = 1, 2*nmax
           p = pb + float(i-1)/float(2*n)*(pm-pb)
           smoth(i) = 0.50d0*(1.0d0+erf((p-pk)/
     &                   (0.01*pm)))
c           write(25,27) p/pm, smoth(i)
       end do

C      Density profile of the Sun
       pbc = 0.71d0*pm
       rhobc = 0.2d-6
       hbc = 0.58d0
       gmm = 5.0d0/3.0d0
c       open(24, file='density.dat',status='unknown')
       do i=1,nmax
         p = pb+float(i-1)/float(n)*(pm-pb)
         do j=1,nmax
           rho(i,j)= rhobc*(1.0d0+(((gmm-1)/gmm)*(pbc/hbc)*
     &                 ((pbc/p)-1.0d0)))**(1/(gmm-1))
       end do
c          write(24,27) p/pm, rho(i,150)*1.0d6
       end do
C      Density profile at mid-grid points       
       do i=1,2*n+1
          p = pb+float(i-1)*(pm-pb)/float(2*n)
          do j=1,2*n+1
              rho2(i,j)= rhobc*(1.0d0+(((gmm-1)/gmm)*(pbc/hbc)*
     &                 ((pbc/p)-1.0d0)))**(1/(gmm-1))
          end do
       end do      
C	
C	Here begin the do loops to calculate profiles of alpha, 
C	diffusivity and differential rotation at all the grid 
C	points. The differential rotation is written in the file 
C	'diffrot.dat'.
C
c	open(25,file='viscosity.dat',status='unknown')   
	do i = 1, nmax
	   ra(i)=pb+float(i-1)/float(n)*(pm-pb) 
          p=pb+float(i-1)/float(n)*(pm-pb)
	   do j = 1, nmax
             q=qm-float(j-1)/float(n)*qm
	      co = dcos(q)
C ALPHA PROFILE
         al(i,j) = 0.25d0*co*(1.00d0+
     &       erf((p-0.95d0*pm)/(0.030d0*pm)))*(1.00d0-erf((p
     &       -pm)/(0.030d0*pm)))
C DIFFUSIVITY PROFILES
          eta(i,j) = 0.000220d0 + (et0/2.00d0)*(1.0d0 +
     &       erf((p-0.7d0*pm)/(0.030d0*pm)))
          etab(i,j) = 0.00022d0 + (et1/2.00d0)*(1.0d0 +
     &       erf((p-0.725d0*pm)/(0.030d0*pm)))+ (et0/2.0d0)
     &       *(1.0d0+erf((p-0.975d0*pm)/(0.030d0*pm))) 
C VISCOSITY COEFFICIENT       
          nuD(i,j) = 80d0*(1.0d0 + (4.5d0 - 1.0d0)*
     &               smoth(2*i-1))*etab(i,j)    
c           nuD(i,j) = 500.0d0*etab(i,j)
C SOLAR DIFFERENTIAL ROTATION PROFILE
          dom(i,j) = 271.9d0 + 0.50d0*(1.000d0 +
     &       erf((p-0.7d0*pm)/(0.030d0*pm)))*(289.5d0 -
     &       39.4d0*co*co - 42.2d0*co*co*co*co -
     &       271.9d0) 
C density rho
c           rho(i,j) = 1.0d-6*(pm/p - 0.95d0)**1.5/2.4573d0   
c          write(25,27)q,p,dom(i,j)

         end do
c         write(25,27)p,nuD(i,125)
      end do
c        close(25)
C  
      do i= 2,n
          p=pb+float(i-1)/float(n)*(pm-pb)
          do j = 2,n
            q=qm-float(j-1)/float(n)*qm
c
            wr(i,j) = (((p+dp)**4)*rho(i+1,j)*nuD(i+1,j)-
     &           ((p-dp)**4)*rho(i-1,j)*nuD(i-1,j))/
     &            (2.0d0*dp*rho(i,j)*p**4)
            wt(i,j) = (rho(i,j+1)*nuD(i,j+1)*(dsin(q+dq))**3- 
     &            rho(i,j-1)*nuD(i,j-1)*(dsin(q-dq))**3)/
     &           (2.0d0*dq*rho(i,j)*(p**2)*(dsin(q))**3)
c            vv0(i,j) = 0.37d0
          end do
       end do
C	The do loops are closed.  We also need the diffusivity at
C	mid-points in the grid for some calculations.  This is 
C	obtained and stored now.
C
	do i=1,2*n+1
	    p=pb+float(i-1)*(pm-pb)/float(2*n)
	   do j=1,2*n+1
	      etab2(i,j) = 0.00022d0 + (et1/2.00d0)*(1.0d0 +
     &       erf((p-0.725d0*pm)/(0.030d0*pm)))+ (et0/2.0d0)
     &       *(1.0d0+erf((p-0.975d0*pm)/(0.030d0*pm))) 
C      Lambda effect coefficient     
             nuT2(i,j) = (0.01d0+(1.0d0 - 0.01d0)*
     &               smoth(i))*etab2(i,j)
c             rho2(i,j) = 1.0d-6*(pm/p - 0.95d0)**1.5/2.4573d0
	   end do
      end do
C	
C	Derivatives of differential rotation are obtained and stored
C	now for future use.
C
      do i=2,n
	   do j=2,n
	       dror(i,j)=(dom(i+1,j)-dom(i-1,j))/(2.0d0*dp)  
	       drot(i,j)=(dom(i,j+1)-dom(i,j-1))/(2.0d0*dq)
         end do
      end do
C
C MERIDIONAL CIRCULATION. Note that this is calculated both at the 
C   grid-points and mid-points. First the stream function psi(i,j)is 
C   calculated. It is written in the file 'psi.dat'.	
C
       beta1=1.5d0
       beta2=1.3d0
       pp=0.73d0*pm
       del=2.0000001d0
       gm=3.47d0
       p0=(pm-pb)/3.50d0        
c	open(82,file='psi.dat',status='unknown')     
	do i=1,2*n+1
          p=pb+float(i-1)/float(2*n)*(pm-pb)
	   do j=1,2*n+1
            q=qm-float(j-1)/float(2*n)*qm
C
            if(q.le.pi/2.0d0) then
              exq0=dexp(-beta1*q**del)
              exqm=dexp(beta2*(q-pi/2.0d0))
              gau=dexp(-1.05*((p-p0)/gm)**2)
c
              psi(i,j)=(p-pp)*dsin(pi*(p-pp)/(pm-pp))*
     &             (1.00d0-exq0)*(1.00d0-exqm)*gau
            end if 
C
            if(q.gt.pi/2.0d0) then
              exq0=dexp(-beta1*(pi-q)**del)
              exqm=dexp(beta2*(pi/2.0d0-q))
              gau=dexp(-1.05*((p-p0)/gm)**2)
              psi(i,j)=-(p-pp)*dsin(pi*(p-pp)/(pm-pp))*
     &             (1.00d0-exq0)*(1.00d0-exqm)*gau
            end if

            if (p.lt.pp) then
              psi(i,j)=0.0d0
            end if
C
c            write(82,34)q,p,psi(i,j)
 34         format(3(f13.5,1x))
          end do
       end do
c       close(82)

C	Now components of velocity are calculated at grid-points and 
C	mid-points from the stream function psi(i,j).
C
c       open(87,file='v_theta.dat',status='unknown')
c
	do i=2,2*n
          p=pb+float(i-1)*(pm-pb)/float(2*n)
	   do j=2,2*n
             q=qm-float(j-1)*qm/float(2*n)
             vp1(i,j)=2.4573d0*v0*(psi(i,j+1)-psi(i,j-1))/
     &            (p**2*dsin(q)*dq*(pm/p-0.95d0)**1.5)
             vq1(i,j)=-2.4573d0*v0*(psi(i+1,j)-psi(i-1,j))/
     &            (p*dsin(q)*dp*(pm/p-0.95d0)**1.5)
c
c             write(87,34)q,p,vq1(i,j)
c
             deta(i,j)=(etab2(i+1,j)-etab2(i-1,j))/(dp)
             vpb(i,j)=(vp1(i,j)-deta(i,j))*dt/(2.d0*p*dp)
             vp(i,j)=vp1(i,j)*dt/(2.d0*p*dp)
             vq(i,j)=vq1(i,j)*dt/(2.d0*p*dsin(q)*dq)
          end do
       end do
c       close(87)
C  computing the required functions needed to compute the coefficients
       do i=1,2*n+1
          p = pb+float(i-1)*(pm-pb)/float(2*n)
          pr = p/pm
          do j = 1,2*n+1
            q = qm-float(j-1)*qm/float(2*n)
c            if (p.ge.pk) then
c             coriolis number function
              c0 = 110.80d0*dexp(-602.3235*pr**3+1496.1566*pr**2 -
     &             1245.933*pr + 344.6711)
c      Here,   J0,J1 = z0,z1 and I0,I1 = y0,y1
              z0 = (1.0d0/(2.0d0*c0**4))*(9.0d0-((2.0d0*c0**2)/
     &              (1.0d0+c0**2))-((c0**2+9.0d0)*datan(c0)/c0))
              z1=(1.0d0/(2.0d0*c0**4))*(45.0d0+c0**2-((4.0d0*c0**2)/
     &         (1.0d0+c0**2))+((c0**4-12.0d0*c0**2-45)*datan(c0)/c0))
              y0 = (1.0d0/(4.0d0*c0**4))*(-19.0d0-(5.0d0/
     &           (1.0d0+c0**2))+((3.0d0*c0**2+24.0d0)*datan(c0)/c0))
              y1 = (3.0d0/(4.0d0*c0**4))*(-15.0d0+(2.0d0*(c0**2)/
     &             (1.0d0+c0**2))+((3.0d0*c0**2+15.0d0)*datan(c0)/c0))
                
              V00(i,j) = (10.0d0**2)*smoth(i)*(z0+4.0d0*y0)
              H1(i,j) = (30.0d0**2)*smoth(i)*(z1+4.0d0*y1)
c              
c            * , p, q, c0, (V00(i,j)-H1(i,j)*dcos(q)**2),
c     &                    (dsin(q)**2)*H1(i,j)
          end do
       end do

       do i=1,2*n+1
          p = pb+float(i-1)*(pm-pb)/float(2*n)
          do j = 1,2*n+1
             q = qm-float(j-1)*qm/float(2*n)
             ft1(i,j) = vq1(i,j)/(p*(dsin(q))**2)
             gt1(i,j) = (dsin(q))**2
             ft2(i,j) = 1.0d0/(rho2(i,j)*(p**2)*(dsin(q)**3))
             gt2(i,j) = rho2(i,j)*nuT2(i,j)*dcos(q)*
     &                  (dsin(q)**4)*H1(i,j)
             fr1(i,j) = vp1(i,j)/p**2
             gr1(i,j) = p**2
             fr2(i,j) = 1.0d0/(rho2(i,j)*p**4)
             gr2(i,j) = rho2(i,j)*nuT2(i,j)*(p**3)*(V00(i,j)-
     &                    H1(i,j)*dcos(q)**2)
c             print* , gr2(i,j), gt2(i,j)
          end do
       end do
C
c------------------------------------------------------------
C
C	PART II-A.  Legendre polynomials are calculated here with the
C	help of a subroutine and stored for future use.
C
       do l=1,96
         cl=1.0d0/(1.0d0+float(l)/float(l+1)*(pm/pw)**(2*l+1))
         ss1_p(l)=float(2*l+1)/(float(l)*pm)*(cl-
     &            float(l)/float(l+1)*(1.0d0-cl))
          do j=1,n+1
             qm=4.0d0*atan(1.0d0)
             q=qm+float(j-1)*dq
             x=dcos(q)
             pl(j,l) = plgndr(l,1,x)
             sn(j) = dsin(q)
          end do
       end do
C
C---------------------d---------------------------------------------
C
C	PART III. This is the part where initial values of the 
C	variables u(i,j) and ub(i,j) are chosen, before beginning 
C	the time advancement.

42    format(5(f13.7,1x))
	if(irelax.eq.0) then
	  do i = 1,n+1
            p=pb+float(i-1)/float(n)*(pm-pb)
	  do j = 1,n+1
            q=qm-float(j-1)/float(n)*qm
	  u(i,j) = 0.0d0
	  ub(i,j) = 1.0*dsin(2.0d0*q)*dsin(pi
     &      *((p-pb)/(pm-pb)))
	  end do
	  end do
	  go to 10
	end if

       open(12,file='init.dat',status='unknown')
       do i=1,n+1
          p=pb+(pm-pb)*float(i-1)/float(n)
         do j=1,n+1
          read(12,42) q,p,uin(i,j),ub(i,j)
           if(j.eq.n+1) then
            u(i,j)=0.0d0
           else
            u(i,j)=uin(i,j)/(p*dsin(q))
           end if
         end do
       end do
       close(12)
C      Initial value for differential rotation
       do i=1,nmax
          do j=1,nmax
              omg(i,j)= omg0
          end do
       end do
44     format(4(f13.5,1x))
C
C----------------------------------------------------------------
C
C	PART IV. For advancing in time, we need some matrix 
C	elements arising out of the difference schemes used. We 
C	calculate and store the matrix elements now.  Look at 
C	Sect. 5 of the 'Guide' for a discussion about these 
C	matrix elements.  Only Level III Users need to bother
C	about the matrix elements.
C
10    do i=2,2*n
         do j=2,2*n
          dvp(i,j)=(vp1(i+1,j)-vp1(i-1,j))/(dp)
         end do
      end do  
C     
C     coefficients for the numerical sceme of rotation profile
      do i= 2,n
         p = pb+float(i-1)*(pm-pb)/float(n)
         do j= 2,n
            q = qm-float(j-1)*qm/float(n)
            aw(i,j) = dt*wt(i,j)/(4.0d0*dq) - dt*nuD(i,j)/
     &         (2.0d0*(p*dq)**2) - (dt*ft1(2*i-1,2*j-1)*
     &         gt1(2*i-1,2*j-2)/(4.0d0*dq))*(1.0d0 + 
     &         (dt*ft1(2*i-1,2*j-2)*gt1(2*i-1,2*j-3))/(2.0d0*dq))
     &         -(dt*ft2(2*i-1,2*j-1)* gt2(2*i-1,2*j-2)/(4.0d0*dq))*
     &         (1.0d0 +(dt*ft2(2*i-1,2*j-2)*gt2(2*i-1,2*j-3))/
     &          (2.0d0*dq))
            bw(i,j) = (dt*nuD(i,j)/(1.0d0*(p*dq)**2)) + 
     &         (dt*ft1(2*i-1,2*j-1)*gt1(2*i-1,2*j)/(4.0d0*dq))*
     &         (1.0d0+(dt*ft1(2*i-1,2*j)*gt1(2*i-1,2*j-1))/
     &         (2.0d0*dq)) - (dt*ft1(2*i-1,2*j-1)*gt1(2*i-1,2*j-2)/
     &         (4.0d0*dq))*(1.0d0-(dt*ft1(2*i-1,2*j-2)
     &         *gt1(2*i-1,2*j-1))/(2.0d0*dq))  +
     &         (dt*ft2(2*i-1,2*j-1)*gt2(2*i-1,2*j)/(4.0d0*dq))*
     &         (1.0d0+(dt*ft2(2*i-1,2*j)*gt2(2*i-1,2*j-1))/
     &         (2.0d0*dq)) - (dt*ft2(2*i-1,2*j-1)*gt2(2*i-1,2*j-2)/
     &         (4.0d0*dq))*(1.0d0-(dt*ft2(2*i-1,2*j-2)
     &         *gt2(2*i-1,2*j-1))/(2.0d0*dq))
            cw(i,j) = -(dt*wt(i,j)/(4.0d0*dq)) - (dt*nuD(i,j)/
     &         (2.0d0*(p*dq)**2)) + (dt*ft1(2*i-1,2*j-1)*
     &         gt1(2*i-1,2*j)/(4.0d0*dq))*(1.0d0 - 
     &         (dt*ft1(2*i-1,2*j)*gt1(2*i-1,2*j+1))/(2.0d0*dq)) +
     &         (dt*ft2(2*i-1,2*j-1)*gt2(2*i-1,2*j)/(4.0d0*dq))*(1.0d0
     &         -(dt*ft2(2*i-1,2*j)*gt2(2*i-1,2*j+1))/(2.0d0*dq))
            dw(i,j) = (dt*wr(i,j)/(4.0d0*dp)) - (dt*nuD(i,j)/
     &          (2.0d0*dp**2))- (dt*fr1(2*i-1,2*j-1)*
     &         gr1(2*i-2,2*j-1)/(4.0d0*dp))*(1.0d0 +
     &         (dt*fr1(2*i-2,2*j-1)*gr1(2*i-3,2*j-1))/
     &          (2.0d0*dp))- (dt*fr2(2*i-1,2*j-1)*
     &         gr2(2*i-2,2*j-1)/(4.0d0*dp))*(1.0d0 +
     &         (dt*fr2(2*i-2,2*j-1)*gr2(2*i-3,2*j-1))/(2.0d0*dp))
            ew(i,j) = (dt*nuD(i,j)/(1.0d0*dp**2)) +
     &         (dt*fr1(2*i-1,2*j-1)*gr1(2*i,2*j-1)/(4.0d0*dp))*
     &         (1.0d0 + (dt*fr1(2*i,2*j-1)*gr1(2*i-1,2*j-1))/
     &         (2.0d0*dp))+(dt*fr2(2*i-1,2*j-1)*gr2(2*i,2*j-1)/
     &         (4.0d0*dp))*(1.0d0+(dt*fr2(2*i,2*j-1)*
     &         gr2(2*i-1,2*j-1))/(2.0d0*dp))-
     &         (dt*fr1(2*i-1,2*j-1)*gr1(2*i-3,2*j-1)/(4.0d0*dp))*
     &         (1.0d0 - (dt*fr1(2*i-3,2*j-1)*gr1(2*i-1,2*j-1))/
     &         (2.0d0*dp))-(dt*fr2(2*i-1,2*j-1)*gr2(2*i-3,2*j-1)/
     &         (4.0d0*dp))*(1.0d0 - (dt*fr2(2*i-3,2*j-1)*
     &         gr2(2*i-1,2*j-1))/(2.0d0*dp))
            fw(i,j) = -(dt*wr(i,j)/(4.0d0*dp)) - (dt*nuD(i,j)/
     &         (2.0d0*dp**2)) + (dt*fr1(2*i-1,2*j-1)*gr1(2*i,2*j-1)/
     &         (4.0d0*dp))*(1.0d0-(dt*fr1(2*i,2*j-1)*
     &         gr1(2*i+1,2*j-1))/(2.0d0*dp)) + (dt*fr2(2*i-1,2*j-1)*
     &         gr2(2*i,2*j-1)/(4.0d0*dp))*(1.0d0-(dt*
     &         fr2(2*i,2*j-1)*gr2(2*i+1,2*j-1))/(2.0d0*dp))
         end do
      end do
C
C     Coefficients for numerical scheme of magnetic field
      do i=2,n
           p=pb+float(i-1)/float(n)*(pm-pb)
           do j=2,n
            q=qm-float(j-1)/float(n)*qm
            a(i,j)=-(eta(i,j)*dt/(2.0d0*(p*dq)**2)-eta(i,j)*dt/
     &         (4.0d0*dtan(q)*p*p*dq)+
     &         vq(2*i-1,2*j-1)*dsin(q-dq/
     &         2.0d0)*(1.0d0+vq(2*i-1,2*j-2)*dsin(q-dq))/2.0d0)
            b(i,j)=-(-eta(i,j)*dt/((p*dq)**2)-
     &         0.0*eta(i,j)*dt/(2.0d0*dtan(q)*p*p*dq)-
     &         eta(i,j)*dt/(4.0d0*(p*dsin(q))**2)-vq(
     &         2*i-1,2*j-1)*(dsin(q+dq/2.0d0)*(1.0d0+vq(2*i-1,2*j)
     &         *dsin(q))-dsin(q-dq/2.0d0)*(1.0d0-
     &         vq(2*i-1,2*j-2)*dsin(q)))/2.0d0)
            c(i,j)=-(eta(i,j)*dt/(2.0d0*(p*dq)**2)+
     &         eta(i,j)*dt/(4.0d0*dtan(q)*p*p*
     &         dq)-vq(2*i-1,2*j-1)*dsin(q+dq/2.0d0)*(1.0d0
     &         -vq(2*i-1,2*j)*dsin(q+dq))/2.0d0)         
            d(i,j)=-(eta(i,j)*dt/(2.0d0*dp**2)-eta(i,j)*dt/
     &         (2.0*p*dp)+vp(2*i-1,2*j-1)*(p-dp/2.0d0)*
     &         (1.0d0+vp(2*i-2,2*j-1)*(p-dp))/2.0d0-
     &         0.0d0*eta(i,j)*dt/(p*dp))
            e(i,j)=-(-eta(i,j)*dt/(dp**2)-
     &         0.0d0*eta(i,j)*dt/(p*dp)-eta(i,j)*dt/(4.0d0*
     &         (p*dsin(q))**2)-vp(2*i-1,2*j-1)*(dp+p*((p+dp/
     &         2.0d0)*vp(2*i,2*j-1)+(p-dp/2.0d0)*
     &         vp(2*i-2,2*j-1)))/2.0d0)
            f(i,j)=-(eta(i,j)*dt/(2.0d0*dp**2)+eta(i,j)*
     &         dt/(2.0*p*dp)-vp(2*i-1,2*j-1)*(p+dp/2.0d0)*
     &         (1.0d0-vp(2*i,2*j-1)*(p+dp))/2.0d0)


            ab(i,j)=-(etab(i,j)*dt/(2.0d0*(p*dq)**2)-etab(i,j)*dt/
     &         (4.0d0*dtan(q)*p*p*dq)+
     &         vq(2*i-1,2*j-2)*dsin(q-dq/
     &         2.0d0)*(1.0d0+vq(2*i-1,2*j-3)*dsin(q-dq))/2.0d0)
            bb(i,j)=-(-etab(i,j)*dt/((p*dq)**2)-
     &         etab(i,j)*dt/(4.0d0*(p*dsin(q))**2)-(vq(
     &         2*i-1,2*j)*dsin(q+dq/2.0d0)*(1.0d0+vq(2*i-1,2*j-1)*
     &         dsin(q))-vq(2*i-1,2*j-2)*dsin(
     &         q-dq/2.0d0)*(1.0d0-vq(2*i-1,2*j-1)*dsin(q)))/2.0d0)
            cb(i,j)=-(etab(i,j)*dt/(2.0d0*(p*dq)**2)+
     &         etab(i,j)*dt/(4.0d0*dtan(q)*p*p*
     &         dq)-vq(2*i-1,2*j)*dsin(q+dq/2.0d0)*(1.0d0-
     &         vq(2*i-1,2*j+1)*dsin(q+dq))/2.0d0)
            db(i,j)=-(etab(i,j)*dt/(2.0d0*dp**2)-etab(i,j)*dt/
     &         (2.d0*p*dp)+vpb(2*i-1,2*j-1)*(p-dp/2.0d0)*
     &         (1.0d0 +vpb(2*i-2,2*j-1)*(p-dp))/(2.0d0))
            eb(i,j)=-(-etab(i,j)*dt/(dp**2)-
     &         0.0d0*etab(i,j)*dt/(p*dp)-etab(i,j)*dt/(4.0d0*
     &         (p*dsin(q))**2)-dvp(2*i-1,2*j-1)*dt/2.0d0-
     &         vpb(2*i-1,2*j-1)*(dp+vpb(2*i,2*j-1)*p*(p+dp/
     &         2.0d0)+(p-dp/2.0d0)*p*vpb(2*i-2,2*j-1))/2.0d0)
            fb(i,j)=-(etab(i,j)*dt/(2.0d0*dp**2)+
     &         etab(i,j)*dt/(2.d0*p*dp)-vpb(2*i-1,2*j-1)*
     &         (p+dp/2.0d0)*(1.0d0 -vpb(2*i,2*j-1)*(p+dp))/2.0d0)
c              print* ,aw(i,j),a(i,j)
          end do
       end do
57     format(8(e15.7, 3x))
C       open(51,file='coeff.dat',status='unknown')
C           do i=1,nmax
C             p=pb+float(i-1)/float(n)*(pm-pb)
C            do j=1,nmax
C             q=qm-float(j-1)/float(n)*qm
C             write(51,57) p,q,aw(i,j),bw(i,j),cw(i,j),
C     &           dw(i,j),ew(i,j),fw(i,j)
C            end do
C           end do
C       close(51)
C
C-------------------------------------------------------------------------
C
C	PART V. Now we come to the central part of the programme where the
C	time advancement takes place. 'k' is the counter to keep tab on
C	the time.  In each step of the do loop 'do k=1,kend', the magnetic
C	fields are advanced through one time step 'dt'.
       kend = tmax/dt

       t=0.0d0
 
       do k=1,kend

C PART V-A. CALCULATING TIME ADVANCED MAGNETIC FIELDS AT INTERIOR GRID 
C POINTS. Sect. 5 of the 'Guide' explains how magnetic fields are 
C advanced by solving tridiagonal matrices with the subroutine 'tridag'              

C FIRST HALF-STEP OF TIME ADVANCEMENT
C     Here we calculate the Lorentz force source term
       do i=2,n
          p=pb+float(i-1)/float(n)*(pm-pb)
          do j=2,n
             q=qm-float(j-1)/float(n)*qm
             mu0 = 4.0d0*pi
             lorz(i,j) = 1.0d0/(mu0*rho(i,j)*dsin(q)*(p**2))*
     &          (((dsin(q+dq)*u(i,j+1)-dsin(q-dq)*u(i,j-1))/
     &          (2.0d0*dq))*(((p+dp)*ub(i+1,j) - (p-dp)*ub(i-1,j))
     &          /(2.0d0*dp))-(((p+dp)*u(i+1,j)-(p-dp)*u(i-1,j))/
     &          (2.0d0*dp))*((dsin(q+dq)*ub(i,j+1)-dsin(q-dq)*
     &           ub(i,j-1))/(2.0d0*dq)))
          end do
       end do
       
   
       do j=2,n
         do i=2,n
           phiw(i,j)= -dw(i,j)*omg(i-1,j)+(1-ew(i,j))*
     &               omg(i,j)-fw(i,j)*omg(i+1,j)+lorz(i,j)*dt/2.0d0
         end do
       end do
   
       do i=2,n
         do j=2,n
           aw1(j-1)= aw(i,j)
           bw1(j-1)= bw(i,j)+1.0d0
           cw1(j-1)= cw(i,j)
           rw(j-1)= phiw(i,j)
         end do
         rw(1)=phiw(i,2)-aw1(1)*omg(i,1)
         rw(n-1)=phiw(i,n)-cw1(n-1)*omg(i,n+1)
         call tridag(aw1,bw1,cw1,rw,ww,n-1)
         do j=2,n
           omg(i,j) = ww(j-1)
           phiw(i,j)=-phiw(i,j)+2.0d0*ww(j-1)
         end do
       end do
c       print* , 'phiw1', phiw(200,150), 'omg1', omg(200,150) 

       do i=2,n
         do j=2,n
           dror(i,j)=(omg(i+1,j)-omg(i-1,j))/(2.0d0*dp)  
           drot(i,j)=(omg(i,j+1)-omg(i,j-1))/(2.0d0*dq)
         end do
       end do

       do i=2,n
          do j=2,n
             phi(i,j)=-d(i,j)*u(i-1,j)+(-e(i,j)+1.0d0)*u(i,j)
     &        -f(i,j)*u(i+1,j)+al0*al(i,j)*ub(i,j)*dt/2.0d0
           end do
       end do

       do i=2,n
           do j=2,n
              a1(j-1)=a(i,j)
              b1(j-1)=b(i,j)+1.0d0
              c1(j-1)=c(i,j)
              r(j-1)=phi(i,j)
           end do
           r(1)=phi(i,2)-a1(1)*u(i,1)
           r(n-1)=phi(i,n)-c1(n-1)*u(i,n+1)
           call tridag(a1,b1,c1,r,uu,n-1)
           do j=2,n
              u(i,j) = uu(j-1)
              phi(i,j)=-phi(i,j)+2.0d0*uu(j-1)
           end do
       end do

       do i=2,n
              p=pb+float(i-1)/float(n)*(pm-pb)
              do j=2,n
                q=qm-float(j-1)/float(n)*qm
                  br = (u(i,j+1)*dsin(q+dq)-u(i,j-1)*dsin(q-dq))*dt
     &       /(4.0d0*dq)
                bt=(u(i-1,j)*(p-dp)-u(i+1,j)*(p+dp))*dsin(q)*dt/
     &        (4.0d0*p*dp)
                phib(i,j)=-db(i,j)*ub(i-1,j)+(-eb(i,j)+1.0d0)*ub(i,j)
     &        -fb(i,j)*ub(i+1,j)+br*dror(i,j)+bt*drot(i,j)
                end do
       end do

       do i=2,n
          do j=2,n
              a1(j-1)=ab(i,j)
              b1(j-1)=bb(i,j)+1.0d0
              c1(j-1)=cb(i,j)
              r(j-1)=phib(i,j)
          end do
          r(1)=phib(i,2)-a1(1)*ub(i,1)
          r(n-1)=phib(i,n)-c1(n-1)*ub(i,n+1)
          call tridag(a1,b1,c1,r,uub,n-1)
          do j=2,n
              ub(i,j) = uub(j-1)
              phib(i,j)=-phib(i,j)+2.0d0*uub(j-1)
              phi(i,j)=phi(i,j)+al0*al(i,j)*ub(i,j)*dt/2.0d0
          end do
       end do

C
C   SECOND HALF-STEP OF TIME ADVANCEMENT 
C
       do i=2,n
          p=pb+float(i-1)/float(n)*(pm-pb)
            do j=2,n
              q=qm-float(j-1)/float(n)*qm
              mu0 = 4.0d0*pi
              lorz(i,j) = 1.0d0/(mu0*rho(i,j)*dsin(q)*(p**2))*
     &          (((dsin(q+dq)*u(i,j+1)-dsin(q-dq)*u(i,j-1))/
     &          (2.0d0*dq))*(((p+dp)*ub(i+1,j) - (p-dp)*ub(i-1,j))
     &          /(2.0d0*dp))-(((p+dp)*u(i+1,j)-(p-dp)*u(i-1,j))/
     &          (2.0d0*dp))*((dsin(q+dq)*ub(i,j+1)-dsin(q-dq)*
     &           ub(i,j-1))/(2.0d0*dq)))
              phiw(i,j)=phiw(i,j)+lorz(i,j)*dt/2.0d0
           end do
       end do

   
c       do i=2,n
c           do j=2,n
c             phiw(i,j)= -aw(i,j)*omg(i,j-1)+(1-bw(i,j))*
c     &         omg(i,j)-cw(i,j)*omg(i,j+1)+lorz(i,j)*dt/2.0d0
c           end do
c       end do

       do j=2,n
           do i=2,n
             aw1(i-1)=dw(i,j)
             bw1(i-1)=ew(i,j)+1.0d0
             cw1(i-1)=fw(i,j)
             rw(i-1)=phiw(i,j)
           end do
           rw(1)=phiw(2,j)-aw1(1)*omg(1,j)
           rw(n-1)=phiw(n,j)-cw1(n-1)*omg(n+1,j)
           call tridag(aw1,bw1,cw1,rw,ww,n-1)
           do i=2,n
             omg(i,j)=ww(i-1)
           end do
       end do  
c       print* ,'t', t , 'omg2', omg(200,150)   
       
       do i=2,n
        do j=2,n
          dror(i,j)=(omg(i+1,j)-omg(i-1,j))/(2.0d0*dp)  
          drot(i,j)=(omg(i,j+1)-omg(i,j-1))/(2.0d0*dq)
        end do
       end do
c 
       do i=2,n
          p=pb+float(i-1)/float(n)*(pm-pb)
          do j=2,n
           q=qm-float(j-1)/float(n)*qm
           br = (u(i,j+1)*dsin(q+dq)-u(i,j-1)*dsin(q-dq))*dt/
     &      (4.0d0*dq)
           bt=(u(i-1,j)*(p-dp)-u(i+1,j)*(p+dp))*dsin(q)*dt/
     &      (4.0d0*p*dp)
           phib(i,j)=phib(i,j)+br*dror(i,j)+bt*drot(i,j)        
          end do
       end do
 
       do j=2,n
           do i=2,n
             a1(i-1)=d(i,j)
             b1(i-1)=e(i,j)+1.0d0
             c1(i-1)=f(i,j)
             r(i-1)=phi(i,j)
            end do
             r(1)=phi(2,j)-a1(1)*u(1,j)
             r(n-1)=phi(n,j)-c1(n-1)*u(n+1,j)
             call tridag(a1,b1,c1,r,uu,n-1)
          do i=2,n
             u(i,j)=uu(i-1)
          end do
       end do
C
C
       do j=2,n
           do i=2,n
             a1(i-1)=db(i,j)
             b1(i-1)=eb(i,j)+1.0d0
             c1(i-1)=fb(i,j)
             r(i-1)=phib(i,j)
            end do
             r(1)=phib(2,j)-a1(1)*ub(1,j)
             r(n-1)=phib(n,j)-c1(n-1)*ub(n+1,j)
             call tridag(a1,b1,c1,r,uub,n-1)
            do i=2,n
             ub(i,j)=uub(i-1)
            end do
       end do

C      BOUNDARY CONDITIONS FOR DIFFERENTIAL ROTATION
C      Boundary Conditions at the bottom and top
       do j=1,nmax
          omg(1,j) = omg0
          omg(nmax,j) = omg(n,j)
       end do
C         Boundary condition at poles
       do i=2,n
          omg(i,1) = omg(i,2)
          omg(i,nmax) = omg(i,n)
       end do
C
C      Printing omega every 200th step
c       if (k/1000*1000.eq.k) then
c              print*, t, omg(200,150)
c       end if
C       
C      Storing rotation speed at a fix radius at constant time interval 
       if (k/200*200.eq.k) then
       open(71,file='omg_lat.dat',status='unknown',access='append')
         ir = idint(1+((0.85*pm)-pb)*float(n)/(pm-pb))
         do j=1,nmax
           q = qm-float(j-1)/float(n)*qm  
           write(71,37) q,t,omg(ir,j)
         end do
         close(71)
       end if 
 
C PART V-B. BOUNDARY CONDITIONS. These are discussed in Sect. 2.1 of .
C the 'Guide'. 
C Boundary condition at bottom
       do j=1,n+1
            u(1,j)=0.0d0
            ub(1,j)=0.0d0
       end do
C Boundary condition at top surface
       do nup=1,11
          do l=1,96
             ss1(l)=ss1_p(l)*ss(n,dq)
          end do
          do j=1,n+1
             call coeff(n,ss1,dq,j,cf)
             u(n+1,j)=u(n,j)+cf*dp
          end do
       end do
       do j=1,n+1
            ub(n+1,j)=0.0d0
       end do
C Boundary conditions at the poles
       do i=2,n+1
            u(i,n+1)=0.0d0
            ub(i,n+1)=0.0d0
            u(i,1)=0.0d0
            ub(i,1)=0.0d0
       end do
C
C PART V-C.  INCORPORATING THE MAGNETIC BUOYANCY by radial transport
C of toroidal field exceeding a specified critical value:
!
47       format(3(d15.6,1x))
c       if(k/40*40.eq.k)then
c          nt = n - 6.
c       do i = int(n/3.), int(n/2.)
c          if (ra(i).gt.0.71*pm)then
c            do j=2,n
c             if (ABS(ub(i,j)).GT.0.8d0)then
c              jer=j
c              ier=i
c              qjer=qm-float(jer-1)*qm/float(n)	
c              ber=ub(ier,jer)
c              bold=ub(nt,jer)
c              ub(nt,jer)=((ra(ier)/ra(nt))*0.5d0*ber) + bold
c              ub(ier,jer)=0.5d0*ber
c               open(25,file='ber.dat',status='unknown',
c     &              access='append')
c               write(25,47) t,qjer,ber
c               close(25)
c             end if
c            end do
c          end if
c         end do
c        end if
C
C PART V-D. WRITING IN THE FILES 'rad.dat' and 'butbot.dat' after 
C certain intervals of time
C
37     format(3(f13.7,1x))
c       if(t.gt.0.0d0) then
c           p20 = pi/20.
c           istep = t/p20
c           tdiff = t - float(istep)*p20
c           if(tdiff.ge.dt) go to 55
c       open(17,file='rad.dat',status='unknown',access='append') 
c       open(18,file='butbot.dat',status='unknown',access='append') 
c      open(19,file='apot.dat',status='unknown',access='append')
c          do j = 2, n
c            q=qm-float(j-1)/float(n)*qm
c            br = -(u(n-1,j-1)*dsin(q-dq) - u(n-1,j+1)
c     &         *dsin(q+dq))/(2.0d0*dq*dsin(q))/pm
c            write(17,37) t, q, br
c            write(18,37) t, q, ub(90,j)
c 	    write(19,37) t, q, u(122,j)
c          end do
c         close(17)
c         close(18)
c 	  close(19)
c55    end if	
C
C WRITING IN THE FILE 'run.dat' after every 40 time steps
C
49     format(2(f13.6,1x))
c        if (k/40*40.eq.k) then
c          open(95,file='series_b.dat',status='unknown',access='append')
c           write(95,49)(t*fac),abs(ub(ita,150))
c           close(95)
c        end if

       if (k/50*50.eq.k) then
         print* , t, omg(200,150)
         open(91,file='diffrot.dat',status='unknown')
         do i=1,nmax
          p=pb+float(i-1)/float(n)*(pm-pb)
           do j=1,nmax
             q=qm-float(j-1)/float(n)*qm
             write(91,44) q,p,omg(i,j),t
           end do
         end do
         close(91)
       end if


       t=t+dt    ! Time advancing
       end do
C
C  Time advancing loop ends here.
C
c       do i=1,nmax
c          do j=1,nmax
c              print* , omg(i,j)
c          end do
c       end do
       open(91,file='diffrot.dat',status='unknown')
       do i=1,nmax
         p=pb+float(i-1)/float(n)*(pm-pb)
         do j=1,nmax
          q=qm-float(j-1)/float(n)*qm
          write(91,37) q,p,omg(i,j)
         end do
       end do
       close(91)

C---------------------------------------------------------------
C	PART VI. Now the final state is written in the file 
C	'final.dat'
C
43      format(6(f13.7,1x))
        open(12,file='final.dat',status='unknown')
        do i=1,n+1
          p=pb+float(i-1)/float(n)*(pm-pb)
          do j=1,n+1
            q=qm-float(j-1)/float(n)*qm
            write(12,42) q,p,omg(i,j),p*dsin(q)*u(i,j),ub(i,j)
           end do
        end do
        close(12)
C       
       stop
       end
 
C
C     MAIN PROGRAM ENDS HERE.
C -----------------------------------------------------
C
C     SUBROUTINE FOR INVERTING TRIDIAGONAL MATRIX :
         subroutine tridag(a,b,c,r,u,n)
         implicit real*8 (a-h,o-z)
         parameter (nmax=257)
         dimension gam(nmax),a(nmax),b(nmax),c(nmax),r(nmax),u(nmax)
         if(b(1).eq.0.d0) pause
         bet=b(1)
         u(1)=r(1)/bet
         do 11 j=2,n
           gam(j)=c(j-1)/bet
           bet=b(j)-a(j)*gam(j)
           if(bet.eq.0.d0) pause
           u(j)=(r(j)-a(j)*u(j-1))/bet
11         continue
         do 12 j=n-1,1,-1
           u(j)=u(j)-gam(j+1)*u(j+1)
12         continue
         return
         end
 
 
C
C     THREE SUBROUTINES NEEDED FOR UPPER BOUNDARY CONDITION :
C       The subroutine plgndr calculates Legendre polynomials.
         subroutine coeff(n,ss1,dq,j,cf)
         implicit real*8(a-h,o-z)
         external plgndr
         parameter (nmax=257,lmax=96)
         dimension ss1(nmax)
         common u(nmax,nmax),t,it
         common/lmx/l
         common /plsincos/pl(nmax,lmax),sn(nmax)
         cf=0.0d0
         do l=1,96
         cf=cf+ss1(l)*pl(j,l)
         end do
         return
         end
 
 
 
         function ss(n,dq)
         implicit real*8(a-h,o-z)
         external plgndr
         parameter (nmax=257,lmax=96)
         dimension y(nmax)
         common u(nmax,nmax),t,it
         common/lmx/l
         common /plsincos/pl(nmax,lmax),sn(nmax)
         do j=1,n+1
         y(j)=pl(j,l)*u(n+1,j)*sn(j)
         end do
         term1=0.0d0
         do k1=2,n,2
         term1=term1+y(k1)
         end do
         term2=0.0d0
         do k2=3,n-1,2
         term2=term2+y(k2)
         end do
         ss=dq/3.0d0*(y(1)+4.0d0*term1+2.0d0*term2+y(n+1))
         return
         end
 
 
         function plgndr(l,m,x)
         implicit real*8(a-h,o-z)
         if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.0d0)pause 'bad arguments'
         pmm=1.0d0
         if(m.gt.0) then
           somx2=dsqrt((1.0d0-x)*(1.0d0+x))
           fact=1.0d0
           do 11 i=1,m
              pmm=-pmm*fact*somx2
              fact=fact+2.0d0
11          continue
         end if
         if(l.eq.m) then
           plgndr=pmm
         else
           pmmp1=x*(2*m+1)*pmm
           if(l.eq.m+1) then
              plgndr=pmmp1
           else
              do 12 ll=m+2,l
                  pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
                  pmm=pmmp1
                  pmmp1=pll
12         continue
            plgndr=pll
           end if
         end if
         return
         end
 
c
       function erf (x)
       implicit real*8 (a-h,o-z)
       data p/0.3275911d0/,a1/0.254829592d0/,a2/-0.284496736d0/,
     &     a3/1.421413741d0/,a4/-1.453152027d0/,a5/1.061405429d0/
       sig=sign(1.0d0,x)
       xabs=dabs(x)
       if (xabs.gt.20) then
c     
c     *** vermeidung von underflow ***
c
        erf=sig
       else
        t=1.0d0/(1.0d0+p*xabs)
        erf=1.0d0-((((a5*t+a4)*t+a3)*t+a2)*t+a1)*t*exp(-xabs*xabs)
        erf=erf*sig
       endif
       return
       end

