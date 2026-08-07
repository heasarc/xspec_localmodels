!------------------------------------------------------------------------
      program test
!
!  The main program and the geobins subroutine are just wrapper programs
!  which take the place of XSPEC for testing purposes. XSPEC directly
!  calls propfluc()
!
      implicit none
      integer nf,ifl,n,i
      real param(21),far(0:4000),photar(4000),photer(4000)
      real dt,c,f,dfbin
      param(1)  = 13.91        !Sigma0
      param(2)  = 7.36         !rbw
      param(3)  = 3.0          !kappa
      param(4)  = 0.9          !lambda
      param(5)  = 0.0          !zeta
      param(6)  = 0.165        !Fvar
      param(7)  = 25.63        !ro
      param(8)  = 3.3          !ri
      param(9)  = 10.0         !Q
      param(10) = 10.0         !Qsub (sub-harmonic)
      param(11) = 0.1638       !n_qpo
      param(12) = 0.049        !n_2qpo
      param(13) = 0.0249       !n_3qpo
      param(14) = 0.0          !n_0.5qpo
      param(15) = 3.738        !gamma_h: hard band emissitivy index
      param(16) = 3.0          !gamma_s: soft band emissivity index
      param(17) = 10.0         !M
      param(18) = 0.5          !a
      param(19) = 35.0         !Ndec
      param(20) = 1.0          !mode: 1 = hard band PSD, 2 = soft band PSD, 3 = lag spectrum
      param(21) = 1.0          !conv: 0 = add QPO, 1 = multiply QPO
      !-----------------------------------------------------------------
      n  = 2**15
      dt = 1.0/256.0
      c  = 1.05
      call geobins(dt,n,c,far,nf)
      call propfluc(far,nf,param,ifl,photar,photer)
      !Write out in fort.99
      do i = 1,nf
        f     = ( far(i) + far(i-1) ) / 2.0
        dfbin = far(i) - far(i-1)
        write(99,*)f,photar(i)/dfbin
      end do
      write(99,*)"log"
      write(99,*)"li s on"
      end
!------------------------------------------------------------------------

!------------------------------------------------------------------------
      subroutine geobins(dt,n,c,far,nf)
      !A subroutine to define gemetric frequency bins
      implicit none
      integer n,i,nf,loop,nptot
      real dt,c,far(0:400),points(400)
      real df,np
      loop=1
      df=1./(real(n)*dt)
      np=1.0
      nptot=0
      far(0)=df/2.0
      i=1
      do while(loop.eq.1)
          np=real(c*np)
          nptot=nptot+int(np)
          if(nptot.gt.n/2)then
             np=np-real(nptot)+n/2
             nptot=n/2
             loop=0
          end if
          points(i)=float( int(np) )
          far(i)=(real(nptot)+0.5)*df
          i=i+1
      end do
      nf=i-1
      return
      end
!------------------------------------------------------------------------



!------------------------------------------------------------------------
      subroutine propfluc(far,nf,param,ifl,photar,photer)
! PROPFLUC (Ingram & Done 2011; 2012; Ingram & van der Klis 2013)
! Parameters:
! 1  - Sigma0 (normalisation of surface density profile)
! 2  - rbw (bending wave radius)
! 3  - kappa
! 4  - lambda
! 5  - zeta (these three parameters govern surface density profile)
! 6  - Fvar (Fractional variability per decade in radius)
! 7  - ro (truncation radius)
! 8  - ri (inner radius of flow)
! 9  - Q (quality factor, Q=centroid/FWHM, of the QPO harmonics)
! 10 - Qsub (quality factor of the subharmonic can be different)
! 11 - n_qpo (Normalisation of fundamental)
! 12 - n_2qpo (Normalisation of 2nd harmonic)
! 13 - n_3qpo (Normalisation of 3rd harmonic)
! 14 - n_05qpo (Normalisation of sub-harmonic)
! 15 - gamma_h (hard band emissivity index)
! 16 - gamma_s (soft band emissivity index)
! 17 - M (BH mass in solar masses)
! 18 - a (dimensionless spin parameter)
! The following are just settings
! 19 - Ndec (number of rings per decade in radius)
! 20 - mode (1 = hard band PSD, 2 = soft band PSD, 3 = lag spectrum)
! 21 - conv (0 = add QPO, 1 = multiply QPO)
      implicit none
      integer nf,n,i,ifl,nmax
      parameter (nmax=2**20)
      real fmax,fmin,far(0:nf),param(21),photar(nf),photer(nf)
      real dt,P(nmax/2),df,f,Pin,dfbin
      logical firstcall
      save firstcall
      data firstcall/.true./
      !Dummy variable for xspec
      ifl  = 1
      !Header
      if(firstcall)then
        write(*,*) '--------------------------------------------------'
        write(*,*) 'This is a power/lag spectrum model by Adam Ingram.'
        write(*,*) 'Details of physical assumptions are in:'
        write(*,*) 'Ingram & Done (2012, MNRAS, 419, pg 2369).'
        write(*,*) 'Details of analytic formalism are in:'
        write(*,*) 'Ingram & van der Klis (2013).'
        write(*,*) 'Please cite these papers if you use this model.'
        write(*,*) '--------------------------------------------------'
        firstcall=.false.
      end if
      !Set n and dt
      fmax = far(nf)
      fmin = far(0)
      dt   = 1.0 / ( 2.0 * fmax )
      n    = 2**( ceiling( -log(fmin*dt) / log(2.) ) )
      if(n.gt.2**15) n = 2**15
      !Call the model with a working grid
      call dopropfluc(param,n,nmax,dt,p)
      !Now bin/interpolate onto the output grid
      df = 1.0 / ( float(n)*dt )
      do i = 1,nf
        f     = ( far(i) + far(i-1) ) / 2.0
        dfbin = far(i) - far(i-1)
        call modbin(dt,n,P,f,dfbin,Pin)
        photar(i) = Pin * dfbin
        photer(i) = 0.0
      end do
      return
      end
!------------------------------------------------------------------------


!=======================================================================
      subroutine dopropfluc(param,n,nmax,dt,pout)
!
! This is the main routine that does the calculation
!
      implicit none
      integer j,i,n,k,m,nmax,mmax,mode
      parameter (mmax=400)
      real p(0:nmax/2),mu,sig,dt,x,df,cp(0:nmax/2),prev(0:nmax/2)
      integer l,il(mmax),conv
      real h(0:mmax),pi,ph(0:nmax/2),dph,ps(0:nmax/2),dps
      real mf(0:mmax,0:nmax/2),phase,muk(mmax)
      real pk(mmax,0:nmax/2),param(*)
      real hard,soft,dronr,fv,a,s(0:mmax)
      real Td,dTd,Jd,dJd,fqpo,ptot(0:nmax/2)
      real pout(nmax/2),rc(0:nmax/2),ic(0:nmax/2),drc,dic,tlag(0:nmax/2)
      !Define pi
      pi = acos(-1.0)
      !Define number of rings and radial grid
      call rgrid(param,m,a,dronr)
      !Define frequency grid
      df  = 1.0 / (float(n)*dt)
      !Define mean and standard deviation of fluctuations
      mu  = 1.0
      sig = param(6) / sqrt( param(19) )
      !Initialize
      do j = 0,n/2
          mf(0,j) = 0.0
          ph(j)   = 0.0
          ps(j)   = 0.0
          rc(j)   = 0.0
          ic(j)   = 0.0
      end do
      h(0) = 0.0
      s(0) = 0.0
      Td   = 0.0
      Jd   = 0.0
      !Loop through rings
      do k = 1,m
          !Define viscous frequency plus hard and soft emissivities
          call ring(k,param,a,dronr,fv,hard,soft,dTd,dJd)
          Td    = Td + dTd
          Jd    = Jd + dJd
          h(k)  = hard
          s(k)  = soft
          il(k) = int( dronr / (fv*dt) )
          !Define input power spectrum
          call mkpwr(n,dt,mu,sig,fv,0.0,p)
          !My FFT convolution
          if(k.eq.1) then
             do j = 0,n/2
                 cp(j) = p(j)
             end do
          else
             do j = 0,n/2
                 prev(j) = cp(j)
             end do
             call fftconv(dt,n,p,prev,cp)
          end if
          !Record useful vectors and do simple sum
          !See equations 20, 21 & 22 from Ingram & van der Klis (2013)
          muk(k) = sqrt( p(0) )
          do j = 0,n/2
              pk(k,j) = p(j)
              mf(k,j) = cp(j)
              ph(j)   = ph(j) + h(k)**2.*mf(k,j)
              ps(j)   = ps(j) + s(k)**2.*mf(k,j)
              rc(j)   = rc(j) + h(k)*s(k)*mf(k,j)
          end do
          !Cross terms
          phase = 0.0
          do i = 1,k-1
              l = k - i
              phase = phase + 2.*pi*float(il(l+1))/float(n)
              do j = 0,n/2
                  x     = mod(phase*float(j),2.0*pi)
                  dph   = 2.*h(l)*h(k)*cos(x)*mf(l,j)
                  ph(j) = ph(j) + dph
                  dps   = 2.*s(l)*s(k)*cos(x)*mf(l,j)
                  ps(j) = ps(j) + dps
                  drc   = ( h(k)*s(l) + h(l)*s(k) )
                  drc   = drc * cos(x)*mf(l,j)
                  rc(j) = rc(j) + drc
                  dic   = ( h(k)*s(l) - h(l)*s(k) )
                  dic   = dic * sin(x)*mf(l,j)
                  ic(j) = ic(j) + dic
              end do
          end do
      end do
      !Include the QPO
      fqpo = Td / Jd
      mode = int( param(20) )
      conv = int( param(21) )
      if(mode.eq.1)then
         call incQPO(dt,param,conv,fqpo,n,ph,ptot)
      else if(mode.eq.2)then
         call incQPO(dt,param,conv,fqpo,n,ps,ptot)
      else if(mode.eq.3)then
         call lincQPO(dt,param,conv,fqpo,n,rc,ic,tlag)
      end if
      !Write output
      do j = 1,n/2
          if(mode.eq.3)then
             pout(j) = tlag(j)
          else
             pout(j) = ptot(j)
          end if
      end do
      return
      end
!=======================================================================


!============Routines relating to physical assumptions==================

!-----------------------------------------------------------------------
      subroutine ring(k,par,a,dronr,fv,h,s,dTd,dJd)
! Calculates the viscous frequency and LT frequency at a given ring
      implicit none
      integer k
      real par(30),Rg,fv,r,a
      real rbw,dronr,ro
      real x,pi,c,Sigma0,gammaH,gammaS,h,s
      real Sigma,kappa,zeta,lambda
      real fk,flt,dTd,dJd
      !Definitions
      pi = 6.0*asin(0.5)
      c  = 3e8
      !Read in neccessary parameters
      Sigma0 = par(1)
      rbw    = par(2)
      kappa  = par(3)
      lambda = par(4)
      zeta   = par(5)
      ro     = par(7)
      gammaH = par(15)
      gammaS = par(16)
      Rg     = par(17) * 1474.8
      !Calculate the radial grid
      r        = ro * ( a**float(k-1) + a**float(k) ) / 2.0
      !Calculate surface density
      x  = r / rbw
      Sigma = x**lambda/( (1.0+x**kappa)**((zeta+lambda)/kappa) )
      Sigma = Sigma * Sigma0
      !Calculate viscous frequency
      fv = (1.+x**kappa)**((lambda+zeta)/kappa)/(x**(lambda+2.))
      fv = fv / (2.*pi*Sigma0*rbw**2.)
      fv = fv * c/Rg
      !Calculate emissivities
      h = dronr * r**(2.0-gammaH) * Sigma
      s = dronr * r**(2.0-gammaS) * Sigma
      !Calculate LT torque and disc angular momentum for the ring  
      fk  = r**(-1.5) / (2.0*pi) * c/Rg
      flt = 1.-sqrt(1.-4.*par(18)*r**(-1.5)+3.*par(18)**2.*r**(-2.))
      flt = flt * fk
      dTd = flt*fk*Sigma*r**4.0*dronr
      dJd =     fk*Sigma*r**4.0*dronr
      return
      end        
!-----------------------------------------------------------------------

!===============QPO subroutines=========================================


!-----------------------------------------------------------------------
      subroutine incQPO(dt,param,conv,fqpo,n,psum,ptot)
      implicit none
      integer n,nmax,conv,j
      parameter (nmax=2**20)
      real dt,param(*),fqpo,psum(0:nmax/2),ptot(0:nmax/2)
      real pqpo(0:nmax/2),mu
      do j = 1,n/2
          psum(j) = psum(j) / psum(0)
      end do
      psum(0) = 1.0
      if(conv.eq.0)then
         mu = 0.0
         call QPOpower(dt,param,mu,fqpo,pqpo,n)
         do j = 0,n/2
             ptot(j) = psum(j) + pqpo(j)
         end do
      else
         mu = 1.0
         call QPOpower(dt,param,mu,fqpo,pqpo,n)
         call fftconv(dt,n,psum,pqpo,ptot)
      end if
      return
      end
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine lincQPO(dt,par,conv,fqpo,n,rc,ic,tlag)
      implicit none
      integer n,nmax,conv,j
      parameter (nmax=2**20)
      real dt,par(30),fqpo,rc(0:nmax/2),ic(0:nmax/2)
      real pqpo(0:nmax/2),tlag(0:nmax/2),cov0,mu,f,pi
      real rcqpo(0:nmax/2),icqpo(0:nmax/2),df
      pi = acos(-1.0)
      df  = 1.0 / (float(n)*dt)
      cov0 = rc(0)**2.0 + ic(0)**2.0
      do j = 0,n/2
          rc(j) = rc(j) / sqrt( cov0 )
          ic(j) = ic(j) / sqrt( cov0 )
      end do
      if(conv.eq.0)then
         mu = 0.0
         call QPOpower(dt,par,mu,fqpo,pqpo,n)
         do j = 1,n/2
             rc(j) = rc(j) + pqpo(j)
             f = float(j) * df
             tlag(j) = atan( ic(j) / rc(j) )
             if(ic(j).gt.0..and.rc(j).lt.0.) tlag(j) = tlag(j) + pi
             if(ic(j).lt.0..and.rc(j).lt.0.) tlag(j) = tlag(j) - pi
             tlag(j) = tlag(j) / ( 2.0 * pi * f )
         end do
      else
         mu = 1.0
         call QPOpower(dt,par,mu,fqpo,pqpo,n)
         call fftconv(dt,n,rc,pqpo,rcqpo)
         call fftconv(dt,n,ic,pqpo,icqpo)
         do j = 1,n/2
             rc(j) = rcqpo(j)
             ic(j) = icqpo(j)
             f = float(j) * df
             tlag(j) = atan( ic(j) / rc(j) )
             if(ic(j).gt.0..and.rc(j).lt.0.) tlag(j) = tlag(j) + pi
             if(ic(j).lt.0..and.rc(j).lt.0.) tlag(j) = tlag(j) - pi
             tlag(j) = tlag(j) / ( 2.0 * pi * f )
         end do
      end if
      return
      end
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine QPOpower(dt,param,mu,fqpo,Pqpo,n)
      implicit none
      integer n,j,nmax,k
      parameter (nmax=2**20)
      real dt,param(*),mu,fqpo,Pqpo(0:nmax/2)
      real Q,fac,width,fcen,Norm,P(0:nmax/2)
      !Initialize
      Pqpo    = 0.0
      Pqpo(0) = mu**2.
      Q       = param(9)
      !Loop through harmonics
      do k = 1,4
        fac   = float(k)
        if( k .eq. 4 )then
          fac = 0.5
          Q   = param(10)
        end if
        fcen  = fac * fqpo
        width = fcen / ( 2. * Q )
        Norm  = param(10+k)
        call mkpwr(n,dt,mu,Norm,width,fcen,P)
        do j = 1,n/2
          Pqpo(j) = Pqpo(j) + P(j)
        end do
      end do
      return
      end
!-----------------------------------------------------------------------


!==================Simple mathematical routines=========================

!-----------------------------------------------------------------------
      subroutine mkpwr(np,dt,mu,r,Delta,f0,ap)
      implicit none
      integer j,np
      real df,dt,f,ap(0:np/2),r,Delta,pi,f0,mu,norm
      pi = acos(-1.0)
      df = 1.0 / ( float(np)*dt )
      ap(0) = mu**2.0
      norm = r**2.0 / ( pi/2.0 + atan(f0/Delta) )
      do j = 1, np/2
          f     = float(j) * df
          ap(j) = Delta / ( Delta**2. + (f-f0)**2. )
          ap(j) = ap(j) * norm
      end do
      return
      end
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine rgrid(par,m,a,dronr)
      implicit none
      real par(30),ro,ri,Ndec,a,dronr
      integer m
      ro     = par(7)
      ri     = par(8)
      Ndec   = par(19)
      m      = int(Ndec*(log10(ro)-log10(ri)))
      a      = exp( log(ri/ro) / float(m) )
      dronr  = 2.0 * (1.0-a) / (1.0+a)
      return
      end
!-----------------------------------------------------------------------


!====================binning routines===================================

!-----------------------------------------------------------------------
      subroutine modbin(dt,n,P,fbin,dfbin,Pbin)
      !Binning P into already defined frequency bins
      !If the bin is going to be empty, it interpolates.
      !Is designed to be called for every point.
      implicit none
      integer n,nmax,j,jmin,jmax
      parameter (nmax=2**20)
      real dt,fbin,dfbin,Pbin,P(nmax/2),df,T
      df   = 1.0 / (float(n)*dt)
      T    = float(n) * dt
      jmax = floor( (fbin + dfbin/2.) / df )
      jmin = ceiling( (fbin - dfbin/2.) / df )
      Pbin = 0.0
      if(jmax.gt.jmin)then
         !Bin
         do j = jmin,jmax
             Pbin = Pbin + P(j)
         end do
         Pbin = Pbin / float(jmax-jmin+1)
      else
         !Interpolate
         call interpolate(fbin,T,df,P,nmax,n,Pbin)
      end if
      return
      end
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine interpolate(fin,T,df,P,npmax,np,Pin)
      implicit none
      integer j,js,np,npmax
      real Pin,P(npmax/2),fr,m,c,fin,T,df
      j   = ceiling( fin*T )
      js  = floor( fin*T )
      if(js.lt.1)then
         j  = 2
         js = 1
      else if(j.gt.np/2)then
         j  = np/2
         js = np/2 - 1
      end if
      if(j.eq.js)then
         Pin = P(j)
      else
         m   = ( P(j)-P(js) ) / df
         fr  = float(j) / float(js)
         c   = ( fr*P(js) - P(j) ) / ( fr - 1.0 )
         Pin = m*fin + c
      end if
      return
      end
!-----------------------------------------------------------------------


!=======Fourier routines that do the analytic calculation===============

!-----------------------------------------------------------------------
      subroutine fftconv(dt,n,ap,bp,cp)
      
      ! to call the fourier transform subroutine
      USE FFTW3_wrap
      
      !Does a convolution using FFTs
      implicit none
      integer n,j,k,i
      real ap(0:n/2),bp(0:n/2),cp(0:n/2)
      real adata(2*n),bdata(2*n),cdata(2*n),dt
      !Fill data vectors
      adata(1) = ap(0) * float(n) * dt
      adata(2) = 0.0
      bdata(1) = bp(0) * float(n) * dt
      bdata(2) = 0.0
      do j = 1, n/2
          k          = 2*j + 1
          adata(k)   = ap(j)
          adata(k+1) = 0.0
          bdata(k)   = bp(j)
          bdata(k+1) = 0.0
          i          = 2*n + 2 - k
          adata(i)   = adata(k)
          adata(i+1) = adata(k+1)
          bdata(i)   = bdata(k)
          bdata(i+1) = bdata(k+1)
      end do
      !Now do the inverse FFT
      call fftw_wrap_four1_sng(adata,n,-1)
      call fftw_wrap_four1_sng(bdata,n,-1)
      !Now multiply together
      do j=1,n
          k = 2*j
          cdata(k)   = adata(k)*bdata(k)
          cdata(k-1) = adata(k-1)*bdata(k-1)
      end do
      !Then transform back
      call fftw_wrap_four1_sng(cdata,n,1)
      !Finally put the +ve real frequencies into cp(n/2)
      do j = 0,n/2
          k = 2*j + 1
          cp(j) = cdata(k) / (float(n**2)*dt)
      end do
      cp(0) = cp(0) / (float(n) * dt)
      return
      end
!-----------------------------------------------------------------------

