      Subroutine sssed(ear,ne,param,ifl,photar,photer)
      
      use xsfortran
      
      implicit none

      integer NPAR
      parameter(NPAR=14)

      integer ne, ifl
      real ear(0:ne), photar(ne), param(NPAR), photer(ne)

!------------------------------------------------------------------------------------ 
! This subroutine checks whether any of the parameters have changed, if they have, 
! then it creates a new energy binning scheme with 5000 bins increasing exponentially
! in size.  The input parameters and new sbinning scheme are used in a call to the 
! subroutine 'asdiskf'.  The returned fluxes are then corrected for redshift and
! shifted back into the original energy binning scheme.
!------------------------------------------------------------------------------------ 

!     ear = Energy array
!     ne = size of energy array (1 more than flux array as gives boundaries)
!     param = array containing model parameters
!     ifl = spectrum number
!     photar = output flux array (ne-1 parameters)

      integer NN
      parameter(NN=5000)

      integer n,i
      real e(0:NN),ph(NN)
      real dloge,zfac,oldpar(NPAR),newemin,newemax
      logical parchange, echange

      REAL fstart(ne), fend(ne)
      integer istart(ne), iend(ne)

      save oldpar,ph,e

c     param(1) mass in solar
c     param(2) luminosity distance
c     param(3) log mass accretion rate in L/LEdd
c     param(4) rin
c     param(5) cosi 0.5 is isotripic  (14)
c     param(6) Te of hard compton (5)
c     param(7) Te of warm compton
c     param(8) gamma hot       (9)
c     param(9) gamma warm      (8)
c     param(10) Rhot
c     param(11) Rwarm (15)
c     param(12) log10 rout/rg (6)
c     param(13) redshift (11)
c     param(14) colour correction on/off 0=iff 1=on

      DO i = 1, ne
         photer(i) = 0.0
      ENDDO

      zfac = 1.0 + param(13)

!     test whether any parameters have changed
      parchange=.false.
      do i=1,NPAR,1
         if (param(i).ne.oldpar(i)) parchange=.true.
      end do
!     test whether own energy grid needs to be changed
      if (ear(0).eq.0.0) then 
         newemin=ear(1)-(ear(2)-ear(1))/10.0
      else
         newemin=min(1.0e-5,ear(0)*zfac)
      end if
      newemax=max(1.0e3,ear(ne)*zfac)

      if ((e(0).ne.newemin).or.(e(NN).ne.newemax)) then
         echange=.true.
      end if

!     if necessary recalculate the energy grid

      if (echange) then

         dloge=log10(newemax/newemin)/float(NN)
         e(0) = newemin
         do n=1,NN,1
!           populate the new energy array and zero flux
            e(n)=10**(log10(e(0))+dloge*float(n))
         end do

      endif
!     call the irradiated disc if parameters or energy grid have changed
      if (parchange .or. echange) then 

         call ssdiskf(e,NN,param,ifl,ph)

!        now redshift the energy bins back
         do n=0,NN,1
            e(n)=e(n)/zfac
         end do

!        and now correct the flux for redshift
         do n=1,NN,1
            ph(n)=ph(n)/zfac
         end do

!        save the new parameters into the oldpar array
         do i=1,NPAR,1
            oldpar(i)=param(i)
         end do

      end if

!     rebin the calculated fluxes back onto original energy grid
      CALL ainibin(NN, e, ne, ear, istart, iend, fstart, fend, 0)
      CALL aerebin(NN, ph, ne, istart, iend, fstart, fend, photar)
     
      return
      end


      subroutine ssdiskf(ear,ne,param,ifl,photar)

      use xsfortran

c     program to integrate the disk equations from shakura-sunyaev disk
c     as given by Novikov and Thorne

      implicit none
      double precision tseedpl
      double precision m,mdot,rsg
      double precision rep
      double precision rin
      double precision pow,rcor,sstemp
      double precision seeddis
      double precision frac
      double precision tnt
      double precision rgcm,pi,t,t0,r,dr,dlogr,dflux,tprev
      double precision dfluxint,dfluxseed !add add add
      double precision en,kkev,h,kevhz,d0,d,logrout
      double precision dllth,dldiskint !add add add
      integer i,icor,imax,iout,n,ne,ifl
      logical first,firstdisk,firstdisk2,firstdisk3,firstpl,firstcomp
      double precision flux(ne),displ,discomp,discompall
      double precision mdotedd,alpha,eff
      double precision fluxint(ne),fluxseed(ne) ! add add add
      double precision cosi
      double precision gammah,gammas,teh,tes
      double precision gparamh, gparams,gfach, gfacs
      double precision ledd
      double precision kappa
      double precision fcol,tprevcol,sw
      double precision tcol,tsgcol
      double precision tcor,tpow
      double precision rmax,tsg !critical temperature[K] for soft compton
      real ear(0:ne),photar(ne),param(*)
      real lpar(5),lphot(ne),lphote(ne),hpar(5),hphot(ne),hphote(ne)
      real lphotall(ne),hphotall(ne) !add add add
      character(255) comment

c     constants
      h=6.62617d-27 ! [erg s]
      kkev=1.16048d7  ![K/keV]
      kevhz=2.417965d17 ![Hz/keV]
      pi=4.0*atan(1.0)      
      sw=dble(param(14))
c     system parameters
      tseedpl=-1.0d0
      m=dble(param(1))         !in solar units
c      Htmax=dble(param(13))  !upper limit of scaleheight
      d0=dble(param(2))     !in kpc
      cosi=dble(param(5))
      mdotedd=dble(10**(param(3)))      !in L/Ledd if -ve plot disc
      rin=dble(param(4))     
      rmax=rin*(7.0/6.0)**2
      alpha=0.1
      rep=0.0
c     corona parameters
      gammah=dble(abs(param(8)))
      gammas=dble(abs(param(9)))  

      teh=dble(abs(param(6)))/511.0
      tes=dble(abs(param(7)))/511.0
      rcor=dble(abs(param(11)))
c      gparamh=dble(param(16)) !0~1 0:SLAB, 1:ISOTROPIC
c     gparams=dble(param(17)) !0~1 0:SLAB, 1:ISOTROPIC
      gparamh = 1
      gparams = 1
      gfach=dble((cosi+(0.5-cosi)*gparamh)/0.5)
      gfacs=dble((cosi+(0.5-cosi)*gparams)/0.5)      
c     d=d0*1d6*3.085677d18 !mod
      d=d0*1d3*3.085677d18 !mod      
      rgcm=1.477d5*m  ![cm]
      eff=1/(rin*2.0)
      mdot=mdotedd*1.39d18*m/8.98755 ![g/s]      
      ledd=1.39d38*m

      logrout=dble(param(12))   !log10 Rout/Rg
c     get reprocess parameter
c      a=0.30d0!albedo
c      hrdc = 0.3d0
c     iv -ve then calculate rout=rsg frmo laor & netzer 1989
      if (param(12).lt.0.0) then
         logrout=((m/1e9)**(-2.0/9.0)) * ((mdotedd)**(4.0/9.0))
         logrout=2150*logrout* (alpha**(2.0/9.0))
         logrout=log10(logrout)
      end if
      rsg=10**logrout !rout
      tsg=sstemp(m,mdot,rin,rsg)
      tsgcol=tsg*fcol(tsg,sw)
      frac=dble(param(10))

c     initialise 
      do n=1,ne,1
         photar(n)=0.0
         flux(n)=0.0
         fluxint(n)=0.0
         fluxseed(n)=0.0
         lphotall(n)=0.0 !add add add 
         hphotall(n)=0.0 !add add add 
      end do
      
      imax=2000
      iout=1000
      icor=200
c      icor=500
      first=.true.
      firstdisk=.true.
      firstdisk2=.true.
      firstdisk3=.true.
      firstpl=.true.
      firstcomp=.true.      
      displ=0.0d0
      discomp=0.0d0
      discompall=0.0d0            
      seeddis=0.0d0
      t0=0.0d0
      tpow=0.0d0
      tcor=0.0d0
      tprev=0.0d0
      dflux=0.0d0
      dfluxint=0.0d0
      dfluxseed=0.0d0

c calculate Redd@L>Ledd and Rpow@L<Ledd first run
      do i=1,imax,1  !outer sc + disk is calculated
         dlogr=log10(rsg/rin)/float(imax)
         r=10**(log10(rin)+float(i-1)*dlogr+dlogr/2.0)
         dr=10**(log10(r)+dlogr/2.0) - 10**(log10(r)-dlogr/2.0)
         t=sstemp(m,mdot,rin,r) !intrinsic effective temperature [K]
         tcol=fcol(t,sw)*t
c---------
         if(rcor.gt.rin)then
            if (firstcomp) then
               discompall=discompall+2*2*pi*r*dr*rgcm*rgcm*5.670367e-5
     &              *t**4
c               displ=displ+2*2*pi*r*dr*rgcm*rgcm*5.670367e-5*t**4
               if (r.ge.rcor)then
                  firstcomp=.false.
               endif
            endif
         endif
        tnt=sstemp(m,mdot,rin,r) !intrinsic effective temperature [K]
        t=tnt
        kappa=fcol(t,sw)
        tcol=kappa*t
c     setting Rcor for soft compton
        tprev=t
        tprevcol=tcol
            
      enddo                     !end of r

c mod mod mod

c-------set rcor------
      if (rcor.le.rin)then
         icor=0
         rcor=rin
        write(*,*)'============================================'
        write(*,*)'R_outer is smaller than R_inner <reset>'
        write(*,*)'============================================'
      endif
c--------------------
c------- set rcor------0109_rev
      if (rcor.gt.rsg)then
         rcor=rsg
        write(*,*)'============================================'
        write(*,*)'R_outer is larger than Rout <reset>'
        write(*,*)'============================================'
      endif
c--------------------
ccccccccccccccccccccccccccccccccc
      do i=1,icor+iout,1
         if (i.le.icor) then
            dlogr=log10(rcor/rin)/float(icor)
            r=10**(log10(rin)+float(i-1)*dlogr+dlogr/2.0)
         else
            dlogr=log10(rsg/rcor)/float(iout)
            r=10**(log10(rcor)+float(i-icor-1)*dlogr+dlogr/2.0)
         end if
c     geometry factor
c         theta0=asin(Ht/r)
         dr=10**(log10(r)+dlogr/2.0) - 10**(log10(r)-dlogr/2.0)
         tnt=sstemp(m,mdot,rin,r) !intrinsic effective temperature [K]
         t=tnt
         kappa=fcol(t,sw)
         tcol=kappa*t
         t=t/kkev               ! [keV]
         tcol=tcol/kkev

c     go over each photon energy - midpoint of bin
         do n=1,ne,1
            en=dble(log10(ear(n))+log10(ear(n-1)))
            en=en/2.0
            en=10**en
           
c     do blackbody spectrum  @r>rcor
            if ((en.lt.30.0*tcol).and.(r.gt.rcor)) then
               dflux=pi*2.0*h*((en*kevhz)**3)/8.98755d20
               dflux=dflux*4.0*pi*r*dr*rgcm*rgcm
               dflux=dflux/(exp(en/(tcol))-1.0)
               dflux=dflux/(kappa**4)
            else
               dflux=0.0d0
            end if

            flux(n)=flux(n)+(dflux)
            
         end do 


c  start ---- mod mod mod mod slice nthcomp
         dllth=0.0
         pow=0.0
         dldiskint=2.0*2.0*pi*r*dr*rgcm*rgcm*5.670367e-5
     &        *((t*kkev)**4) !mod mod mod
         dldiskint=dldiskint/(4.0*pi*d*d) ! erg/s/cm^2

ccccccccccccccc

         if (i.le.icor) then
            hpar(1)=sngl(gammah)
            hpar(2)=abs((param(6))) !Te
            hpar(3)=SNGL(tcol) !keV
            hpar(4)=0.0
            hpar(5)=0.0
            lpar(1)=abs(sngl(gammas))
            lpar(2)=abs(param(7)) !kTe if - plot comp
            lpar(3)=sngl(tcol) !set to reprocessed effective temperature mod mod
            lpar(4)=0.0         ! int type
            lpar(5)=0.0
            
            ifl=1
            call donthcomp(ear,ne,hpar,ifl,hphot,hphote)
            pow=0.0d0
            do n=1,ne,1
               pow=pow+hphot(n)*ear(n)*kevhz*h !hard compton luminosity
            end do

               ifl=1
c thcomp is called but scale may be wrong yet               
               call donthcomp(ear,ne,lpar,ifl,lphot,lphote)
               do n=1,ne,1
                  dllth=dllth+lphot(n)*ear(n)*kevhz*h !photons/s/cm2/Hz *Hz *keV
               end do
         end if

         do n=1,ne,1
            if (pow.eq.0) then
               hphot(n)=0.0
            else 
               hphot(n)=hphot(n)*sngl(dldiskint*frac/pow)
            endif 
            hphotall(n)=hphotall(n)+hphot(n)               
         end do 
cccccccccccccc
         
         do n=1,ne,1
            if (dllth.eq.0) then
               lphot(n)=0.0
            else 
               lphot(n)=lphot(n)*sngl(dldiskint*(1-frac)/dllth)
            endif 
            lphotall(n)=lphotall(n)+lphot(n)               
         end do 

c  end ----- mod mod mod mod slice nthcomp
         
      end do                    !end of r

      do n=1,ne,1
c        this is ergs cm^2 s-1 Hz^-1 
         flux(n)=flux(n)/(4.0*pi*d*d)
         fluxseed(n)=fluxseed(n)/(4.0*pi*d*d) 
c        photons is flux/hv - photons cm^2 s-1 Hz^-1
         flux(n)=flux(n)/(h*kevhz*ear(n))
         fluxseed(n)=fluxseed(n)/(h*kevhz*ear(n))!intrinsic disk flux at r>rpow

c        now multiply by energy band in Hz
         photar(n) = sngl(flux(n)*kevhz) * (ear(n)-ear(n-1)) !kept 

      end do

      write(comment,*)'-----------'
      CALL xwrite(comment, 20)
      write(comment,*)'T(Rcor)=',
     & sstemp(m,mdot,rin,rcor)/kkev,'keV'
      CALL xwrite(comment, 20)
      write(comment,*)'Rout=', rsg,'Rg'
      CALL xwrite(comment, 20)
      write(comment,*) 'Tout=',tsg/kkev,'keV'
      CALL xwrite(comment, 20)
      write(comment,*)'-----------'



      do n=1,ne,1
         if ((param(6).lt.0.0).or.(param(7).lt.0.0).or.
     &        (param(9).lt.0.0)) then
            if (param(6).lt.0.0)  photar(n)=
     &           SNGL(hphotall(n)*gfach) !inner corona    
            if (param(7).lt.0.0)  photar(n)=   
     &           SNGL(lphotall(n)*gfacs) !disc-corona
            if (param(9).lt.0.0) photar(n)=photar(n)
     &           *sngl(cosi/0.5) ! outer disc
         else
            photar(n)=SNGL(photar(n)*cosi/0.5
     &           +lphotall(n)*sngl(gfacs)
     &           +hphotall(n)*sngl(gfach))
         end if


      end do

      return
      end


      function sstemp(m0,mdot0,rin,r0)
c     find temperature as a function of mass, spin and r in units of rg

        implicit none
        double precision m0,mdot0,r,r0,rgcm,y,yms
        double precision pi,rin,c,b
        double precision sstemp

        pi=4.0*atan(1.0)
        rgcm=1.477d5*m0

        r=r0
        y=sqrt(r0)
        yms=sqrt(rin)
        c=1.0-yms/y
        b=1
        sstemp=3.0*6.67e-8*1.989d33*m0*mdot0
        sstemp=sstemp/(8.0*pi*5.670367e-5*(r*rgcm)**3)
        sstemp=(sstemp*c/b)**0.25
        return
        end


      function fcol(t,sw)
c     find temperature as a function of mass, spin and r in units of rg

        implicit none
        double precision t,kkev,fcol,sw
        kkev=1.16048d7          ![K/keV]

        if (sw.ne.0.0)then
           if (t.gt.1.0e5) then
              fcol=(72.0/(t/kkev))
              fcol=fcol**(1.0/9.0)
           else if (t.gt.3.0e4) then
c     linear from 1 to 2.7 between 4e4 and 1e5
              fcol=(t/3.0e4)**0.82
           else
              fcol=1.0
           end if
        else
           fcol=1.0
        endif
        return
        end





