c Local copy of installed TRANSOPTR RFQ.f; routines renamed.
c Fix: use COMMON/MOM/ CURRENT (bunch charge in mode 5).
c Original SHA256: 233231dd7604ce4133e3df6fc07545d017c0622e8b30c7338cb67c49d9bac8e1
      SUBROUTINE RFQSC(IUNIT,IPO,VSCALE,ELu,freq,phasedeg,nsubdr)
      CHARACTER*16 QUA

c    	          Olivier Shelbaya, TRIUMF
c TRANSOPTR subroutine for an RFQSC, described by a potential 
c expansion. Currently supports 2-term potential. Subroutine expects
c an input file containing:
c
c 	           z      A01	A10         k
c
c Where z is measured from the start of the vanes, A01 is the transverse
c focusing strength parameter and A10 the longitudinal. Units of k are 
c inverse units of z, scaled w/ /CONS/conx(8). While in the field, momenta 
c are referred to the INITIAL referencec momenta, so momenta are no 
c longer angles. This simplifies calc. At the end, momenta are suddenly 
c changed to reference them to the new central momentum.
c
c     Only works in space charge mode.
c     ARGUMENTS...
c     IUNIT         unit # from which axial potential values are read
c     IPO          =# of z,A01,A10,k lines to be read.
c     VSCALE        vane voltage [MV]
c     ELu           The IPO points cover this real length
c     freq          frequency (omega/2/pi) in Hz
c     phasdeg       phase in degrees
c     nsubdr        = not used any more (the number of subdivision is now controlled by CMPS)


      
      QUA="RFQSC"
      CALL RFQSCn(IUNIT,IPO,VSCALE,ELu,freq,phasedeg,nsubdr,qua)
      RETURN
      END
      
      
      
      SUBROUTINE RFQSCn(IUNIT,IPO,VSCALE,ELu,freq,phasedeg,nsubdr,sol)
      REAL VANEVOLT
      EXTERNAL SCRFQSC
      COMMON/ZED/ZINIT,Z
      COMMON/EMIT/ Emitx,Emity
      COMMON/PRINT/IPRINT,IQ1,JQ1,IQ2,JQ2,IQ3,JQ3,IQ4,JQ4
      COMMON/FMATRX/F(6,6),FT(6,6)
      COMMON/SCPARM/QSC,ISC,CMPS
      COMMON/MOM/P,BRHO,PMASS,ENERGK,GSQ,ENERGKI,CHARGE,CURRENT,ENERGKS
      COMMON/VECTORS/VECI(6),VECTOR(6),RVECT(6),NPARAM,IVOPT
      REAL A(6,6)
      COMMON/rfqvars/HATBL(4,9999),EMTBL(4,9999),CAYTBL(4,9999),VANEVOLT
      COMMON/axezrf/omoc,phase,etot,
     $     ZTBLE(9999),EZTBL(4,9999),zse,iez,iaxezrf,spphase,ezscl,vzscl
      COMMON/CALLS/NSC,nagflag
c
c     etot is energy ((gamma-1)*mass) at start of linac
c     omoc is omega/c = 2pi nu/c
c     zse is z at start of linac
c
      COMMON/CONS/CONX(8),unitu(8)
      COMMON/axez2/pratio,XYASYM
      save
      COMMON/LABELS/start,QNAR,FFIT,skp,Se4,QDOT
      CHARACTER*16 start,QNAR,FFIT,skp,Se4,QDOT,ent,exi
      character (*) :: sol
    
      vzscl=VSCALE
      ent="entr"
      exi="exit"
      IF (IVOPT.EQ.5)GO TO 50
      write(6,*)'Must be in mode 5 (bunched) for this element'
      stop
 50   continue
      iaxezrf=2!option 2 added for RFQSC. 1 is for SCLINAC, 0 for SC.
      EL=ELu/conx(8)

      nsubd=abs(el)/CMPS+1 !number of subdivisions controled by CMPS (cm per step)
      nsubd=((nsubd+1)/2)*2!make it an even number. Important for the XML output to be do right in the middle

      if(iprint.eq.0)nsubd=1 !don't need intermediate pts. when optimizing; only when plotting
      zse=z
      omoc=2.*3.141593*freq/(2.997925e10)
      RFREQ=freq
      phase=phasedeg*3.141593/180.! in radians
      energks = energk ! note kinetic energy at start of accel. element
      etot=energk
      IEZ=ipo

      DO 31 I=1,IEZ

         READ(iunit,*)ZTBLE(I),HATBL(1,I),EMTBL(1,I),CAYTBL(1,I)
         HATBL(2,I)=0.0
         HATBL(3,I)=0.0
         HATBL(4,I)=0.0
         EMTBL(2,I)=0.0
         EMTBL(3,I)=0.0
         EMTBL(4,I)=0.0
         CAYTBL(2,I)=0.0
         CAYTBL(3,I)=0.0
         CAYTBL(4,I)=0.0

 31     CONTINUE
        rewind iunit  
      CLOSE(iunit)

         VANEVOLT=vscale!vane voltage [MV]



C	Here I convert the parameters into an on-axis energy gain.
C	The value EINCREMENT is incremented by negative EINCREMENT
C	which is equal to q(pi)^2A_10 V_0 I0(KR)cos(phi) / 8K, for
C	R = 0 (bessel function ~1), q is the charge and V_0 the vane
C       voltage in MV
C	The computed energy per Z is stored in the array ERFQ(IEZ),
C	which is then interpolated into the array RFQEBTL(IEZ) for
C	use in the subroutine SCRFQSC

c  and spline. 
c first der.=0 at ends
      CAYTBL(2,1)=0.
      CAYTBL(2,IEZ)=0.
      HATBL(2,1)=0.
      HATBL(2,IEZ)=0.
      EMTBL(2,1)=0.
      EMTBL(2,IEZ)=0.
      call cubspl(ZTBLE,CAYTBL,IEZ,1,1)
      call cubspl(ZTBLE,HATBL,IEZ,1,1)
      call cubspl(ZTBLE,EMTBL,IEZ,1,1)


C      Preset THE F-MATRIX. All non-zero elements calculated in SCLINAC
         DO 21 I=1,6
            DO 21 J=1,6
   21          F(I,J)=0.
C         IAXEZ=1           ! z-varying elements are in SC
         if(nsubd.le.0)nsubd=1
         elus=elu/nsubd
         els=el/nsubd
         do i=1,nsubd
      CALL INTEGRATE(SCRFQSC,ELs)
         if(i.eq.1)then
            call velmer(ent,energk,elus,0.,3,
     >-ft(2,1)*pratio,-ft(4,3)*pratio)
         elseif(i.eq.nsubd/2)then
            call velmer(sol,energk,elus,0.,3,
     >-ft(2,1)*pratio,-ft(4,3)*pratio)
         elseif(i.eq.nsubd)then
            call velmer(exi,energk,elus,0.,3,
     >-ft(2,1)*pratio,-ft(4,3)*pratio)
         else
            call velmer(Se4,energk,elus,0.,3,
     >-ft(2,1)*pratio,-ft(4,3)*pratio)
         endif
c velmer occurs inside the loop to get the intermediate results
         IF(iprint.ne.0)write(20,*)z-zse,spphase-phase
         enddo
C
C     Exit matrix converts momenta back into angles
C
         DO 19 I=1,6
         DO 19 J=1,6
         A(I,J)=0.
19       IF(I.EQ.J)A(I,I)=1.
         A(2,2)=PRATIO
         A(4,4)=PRATIO
         A(6,6)=PRATIO

      CALL EDGEUP(A)

c reset space charge parameter for new momentum definition
      QSC=QSC*pratio
      emitx=emitx*pratio
      emity=emity*pratio
      pratio=1.
      energks = energk ! set kinetic energy at end of accel. element
      
      IF(IPRINT.ne.0.and.energk.ne.etot)then

            call velmer(Se4,energk,0.,0.,3
     &,-ft(2,1)*pratio,-ft(4,3)*pratio)

         call energ_ch_out

         endif
      IF(iprint.ne.0)write(22,90)sol,iunit,ezmax,phasedeg,el
      iaxezrf=0
      RETURN

90    FORMAT(' Element: ',A16/9X,
     ,' LINAC: field file unit #',
     $     i4,', max. field = ',
     $     G12.6,' MV/m, Phase = ',
     $     F9.4,' degrees, Length = ',
     $     G12.6,' cm',
     $     /1X,130('*'))

      END


      SUBROUTINE SCRFQSC(Z,SX,DSX)
c
c     Generates DSX from SX. DSX is effectively d/dz of SX.
c      EXTERNAL BESSI0 ! removed, used to call num.rec.
      REAL HA,EM,A01,A10,CAY
      REAL T1QUOT,T1BRACKET,T2QUOT,T2BRACKET,BQUOT,BBRACKET,
     $CQUOT,CBRACKET,cspace,sspace
      COMMON/PRINT/IPRINT,IQ1,JQ1,IQ2,JQ2,IQ3,JQ3,IQ4,JQ4
      COMMON/VELEM/VELM(21,9999),NELM
      COMMON/LABELS/start,QNAR,FFIT,skp,Se4,QDOT
      CHARACTER*16 start,QNAR,FFIT,skp,Se4,QDOT,LABEL
      DIMENSION SX(13,6),DSX(13,6)
      DIMENSION TEMP(6,6),ROT3(3,3),ROT6(6,6),SIG(3,3)
      DIMENSION WRK(3),EVAL(3),RF(6,6)
      COMMON/CALLS/NSC,nagflag
      COMMON/FMATRX/F(6,6),FT(6,6)
      COMMON/SCPARM/QSC,ISC,CMPS
      COMMON/VECTORS/VECI(6),VECTOR(6),RVECT(6),NPARAM,IVOPT
      COMMON/MOM/P,BRHO,PMASS,ENERGK,GSQ,ENERGKI,CHARGE,CURRENT,ENERGKS
      COMMON/axezrf/omoc,phase,etot,
     $     ZTBLE(9999),EZTBL(4,9999),zse,iez,iaxezrf,spphase,ezscl,vzscl
      COMMON/rfqvars/HATBL(4,9999),EMTBL(4,9999),CAYTBL(4,9999),VANEVOLT
      COMMON/axez2/pratio,XYASYM
      COMMON/CONS/CONX(8),unitu(8)

      save
      do i=1,3
      dsx(13,i)=0.
      enddo



      NSC=NSC+1
      pratio=1.
      gamm1=sx(13,6)
      energk=gamm1*PMASS

         if (energk.le.0.)then
            write(6,*)'*********BEAM DECELERATED TO A STOP'
            write(6,*)'Check parameter values or limits.'
            call errout
            endif

      CAY   = ppvalu(ZTBLE,CAYTBL,IEZ,4,Z-zse,0)*conx(8) 
      HA   = ppvalu(ZTBLE,HATBL,IEZ,4,Z-zse,0)
      EM   = ppvalu(ZTBLE,EMTBL,IEZ,4,Z-zse,0)


c OLD/DEPRECATED; assumed input fort file provided (a,m). This has
c been changed so that the user now supplies A01,A10 directly. 
c Allows for avoidance of reliance upon numerical recipes, in 
c particular the bessel function approximation.
c***********Compute Bessel functions approximations*****************
c 2nd order approximation
c      BESSKA = 1. + (CAY*HA/2.)**2
c      BESSMKA = 1. + (EM*CAY*HA/2.)**2
c numerical recipes for fortran 77 below (higher precision)
c      BESSKA = bessi0(CAY*HA)
c      BESSMKA = bessi0(EM*CAY*HA)
c
c*****Compute focussing and acceleraiting factors A01 & A10*********
c      A10 = (EM**2 - 1.)/(BESSKA*EM**2 + BESSMKA)
c      A01 = (1. - A10*BESSKA)/HA**2
c*******************************************************************
      A01 = HA
      A10 = EM

          gamma=1.+gamm1
      ETA=SQRT(gamm1*(gamm1+2.))
      BRHO=ETA*PMASS/CHARGE*3.3356397
      BRHO=BRHO*.001
      P=ETA*PMASS !this P is actually Pc
      gsq=GAMMA**2
      BETA=eta/gamma

          dsx(13,5)=1./beta ! dctime/ds=1./beta
      gamm1old=etot/pmass
      ETAold=SQRT(gamm1old*(gamm1old+2.))
      Pold=ETAold*PMASS
      pratio=etaold/eta

c compute time and spatial sine & cosine terms for matrix operations
      spphase=sx(13,5)*omoc+phase-sx(13,4)
c      write(21,*)z-zse,cay,omoc/beta,cay-omoc/beta
c take me out later. Write statements should not occur in RK integration
           ccc=cos(sx(13,5)*omoc+phase)
           sss=sin(sx(13,5)*omoc+phase)
           cspace=cos(sx(13,4))
           sspace=sin(sx(13,4))
C The argument of what wangler referes to as sin("kz") is actually
C sin[integral( k(s)ds)]. We set dsx(13,4) = K and this is 
c integrated each step. The spatial sin term in dsx(13,6) below 
c evaluates sx(13,4), which is the integral.
           dsx(13,4) = CAY



C**********Compute energy gain (2-term potential)*******************
      dsx(13,6)=(0.5)*(A10*CAY*VANEVOLT*charge)
     +*sspace*sss/PMASS
!*******************************************************************

c Note: Have to divide by pratio as space charge parameter is different
c when x' is momentum compared with p_x/p_i.

C Compute 2-term potential expansion matrix elements.


      T1QUOT = (charge*VANEVOLT)/(4*beta)
      T1BRACKET= 4.*A01 + cspace*A10*CAY**2
      T2BRACKET = 4.*A01 - cspace*A10*(CAY**2)
      BQUOT = (charge*VANEVOLT*A10)/(2.*BETA**2*gamma**3*PMASS)
      BBRACKET = (omoc/BETA)*cspace*ccc + 
     +CAY*sspace*sss

      CQUOT = (A10*(omoc**2)*cspace)/(4.*BETA**5*gamma**3
     +*PMASS)

      CBRACKET=(2*abs(charge)*VANEVOLT*BETA**2*gamma**3*PMASS*sss
     $-abs(charge)**2*VANEVOLT**2*A10*cspace*ccc**2)

      QSC=abs(CURRENT*CHARGE*898755.2/(ETA**2*PMASS)/pratio)
         F(1,2)=pratio
         F(2,1)=-T1QUOT*T1BRACKET*sss/Pold
         F(3,4)=pratio
         F(4,3)=T1QUOT*T2BRACKET*sss/Pold
         F(5,5)=BQUOT*BBRACKET
         F(5,6)=pratio/gsq
         F(6,5)=CQUOT*CBRACKET/Pold
         F(6,6)=-1.*F(5,5)

!****K MAPPING GENERATIOR*******************************************
! leave commented out, except for development.
!      write(20,'(5(E20.10E2))')Z-zse,A01,A10,CAY,
!     +CAY/(sx(13,5)*omoc/Z)
!USE WITH CASEVAR 3 FOR SMOOTHNESS**********************************

C output s,phase to fort.20 - observe: here we offset phase
C by +25deg, such that for phi = -25 (ISAC RFQSC design), phase
C will be zero.
c      write(20,'(8(E20.10E2))')Z-zse,
c     +(sx(13,4)-omoc*sx(13,5))*180/3.1415926

c     It is called even if QSC=0 because it changes F (the no-space-charge F matrix) to FT (the one with space charge).
      CALL F2FT_SC(z,SX,DSX,F,FT)
      CALL COMPUTE_DSX(SX,DSX)
      
      RETURN
      END
      
