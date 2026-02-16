c
c$debug
      program diffrad
c
c      version 1.0  01.06.1998
c
c      Igor Akushevich	aku@hep.by  or	akush@hermes.desy.de
c      National Center of Particle and High Energy Physics
c      Minsk, Belarus
c
c
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/vv1vv2/aa1,aa2,bb,sib
      common/ivv/vcurr,cutv
      external sigphi
      real*8 sigphi

      dimension xmas(200)
      dimension ymas(200)
      dimension tmas(200)
      dimension phimas(200)
      parameter (netx=1)
      parameter (nety=6)
      parameter (nett=2)
      dimension xnet(netx),ynet(nety),tnet(nett)
      data xnet/15d0/
     .	   ynet/-0.7d0,-1.5d0,-3d0,-5d0,-8d0,-11d0/
     .	   tnet/-.1d0,-.3d0/
c      data xnet/-10.1d0,-14.9d0,-19.3d0,-24.1d0/
c     .     ynet/-0.24d0,-.61d0,-1.84d0,-5.69d0/
c     .     tnet/-.25d0/
	 open(unit=8,file='input.dat',status='old')
	 read(8,*)bmom
	 read(8,*)tmom
	 read(8,*)lepton
	 read(8,*)ivec
	 read(8,*)intphi
	 read(8,*)cutv
	 read(8,*)npoi
	 if(npoi.ge.1)then
	  read(8,*)(xmas(i),i=1,npoi)
	  read(8,*)(ymas(i),i=1,npoi)
	  read(8,*)(tmas(i),i=1,npoi)
	  if(intphi.ne.1)read(8,*)(phimas(i),i=1,npoi)
	 else
	  do it=1,nett
	  do ix=1,netx
	  do iy=1,nety
	    npoi=npoi+1
	    xmas(npoi)=-xnet(ix)**2
	    ymas(npoi)=ynet(iy)
	    tmas(npoi)=tnet(it)
	  enddo
	  enddo
	  enddo
	 endif
	 close(8)


	 call titout('01.06.1998')
	 write(9,'(a8,f8.3)')' bmom  =',bmom
	 write(9,'(a8,f8.3)')' tmom  =',tmom
	 write(9,'(a8,i2)')' lepton=',lepton
	 write(9,'(a8,i2)')' ivec  =',ivec
	 write(9,'(a8,i2)')' intphi=',intphi
	 write(9,'(a8,f6.3)')' cutv  =',cutv
	 write(9,'(a8,i3)')' npoi  =',npoi

	 call setcon(ivec,lepton)

      snuc=2.*(sqrt(tmom**2+amp*amp)*sqrt(bmom**2+aml2)+bmom*tmom)

      do 1 i=1,npoi
	tdif=tmas(i)
	if(intphi.eq.0)phih=phimas(i)
	if(xmas(i).ge.0d0.and.ymas(i).ge.0d0)then
	  xs=xmas(i)
	  ys=ymas(i)
	  q2=snuc*xs*ys
	  w2=ys*snuc-q2+amp2
	elseif(xmas(i).ge.0d0.and.ymas(i).le.0d0)then
	  xs=xmas(i)
	  q2=-ymas(i)
	  ys=q2/(snuc*xs)
	  w2=ys*snuc-q2+amp2
	elseif(xmas(i).le.0d0.and.ymas(i).ge.0d0)then
	  w2=-xmas(i)
	  ys=ymas(i)
	  q2=ys*snuc-w2+amp2
	  xs=q2/snuc/ys
	elseif(xmas(i).le.0d0.and.ymas(i).le.0d0)then
	  w2=-xmas(i)
	  q2=-ymas(i)
	  ys=(w2+q2-amp2)/snuc
	  xs=q2/snuc/ys
	endif

c      write(*,'(4h x =,f6.3,7h    y =,f6.3,7h    Q2=,f8.3
c     . ,7h    W2=,f8.3,7h    t =,f8.3)')xs,ys,Q2,W2,tdif
      write(9,'(60(1h*))')
      write(9,'(4h x =,f6.3,7h    y =,f6.3,7h    Q2=,f8.3,7h    W2=,f8.3
     . ,7h    t =,f8.3)')xs,ys,Q2,W2,tdif
      yma=1d0/(1d0+amp**2*xs/snuc)
      ymi=(amc2-amp**2)/(snuc*(1d0-xs))
      if(ys.gt.yma.or.ys.lt.ymi.or.xs.gt.1d0.or.xs.lt.0d0)then
	write(9,*)' kinematics'
	goto 1
      endif
	 call conkin(snuc)

      if(sqrt(W2).lt.amp+amv)then
	write(9,*)' W < M + M_v; W =',sqrt(W2),' M + M_v =',amp+amv
	 goto 1
      endif

	 tt1=w2-q2-amp2
	 tt2=w2-amp2+amv**2
	 tdmin=-q2+amv**2-.5d0/W2*(tt1*tt2
     .	 +sqrt(tt1**2+4d0*q2*w2)*sqrt(tt2**2-4d0*amv**2*w2))
	 tdmax=-q2+amv**2-.5d0/W2*(tt1*tt2
     .	 -sqrt(tt1**2+4d0*q2*w2)*sqrt(tt2**2-4d0*amv**2*w2))


      if(tdif.lt.tdmin.or.tdif.gt.tdmax)then
	write(9,*)' t < t_min  or  t > t_max'
	write(9,'(6h tmin=,f8.3,4h t =,f8.3,6h tmax=,f8.5)')
     .				tdmin,tdif,tdmax
	goto 1
      else
	write(9,'(6h tmin=,f8.3,6h tmax=,f8.5)')tdmin,tdmax
      endif

c      vmax=2d0*amp*1d0

      sxt=sx+tdif
      tq=q2+tdif-amv**2
      aa1=(q2*Sxp*Sxt-(S*Sx+2d0*amp2*q2)*tq)/2.d0/aly
      aa2=(q2*Sxp*Sxt-(X*Sx-2d0*amp2*q2)*tq)/2.d0/aly
      sqbb1=sqrt(q2*Sxt**2-Sxt*Sx*tq-amp2*tq**2-amv**2*aly)
      sqbb2=sqrt(q2*(S*X-amp2*q2)-aml2*aly)
      bb=sqbb1*sqbb2/aly

      vmax=tt2+.5d0/q2*(-tt1*tq
     .	-sqrt(tt1**2+4d0*q2*w2)*sqrt(tq**2+4d0*amv**2*q2))
     -	 - 1d-8

      if(cutv.gt.1d-12)vmax=min(vmax,cutv)

	 call bornin(sib)

      if(cutv.lt.-1d-12)then
	 vcurr=0d0
	 if(intphi.eq.0)sig0=sigphi(phih)
	 if(intphi.eq.1)then
	  call dqg32(0d0,2d0*pi,sigphi,sig0)
	 endif
       ivvpo=15
c	vvvv=1.d0
       vvvv=vmax
       step=(vvvv+cutv)/ivvpo
       write(21,'(i4,3f8.4)')ivvpo,-cutv,step,sig0/sib/2d0/pi
       do ivv=1,ivvpo
	 vcurr=-cutv+step*(dble(ivv)-.5d0)
c	  vcurr=.15+0.1*(ivv-1)
	 if(intphi.eq.0)vdi=sigphi(phih)/sig0
	 if(intphi.eq.1)then
	  call dqg32(0d0,2d0*pi,sigphi,sig)
	  vdi=sig/sib/2d0/pi
	 endif
	 write(21,'(2f6.3,4f8.4)')xs,ys,tdif,vcurr,vdi
	 write(*,'(2f6.3,3f8.4)')xs,ys,tdif,vcurr,vdi
       enddo
       stop
      endif

	wei=1d0
       call approx(wei,wei0)

       sig=0.d0
       if(intphi.eq.0)sig=sigphi(phih)
       if(intphi.eq.1)then
	  call dqg32(0d0,2d0*pi,sigphi,sig)
	  sib=2d0*pi*sib
       endif

       write(9,'(2f6.3,3f8.3,4f9.3)')
     .	xs,ys,q2,w2,tdif,sig/sib,wei,abs(wei-sig/sib)
       write(21,'(2f6.3,f8.3,f8.0,f8.3,3f9.3,2g12.4)')
     .	xs,ys,q2,w2,tdif,sig/sib,wei,abs(wei-sig/sib),sig,sib
       write(*,'(2f6.3,3f8.3,4f9.3)')
     .	xs,ys,q2,w2,tdif,sig/sib,wei,abs(wei-sig/sib)


   1  continue
1000  end



****************** setcon *************************************

      subroutine setcon(ivec,lepton)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      dimension amhad(4)
      data amhad/0.7683d0,0.78195d0,1.019412d0,3.0969d0/

      if(lepton.eq.1)aml2=.261112d-6
      if(lepton.eq.2)aml2=.111637d-1
      pi=3.1415926d0
      alpha=.729735d-2
      barn=.389379d6
      amv=amhad(ivec)

      amp=.938272d0
      amp2=amp**2
      ap=2d0*amp
      ap2=2d0*amp2

      amc2=amp2

c	amc2=2.*amp*demin+amp2
      end

****************** sigphi *************************************

      double precision function sigphi(phih)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
       common/vv1vv2/aa1,aa2,bb,sib
      common/ivv/vcurr,cutv
      common/phi/phirad,tdif,phidif,tq,vmax,ivec

      phidif=phih

      vv1=(aa1+bb*cos(phidif))/2.d0
      vv2=(aa2+bb*cos(phidif))/2.d0

      sum=vacpol(q2)

      ssh=x+q2-vv2
      xxh=s-q2-vv1

       dlm=log(q2/aml2)

      deltavr=(1.5d0*dlm-2.d0-.5d0*log(xxh/ssh)**2
     .		    +fspen(1d0-amp2*q2/ssh/xxh)-pi**2/6.d0)

      if(cutv.gt.-1d-12)then
       call qqt(tai)
       delinf=(dlm-1.d0)*log(vmax**2/ssh/xxh)
       extai1=exp(alpha/pi*delinf)
       sigphi=sib*extai1*(1.d0+alpha/pi*(deltavr+sum))+tai

c	write(*,*)delinf,dlm,log(vmax**2/ssh/xxh)
c	write(*,'(f7.3,2g12.5,2f9.3)')phidif,sigphi,sib
c     .  ,tai/sib,sigphi/sib
       write(9,'(f7.3,2g12.5,2f9.3)')phidif,sigphi,sib
     .	,tai/sib,sigphi/sib

      elseif(vcurr.lt.1d-12)then
       delinf=(dlm-1.d0)*log(cutv**2/ssh/xxh)
       extai1=exp(alpha/pi*delinf)
       sigphi=sib*(1.d0+alpha/pi*(deltavr+sum))
      else
	call qqt(tai)
	sigphi=tai
      endif

      end


****************** conkin *************************************

      subroutine conkin(snuc)
c set of kinematical constants
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys

      s=snuc
      x=s*(1.-ys)
      sx=s-x
      sxp=s+x
      if(abs(w2-(amp2+s-q2-x)).gt.1d-10)stop 'w2'
      w2=amp2+s-q2-x
      aly=sx**2+4.*amp2*q2
      sqly=dsqrt(aly)
      anu=sx/ap
      axy=pi*(s-x)
c      an=2.*alpha**2/s*axy*barn
      an=alpha*ys/(8d0*pi**2)*barn
      tamax=(sx+sqly)/ap2
      tamin=-q2/amp2/tamax

      return
      end

****************** bornin *************************************

      subroutine bornin(sibor)
c
c     sibor is born cross section with polarized initial
c     lepton and polarized target
c     siamm is contribution of anomalous magnetic moment.
c
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/sigsig/sigmat,sigmal
      common/pri/ipri

       ga2=q2/anu**2

       ipri=1
       call difflt(q2,w2,tdif,sigmal,sigmat)
       ipri=0

       sibor=2.d0*an/(xs*ys**2)*(ys**2*sigmat+2d0*(1d0-ys-
     .	.25d0*ys**2*ga2)*(sigmal+sigmat))

       write(9,'(a8,3g12.4)')' sigmat ',sigmat*.389
       write(9,'(a8,3g12.4)')' bornin ',sigmal,sigmat,sibor
       write(9,'(a8,3g12.4)')' bornin2',an,xs,ga2

       end


****************  approx ************************************

      subroutine approx(wei,wei0)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/vv1vv2/aa1,aa2,bb,sib
      data slope/5d0/

      vv1=(aa1)/2.d0
      vv2=(aa2)/2.d0

      sum=vacpol(q2)

      ssh=x+q2-vv2
      xxh=s-q2-vv1

       dlm=log(q2/aml2)

      delinf=(dlm-1.d0)*log(vmax**2/ssh/xxh)

      delta=delinf+sum+(1.5d0*dlm-2.d0-.5d0*log(xxh/ssh)**2
     .		    +fspen(1d0-amp2*q2/ssh/xxh)-pi**2/6.d0)

      qv=q2+amv**2

      del_rad=2.*slope*vmax*(dlm-1.d0)*qv/(sx-qv)

      extai1=exp(alpha/pi*delinf)
      wei=extai1*(1.+alpha/pi*(delta-delinf+del_rad))
      wei0=alpha/pi*del_rad

      write(9,'(a12,2f7.4)')' wei, wei0 :',wei,wei0

      end


****************** qqt ****************************************

      subroutine qqt(tai)
      implicit real*8(a-h,o-z)
      external qqtphi
      real*8 qqtphi
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn

	call simpsx(0d0,2.d0*pi,150,1d-2,qqtphi,tai)
	tai=tai/2.d0

      end

****************** qqtphi *************************************

      double precision function qqtphi(phi)
      implicit real*8(a-h,o-z)
      external rv2ln
      real*8 rv2ln
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      dimension tlm(4)
      phirad=phi



      tlm(1)=log(xs+tamin)
      tlm(4)=log(xs+tamax)
      tlm(2)=log(xs-q2/s)
      tlm(3)=log(xs+q2/x)

      res=0d0

      do ii=1,3

      ep=1d-10
      call simptx(tlm(ii)+ep,tlm(ii+1)-ep,100,5d-3,rv2ln,re)
	     tai=an*alpha/pi*re
	     res=res+tai

       enddo

      qqtphi=res

      end



****************** rv2ln **************************************

      double precision function rv2ln(taln)
c
c     integrand (over ta )
c
      implicit real*8(a-h,o-z)
      external podinl
      real*8 podinl
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/amf2/taa,tm(8,6),sfm0(8)
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/ivv/vcurr,cutv


      ta=exp(taln)-q2/sx
      taa=ta

       sqrtmb=sqrt((ta-tamin)*(tamax-ta)*(s*x*q2-q2**2*amp2-aml2*aly))
       z1=(q2*sxp+ta*(s*sx+ap2*q2)-ap*cos(phirad)*sqrtmb)/aly
       z2=(q2*sxp+ta*(x*sx-ap2*q2)-ap*cos(phirad)*sqrtmb)/aly
       bb=1./sqly/pi
       bi12=bb/(z1*z2)
       bi1pi2=bb/z2+bb/z1
       bis=bb/z2**2+bb/z1**2
       bir=bb/z2**2-bb/z1**2
       hi2=aml2*bis-q2*bi12

      tm(1,1)=4.d0*q2*hi2
      tm(1,2)=4.d0*hi2*ta
      tm(1,3)=-2.d0*(2.d0*bb+bi12*ta**2)
      tm(2,1)=2d0*(S*X-amp2*q2)*hi2/amp2
      tm(2,2)=(2.*aml2*bir*sxp-4.*amp2*hi2*ta-bi12*sxp**2*ta+
     . bi1pi2*sxp*sx+2.*hi2*sx)/(2.*amp2)
      tm(2,3)=(2.*(2.*bb+bi12*ta**2)*amp2-bi12*
     . sx*ta-bi1pi2*sxp)/(2.*amp2)


       if(cutv.lt.-1d-12)then
	 res=podinl(vcurr)
       else
	 vmin=1d-4
	 call simpux(vmin,vmax,100,1d-3,podinl,res)
       endif

      rv2ln=res*(q2/sx+ta)

      end

****************** podinl *************************************

      double precision function podinl(v)
c
c     integrand (over r )
c
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/amf2/ta,tm(8,6),sfm0(8)
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/sigsig/sigmat0,sigmal0
      common/ivv/vcurr,cutv

      dimension sfm(8)

      sxtm=sx+tdif-v
      aak=((2d0*q2+ta*Sx)*sxtm-(sx-2d0*ta*amp2)*tq)/aly/2d0
      sqbbk1=sqrt(q2*Sxtm**2-Sxtm*Sx*tq-amp2*tq**2-amv**2*aly)
      sqbbk2=amp*sqrt((tamax-ta)*(ta-tamin))
      bbk=sqbbk1*sqbbk2/aly
      d2kvir=2d0*(aak+bbk*cos(phirad-phidif))
      factor=1d0+ta-d2kvir


      r=v/factor

      tldq2=q2+r*ta
      tldw2=w2-r*(1.d0+ta)
      tldtd=tdif-r*(ta-d2kvir)

       call difflt(tldq2,tldw2,tldtd,sigmal,sigmat)

      sfm(1)=(sx-r)*sigmat	     !*fg
      sfm(2)=2.d0*ap2/(sx-r)*tldq2*(sigmat+sigmal) !*fg
      sfm0(1)=sx*sigmat0	   !*fg
      sfm0(2)=2.d0*ap2/sx*q2*(sigmat0+sigmal0) !*fg

      podinl=0.
      do 11 isf=1,2
      do 1 irr=1,3
	pp=sfm(isf)
      if(irr.eq.1.and.cutv.gt.-1d-12)pp=pp-sfm0(isf)*(1.+r*ta/q2)**2
      pres=pp*r**(irr-2)/(q2+r*ta)**2/2.
      podinl=podinl-tm(isf,irr)*pres
    1 continue
   11 continue

      podinl=podinl/factor

      end


****************** difflt *************************************

      subroutine difflt(q2,w2,t,sigl,sigt)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/pri/ipri
      dimension ggam(4)
      data p02/.5d0/al_s/.25d0/
      data ggam/6.77d-6,0.6d-6,1.37d-6,5.36d-6/

      sx=w2+q2-amp2
      anu=sx/ap
      Sxt=sx+t
      eh=sxt/2d0/amp
      amv2=amv**2

      eta=1d0
      ff2=exp(5d0*t)
c      ff2=1d0/(1d0-t/0.71d0)**4

      pt20=eh**2-amv2-(t+q2-amv2+2d0*anu*eh)**2/4d0/(anu**2+q2)
      tqt=t+q2-amv2
      PT2=(-(4d0*(ANU**2+Q2)*AMV2+4d0*ANU*EH*TQt-4d0*EH**2*Q2+TQt**2)
     . )/(4d0*(ANU**2+Q2))

      if(pt2.lt.0d0.or.w2.lt.(amp+amv)**2)then
	write(*,'(7g11.4)')q2,w2,t,pt2,pt20
	if(pt2.lt.0d0)pt2=0d0
c	stop
      endif
      if(ipri.eq.1)write(9,*)pt2,sqrt(pt2)
      xsb=(q2+amv2+pt2)/w2
      q2b=(q2+amv2+pt2)/4d0

      if(pt2.le.p02)then
	 fm=log((4d0*q2b-pt2+p02)/(pt2+p02))
      else
	 fm=log((pt2+p02)/(4d0*q2b-pt2+p02)*4d0*q2b**2/pt2**2)
      endif

      xsbgm=3d0*(1d0-xsb)**5
      sk=xsbgm*fm/(2d0*q2b*(2d0*q2b-pt2)*log(8d0*q2b/p02))
      sigt=al_s**2*ggam(ivec)*amv**3/3d0/alpha*pi**3*sk**2*ff2*eta**2
      sigl=q2/amv2*sigt
      end



****************** vacpol *************************************

      double precision function vacpol(t)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      dimension am2(3)
c
c    am2 : squared masses of charge leptons
c
      data am2/.26110d-6,.111637d-1,3.18301d0/

      suml=0.
      do 10 i=1,3
	 a2=2.*am2(i)
	 sqlmi=dsqrt(t*t+2.*a2*t)
	 allmi=dlog((sqlmi+t)/(sqlmi-t))/sqlmi
  10  suml=suml+2.*(t+a2)*allmi/3.-10./9.+4.*a2*(1.-a2*allmi)/3./t
      if(t.lt.1.d0)then
	aaa = -1.345d-9
	bbb = -2.302d-3
	ccc = 4.091
      elseif(t.lt.64d0)then
	aaa = -1.512d-3
	bbb =  -2.822d-3
	ccc = 1.218
      else
	aaa = -1.1344d-3
	bbb = -3.0680d-3
	ccc = 9.9992d-1
      endif
      sumh = -(aaa+bbb*log(1.+ccc*t)) *2*pi/alpha

      vacpol=suml+sumh

      end

****************** fspens *************************************

      double precision function fspens(x)
c
c    spence function
c
      implicit real*8(a-h,o-z)
      f=0.d0
      a=1.d0
      an=0.d0
      tch=1.d-16
  1   an=an+1.d0
      a=a*x
      b=a/an**2
      f=f+b
      if(b-tch)2,2,1
  2   fspens=f
      return
      end

      double precision function fspen(x)
      implicit real*8(a-h,o-z)
      data f1/1.644934d0/
      if(x)8,1,1
  1   if(x-.5d0)2,2,3
    2 fspen=fspens(x)
      return
    3 if(x-1d0)4,4,5
    4 fspen=f1-dlog(x)*dlog(1d0-x+1d-10)-fspens(1d0-x)
      return
    5 if(x-2d0)6,6,7
    6 fspen=f1-.5*dlog(x)*dlog((x-1d0)**2/x)+fspens(1d0-1d0/x)
      return
    7 fspen=2d0*f1-.5d0*dlog(x)**2-fspens(1d0/x)
      return
    8 if(x+1d0)10,9,9
   9  fspen=-.5d0*dlog(1d0-x)**2-fspens(x/(x-1d0))
      return
  10  fspen=-.5*dlog(1.-x)*dlog(x**2/(1d0-x))-f1+fspens(1d0/(1d0-x))
      return
      end


************ titout **************************
      subroutine titout(date)
      implicit real*8(a-h,o-z)
      character*10 date

      open(unit=9,file='all.dat')
      open(unit=21,file='allu.dat')
	 write(9,1)date
	 write(9,3)
3     format(/' the file contains work information about'/
     .' contribution of tails')
1     format(1x,'program DIFFRAD 1.0 version from ',a10)
      end


      subroutine simpsx(a,b,np,ep,func,res)
      implicit real*8(a-h,o-z)
      external func
      step=(b-a)/np
      call simps(a,b,step,ep,1d-18,func,ra,res,r2,r3)
      end
      subroutine simptx(a,b,np,ep,func,res)
      implicit real*8(a-h,o-z)
      external func
      step=(b-a)/np
      call simpt(a,b,step,ep,1d-18,func,ra,res,r2,r3)
      end
      subroutine simpux(a,b,np,ep,func,res)
      implicit real*8(a-h,o-z)
      external func
      step=(b-a)/np
      call simpu(a,b,step,ep,1d-18,func,ra,res,r2,r3)
      end

      subroutine dqg32(xl,xu,fct,y)
c
c  computation of integrals by means of 32-point gauss quadrature
c  formula, which integrates polynomials up to degree 63.
c
c
      double precision xl,xu,y,a,b,c,fct
c      external fct
c
      a=.5d0*(xu+xl)
      b=xu-xl
      c=.49863193092474078d0*b
      y=.35093050047350483d-2*(fct(a+c)+fct(a-c))
      c=.49280575577263417d0*b
      y=y+.8137197365452835d-2*(fct(a+c)+fct(a-c))
      c=.48238112779375322d0*b
      y=y+.12696032654631030d-1*(fct(a+c)+fct(a-c))
      c=.46745303796886984d0*b
      y=y+.17136931456510717d-1*(fct(a+c)+fct(a-c))
      c=.44816057788302606d0*b
      y=y+.21417949011113340d-1*(fct(a+c)+fct(a-c))
      c=.42468380686628499d0*b
      y=y+.25499029631188088d-1*(fct(a+c)+fct(a-c))
      c=.39724189798397120d0*b
      y=y+.29342046739267774d-1*(fct(a+c)+fct(a-c))
      c=.36609105937014484d0*b
      y=y+.32911111388180923d-1*(fct(a+c)+fct(a-c))
      c=.33152213346510760d0*b
      y=y+.36172897054424253d-1*(fct(a+c)+fct(a-c))
      c=.29385787862038116d0*b
      y=y+.39096947893535153d-1*(fct(a+c)+fct(a-c))
      c=.25344995446611470d0*b
      y=y+.41655962113473378d-1*(fct(a+c)+fct(a-c))
      c=.21067563806531767d0*b
      y=y+.43826046502201906d-1*(fct(a+c)+fct(a-c))
      c=.16593430114106382d0*b
      y=y+.45586939347881942d-1*(fct(a+c)+fct(a-c))
      c=.11964368112606854d0*b
      y=y+.46922199540402283d-1*(fct(a+c)+fct(a-c))
      c=.7223598079139825d-1*b
      y=y+.47819360039637430d-1*(fct(a+c)+fct(a-c))
      c=.24153832843869158d-1*b
      y=b*(y+.48270044257363900d-1*(fct(a+c)+fct(a-c)))
      return
      end

      subroutine simps(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
c simps
c a1,b1 -the limits of integration
c h1 -an initial step of integration
c reps1,aeps1 - relative and absolute precision of integration
c funct -a name of function subprogram for calculation of integrand +
c x - an argument of the integrand
c ai - the value of integral
c aih- the value of integral with the step of integration
c aiabs- the value of integral for module of the integrand
c this subrogram calculates the definite integral with the relative or
c absolute precision by simpson+s method with the automatical choice
c of the step of integration
c if aeps1    is very small(like 1.e-17),then calculation of integral
c with reps1,and if reps1 is very small (like 1.e-10),then calculation
c of integral with aeps1
c when aeps1=reps1=0. then calculation with the constant step h1
c
      implicit real*8(a-h,o-z)
      dimension f(7),p(5)
      h=dsign(h1,b1-a1)
      s=dsign(1.d0,h)
      a=a1
      b=b1
      ai=0.d0
      aih=0.d0
      aiabs=0.d0
      p(2)=4.d0
      p(4)=4.d0
      p(3)=2.d0
      p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=dabs(reps1)
      aeps=dabs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a
      c=0.d0
      f(1)=funct(x)/3.d0
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.
      if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1)
      di3=dabs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.*f(3)+f(5))*2.*h
      di2=di2*h
      di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=dabs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=dabs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.) 17,14,14
   17 h=2.*h
      f(1)=f(5)
      f(2)=f(6)
      f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5)
      f(3)=f(6)
      f(5)=f(7)
      f(2)=10.d16
      f(4)=10.d16
      f(6)=10.d16
      f(7)=10.d16
   18 di1=di2+(di2-di1)/15.
      ai=ai+di1
      aih=aih+di2
      aiabs=aiabs+di3
      go to 22
   21 h=h/2.
      f(7)=f(5)
      f(6)=f(4)
      f(5)=f(3)
      f(3)=f(2)
      f(2)=10.d16
      f(4)=10.d16
      x=x0
      c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end


      subroutine simpt(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
      implicit real*8(a-h,o-z)
      dimension f(7),p(5)
      h=dsign(h1,b1-a1)
      s=dsign(1.d0,h)
      a=a1
      b=b1
      ai=0.d0
      aih=0.d0
      aiabs=0.d0
      p(2)=4.d0
      p(4)=4.d0
      p(3)=2.d0
      p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=dabs(reps1)
      aeps=dabs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a
      c=0.d0
      f(1)=funct(x)/3.
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.
      if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1)
      di3=dabs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.*f(3)+f(5))*2.*h
      di2=di2*h
      di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=dabs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=dabs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.) 17,14,14
   17 h=2.*h
      f(1)=f(5)
      f(2)=f(6)
      f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5)
      f(3)=f(6)
      f(5)=f(7)
      f(2)=10.d16
      f(4)=10.d16
      f(6)=10.d16
      f(7)=10.d16
   18 di1=di2+(di2-di1)/15.
      ai=ai+di1
      aih=aih+di2
      aiabs=aiabs+di3
      go to 22
   21 h=h/2.
      f(7)=f(5)
      f(6)=f(4)
      f(5)=f(3)
      f(3)=f(2)
      f(2)=10.d16
      f(4)=10.d16
      x=x0
      c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end

      subroutine simpu(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
      implicit real*8(a-h,o-z)
      dimension f(7),p(5)
      h=dsign(h1,b1-a1)
      s=dsign(1.d0,h)
      a=a1
      b=b1
      ai=0.d0
      aih=0.d0
      aiabs=0.d0
      p(2)=4.d0
      p(4)=4.d0
      p(3)=2.d0
      p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=dabs(reps1)
      aeps=dabs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a
      c=0.d0
      f(1)=funct(x)/3.
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.
      if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1)
      di3=dabs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.*f(3)+f(5))*2.*h
      di2=di2*h
      di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=dabs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=dabs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.) 17,14,14
   17 h=2.*h
      f(1)=f(5)
      f(2)=f(6)
      f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5)
      f(3)=f(6)
      f(5)=f(7)
      f(2)=10.d16
      f(4)=10.d16
      f(6)=10.d16
      f(7)=10.d16
   18 di1=di2+(di2-di1)/15.
      ai=ai+di1
      aih=aih+di2
      aiabs=aiabs+di3
      go to 22
   21 h=h/2.
      f(7)=f(5)
      f(6)=f(4)
      f(5)=f(3)
      f(3)=f(2)
      f(2)=10.d16
      f(4)=10.d16
      x=x0
      c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end

