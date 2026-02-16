c$debug
      program mdiffrad
c
c      version 1.0  01.09.1998
c
c      Igor Akushevich	aku@hep.by  or	akush@hermes.desy.de
c      National Center of Particle and High Energy Physics
c      Minsk, Belarus
c
c
      implicit real*8(a-h,o-z)
      parameter(ntmax=100)
      dimension xmas(ntmax,2),ymas(ntmax,2),tmas(ntmax,2)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn,bslope
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/phi/phirad,tdif,phidif,tq,vmax,vmmm,cutv,ivec
      parameter(maxphi=3)
      parameter(maxta=4)
      dimension phrmas(maxphi+1)
      dimension tamas(maxta+1)
      parameter(nqua=1)
      dimension resu(0:nqua,3),resu2(0:nqua,3)
      dimension resig(0:nqua,3),reerr(0:nqua,3)
      real*4 urand,rndm
      data phrmas/0d0,0.001d0,0.999d0,1d0/
c      urand(nn)=rndm(-1)
      open(unit=9,file='ALLmc.dat')
      open(unit=19,file='etamc.dat')
      open(unit=20,file='res20.dat')
      open(unit=21,file='res21.dat')
	 open(unit=8,file='inmdi.dat',status='old')
	 read(8,*)bmom
	 read(8,*)tmom
	 read(8,*)lepton
	 read(8,*)ivec
	 read(8,*)ann1,ann2,ann3
	 read(8,*)cutv
	 read(8,*)nev
	 read(8,*)iy
	 read(8,*)npoi
	 if(npoi.gt.ntmax)stop 'npoi.gt.ntmax'
      write(19,'(a7,f8.3)')' bmom =',bmom
      write(19,'(a7,f8.3)')' tmom =',tmom
      write(19,'(a7,i3)')' lepton =',lepton
      write(19,'(a7,i3)')' ivec =',ivec
      write(19,'(a7,f8.3)')' cutv =',cutv
      write(19,'(a7,i5)')' nev =',nev
      write(19,'(a7,i8)')' iy :  ',iy
      write(19,'(a7,i3)')' npoi =',npoi
	  read(8,*)(xmas(i,1),i=1,npoi)
	  read(8,*)(xmas(i,2),i=1,npoi)
	  read(8,*)(ymas(i,1),i=1,npoi)
	  read(8,*)(ymas(i,2),i=1,npoi)
	  read(8,*)(tmas(i,1),i=1,npoi)
	  read(8,*)(tmas(i,2),i=1,npoi)


      call setcon(ivec,lepton)
      s=2.*(sqrt(tmom**2+amp*amp)*sqrt(bmom**2+aml2)+bmom*tmom)


      do i1=1,npoi

	sumnev=0d0
	sumnev2=0d0
	err=0d0

	do jjiy=1,nev

	do ii1=0,nqua
	do ii2=1,3
	    resu(ii1,ii2)=0.
	    resu2(ii1,ii2)=0.
	    resig(ii1,ii2)=0.
	    reerr(ii1,ii2)=0.
	enddo
	enddo

	keyrad=0
	nmax1=int(ann1+0.01)
	do i=1,nmax1
	   call rankin(i1,ntmax,xmas,ymas,tmas,iy,sg,iacc)
	   if(iacc.eq.1)then
	     resadd=sg*sig_ep(keyrad)
	   else
	     resadd=0d0
	   endif
	   resu(0,1)=resu(0,1)+resadd
	   resu2(0,1)=resu2(0,1)+resadd**2
	   resu(1,1)=resu(1,1)+sg
	   resu2(1,1)=resu2(1,1)+sg**2
c	  if(mod(i,100000).eq.0)print *,'  --> 1:',i
	enddo
	do ii1=0,nqua
	  resig(ii1,1)=resu(ii1,1)/nmax1
	  errim=resu2(ii1,1)/nmax1-resig(ii1,1)**2
	  if(errim.lt.0d0)errim=0d0
	  reerr(ii1,1)=sqrt(errim/nmax1)
	enddo

	keyrad=1
	nmax1=int(ann2+0.01)
	res=0d0
	res2=0d0
	do i=1,nmax1
	   call rankin(i1,ntmax,xmas,ymas,tmas,iy,sg,iacc)
	   if(iacc.eq.1)then
	     resadd=sg*sig_ep(keyrad)
	   else
	     resadd=0d0
	   endif
	   resu(0,2)=resu(0,2)+resadd
	   resu2(0,2)=resu2(0,2)+resadd**2
	   resu(1,2)=resu(1,2)+sg
	   resu2(1,2)=resu2(1,2)+sg**2
cc	  if(mod(i,100000).eq.0)print *,'  --> 2:',i
	enddo
	do ii1=0,nqua
	  resig(ii1,2)=resu(ii1,2)/nmax1
	  errim=resu2(ii1,2)/nmax1-resig(ii1,2)**2
	  if(errim.lt.0d0)errim=0d0
	  reerr(ii1,2)=sqrt(errim/nmax1)
	enddo

	nmax1=int(ann3+0.01)
	do iphir=1,maxphi
	do ita=1,maxta
	   do ii1=0,nqua
	   do ii2=1,3
	      resu(ii1,ii2)=0.
	      resu2(ii1,ii2)=0.
	   enddo
	   enddo
	 do i=1,nmax1
	   call rankin(i1,ntmax,xmas,ymas,tmas,iy,sg,iacc)
c	   phirad=2d0*pi*urand(iy)
c	   ta=(sx-(2d0*urand(iy)-1d0)*sqly)/ap2

	   phrmin=2d0*pi*phrmas(iphir)
	   phrmax=2d0*pi*phrmas(iphir+1)
	   phirad=phrmin+(phrmax-phrmin)*urand(iy)

	   tamas(1)=tamin
	   tamas(2)=(tamin-ys*xs)/2d0
	   tamas(3)=0d0
	   tamas(4)=1d0
	   tamas(5)=tamax
	   taumin=tamas(ita)
	   taumax=tamas(ita+1)
	   ta=taumin+(taumax-taumin)*urand(iy)

	   v=vmax*urand(iy)
	   sg=sg * (taumax-taumin) * vmax * (phrmax-phrmin)
	   if(iacc.eq.1)then
	     resadd=sig_rad(v,ta)*sg
	   else
	     resadd=0d0
	   endif
	   resu(0,3)=resu(0,3)+resadd
	   resu2(0,3)=resu2(0,3)+resadd**2
	   resu(1,3)=resu(1,3)+sg
	   resu2(1,3)=resu2(1,3)+sg**2
c	  if(mod(i,100000).eq.0)print *,'  --> 3:',i
	enddo

	do ii1=0,nqua
	  sirii=resu(ii1,3)/nmax1
	  resig(ii1,3)=resig(ii1,3)+sirii
	  errim=resu2(ii1,3)/nmax1-sirii**2
	  if(errim.lt.0d0)errim=0d0
	  errii=sqrt(errim/nmax1)
	  reerr(ii1,3)=sqrt(reerr(ii1,3)**2+errii**2)
	enddo

       enddo
       enddo

	dev=resig(0,2)/resig(0,1)
	der=resig(0,3)/resig(0,1)
	de=dev+der

	sumnev=sumnev+de
	sumnev2=sumnev2+de**2
	summea=sumnev/jjiy
	if(jjiy.ne.1)err=sqrt((sumnev2/jjiy-summea**2)/jjiy)

c	 write(19,'(2f9.4)')sir
c	 write(*,'(2f9.4)')sir
	write(*,'(a5,g11.4,6f7.3)')' main',
     .	resig(0,1)/resig(1,1)*2d0*pi,dev,der,de,sumnev/jjiy,err
	write(9,'(a5,g11.4,6f7.3)')' main',
     .	resig(0,1)/resig(1,1)*2d0*pi,dev,der,de,sumnev/jjiy,err

	do ii1=0,nqua
	  resdel=(resig(ii1,2)+resig(ii1,3))/resig(ii1,1)
	  errdel=sqrt(reerr(ii1,2)**2+reerr(ii1,3)**2
     .		 +resdel**2*reerr(ii1,3)**2)/resig(ii1,1)
	  lun=20+ii1
	  write(lun,'(2f7.3,6g11.4)') resdel,errdel
     .	  ,resig(ii1,1),reerr(ii1,1)
     .	  ,resig(ii1,2),reerr(ii1,2)
     .	  ,resig(ii1,3),reerr(ii1,3)
	enddo

      enddo
	write(19,'(i5,5f7.3)')i1,summea,err

10    enddo

      end





****************** setcon *************************************

      subroutine setcon(ivec,lepton)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn,bslope
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
      bslope=5d0
c	amc2=2.*amp*demin+amp2
      end

****************** rankin *************************************

      subroutine rankin(i,npoi,xmas,ymas,tmas,iy,sg,iacc)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn,bslope
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/phi/phirad,tdif,phidif,tq,vmax,vmmm,cutv,ivec
      dimension xmas(npoi,2),ymas(npoi,2),tmas(npoi,2)
      common/aa12b/aa1,aa2,bb
      real*4 urand,rndm

c     urand(nn)=rndm(-1)
      data qslope/1d0/

      iacc=1
      xmasmin=min(abs(xmas(i,1)),abs(xmas(i,2)))
      xmasmax=max(abs(xmas(i,1)),abs(xmas(i,2)))
      ymasmin=min(abs(ymas(i,1)),abs(ymas(i,2)))
      ymasmax=max(abs(ymas(i,1)),abs(ymas(i,2)))
	if(xmas(i,1).le.0d0 .and. xmas(i,2).le.0d0.and.
     .	   ymas(i,1).le.0d0 .and. ymas(i,2).le.0d0)then
	   w2min=(amp+amv)**2
	   w2min=max(w2min,xmasmin)
	   smm=s+aml2+amp2
	   w2max=(sqrt(smm)-sqrt(aml2))**2
	   if(xmasmax.gt.1d-12)w2max=min(w2max,xmasmax)
	   if(w2max.lt.w2min)then
	      iacc=0
c	      write(*,'(1x,2hw ,2g11.2)')w2max,w2min
	      sg=0d0
	      return
	   endif
	   absx=w2max-w2min
	   w2=w2min+urand(iy)*absx

	   sw=s+amp2-w2
	   dlsw2=sw**2-4d0*w2*aml2
	   q2max=.5d0/smm*(s*sw-2d0*aml2*(w2+amp2)+
     .	     sqrt((s**2-4d0*aml2*amp2)*dlsw2))
	   q2min=(w2-amp2)**2*aml2/smm/q2max
	   q2min=max(q2min,ymasmin)
	   if(ymasmax.gt.1d-12)q2max=min(q2max,ymasmax)
	   if(q2max.lt.q2min)then
	      iacc=0
c	      write(*,'(1x,2hq ,2g11.2)')q2max,q2min
	      sg=0d0
	      return
	   endif
c	    q1exp=exp(-qslope*q2min)
c	    q2exp=exp(-qslope*q2max)
c	    absy=q1exp-q2exp
c	    ranlog=q2exp+urand(iy)*absy
c	    q2=-1d0/qslope*log(ranlog)
	  if(min(ymasmin,ymasmax).lt.1d-8)then
	     q1log=log(q2min)
	     q2log=log(q2max)
c	     q2=exp(q1log+urand(iy)*(q2log-q1log))
	    q2=q2min*(q2max/q2min)**urand(iy)
	    absydep=q2*(q2log-q1log)
	  else
	    absydep=q2max-q2min
	    q2=q2min+urand(iy)*absydep
	  endif
	  ys=(w2+q2-amp2)/s
	  xs=q2/s/ys
c	  cojk=ys*s**2
	  cojk=1.
	else
	 absx=abs(xmas(i,1)-xmas(i,2))
	 absy=abs(ymas(i,1)-ymas(i,2))
	 xmasi=min(abs(xmas(i,1)),abs(xmas(i,2)))+urand(iy)*absx
	 ymasi=min(abs(ymas(i,1)),abs(ymas(i,2)))+urand(iy)*absy
	endif
	if(xmas(i,1).gt.0d0.and.ymas(i,1).gt.0d0)then
	  xs=xmasi
	  ys=ymasi
	  q2=s*xs*ys
	  w2=ys*s-q2+amp2
	  cojk=1d0
	elseif(xmas(i,1).gt.0d0.and.ymas(i,1).lt.0d0)then
	  xs=xmasi
	  q2=ymasi
	  ys=q2/(s*xs)
	  w2=ys*s-q2+amp2
c	  cojk=s*xs
	  cojk=1d0
	elseif(xmas(i,1).lt.0d0.and.ymas(i,1).gt.0d0)then
	  w2=xmasi
	  ys=ymasi
	  q2=ys*s-w2+amp2
	  xs=q2/s/ys
c	  cojk=ys*s
	  cojk=1d0
	endif

	 tt1=w2-q2-amp2
	 tt2=w2-amp2+amv**2
	 tdmink=-q2+amv**2-.5d0/W2*(tt1*tt2
     .	 +sqrt(tt1**2+4d0*q2*w2)*sqrt(tt2**2-4d0*amv**2*w2))
	 tdmaxk=(q2+amv**2)**2*amp2/w2/tdmink

	       !   tdmink   <  tdmaxk	<  0
	       !   tmasmink <= tmasmaxk <= 0

	tmasmin=-max(abs(tmas(i,2)),abs(tmas(i,1)))
	tmasmax=-min(abs(tmas(i,1)),abs(tmas(i,2)))
	tdmax=min(tmasmax,tdmaxk)
	if(tmasmin.lt.-1d-10)then
	   tdmin=max(tmasmin,tdmink)
	else
	   tdmin=tdmink
	endif
	if(tdmax.lt.tdmin)then
	  iacc=0
c	      write(*,'(1x,2ht ,2g11.2)')tdmax,tdmin
	      sg=0d0
	  return
	endif


	bt2exp=exp(bslope*tdmin)
	bt1exp=exp(bslope*tdmax)
	abst=(bt1exp-bt2exp)
	ranexp=bt2exp+urand(iy)*abst

	tdif=log(ranexp)/bslope

c	write(*,'(6g12.4)')tdmin,tdmax,tdif,log(-tdmin),log(-tdmax)
c    .,log(-tdif)
c	write(*,'(6g12.4)')q2min,q2max,q2,log(q2min),log(q2max)
c     .,log(q2)

	phidif=2d0*pi*urand(iy)


      if(tdif.lt.tdmin.or.tdif.gt.tdmax)then
c	print *,tmas(i,1),tmas(i,2)
       print *,tdmin,tdmax,tdif
c	print *,max(min(tmas(i,2),tmas(i,1)),tdmin)
c	print *,min(max(tmas(i,1),tmas(i,2)),tdmax)
c	print *,expbt2,expbt1,expran
c	pause
	write(9,*)' t < t_min  or  t > t_max'
	write(9,'(6h tmin=,f8.3,4h t =,f8.3,6h tmax=,f8.5)')
     .				tdmin,tdif,tdmax
       if(tdif.lt.tdmin)tdif=tdmin+1d-10
       if(tdif.gt.tdmax)tdif=tdmax-1d-10
c	 stop
      else
c	 write(9,'(6h tmin=,g11.3,6h tmax=,g11.3)')tdmin,tdmax
c	 write(9,'(6h tmin=,f8.3,4h t =,f8.3,6h tmax=,f8.5)')
c     . 			 tdmin,tdif,tdmax
      endif

      sg=2d0*pi
      if(absx.gt.1d-10)sg=sg * absx
      if(abs(ymas(i,1)-ymas(i,2)).gt.1d-10.or.
     .	abs(ymas(i,1)+ymas(i,2)).lt.1d-10)
     . sg=sg * absydep	! /qslope /exp(-qslope*q2)
      if(abs(tmas(i,2)-tmas(i,1)).gt.1d-10)
     . sg=sg * abst/bslope /exp(bslope*tdif)

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
      an=cojk*alpha*ys/(8d0*pi**2)*barn
      tamax=(sx+sqly)/ap2
      tamin=-q2/amp2/tamax

      sxt=sx+tdif
      tq=q2+tdif-amv**2
      aa1=(q2*Sxp*Sxt-(S*Sx+2d0*amp2*q2)*tq)/2.d0/aly
      aa2=(q2*Sxp*Sxt-(X*Sx-2d0*amp2*q2)*tq)/2.d0/aly
      ttptt=(tdmaxk-tdif)*(tdif-tdmink)
      bb1_old=q2*Sxt**2-Sxt*Sx*tq-amp2*tq**2-amv**2*aly
      bb1=w2*ttptt
c      print *,bb1,bb1_old
c      pause
      if(bb1.lt.0d0)then
	write(*,*)' bb1 =',bb1
	bb1=0d0
      endif
      sqbb1=sqrt(bb1)
      bb2=q2*(S*X-amp2*q2)-aml2*aly
      if(bb2.lt.0d0)then
	write(*,*)' bb2 =',bb2
	write(*,*)' sxy',s,x,q2
	bb2=0d0
      endif
      sqbb2=sqrt(bb2)
      bb=sqbb1*sqbb2/aly


      vmmm=tt2+.5d0/q2*(-tt1*tq
     .	+sqrt(tt1**2+4d0*q2*w2)*sqrt(tq**2+4d0*amv**2*q2))

      vmax=w2/q2*ttptt/vmmm

      vmax_old=tt2+.5d0/q2*(-tt1*tq
     .	-sqrt(tt1**2+4d0*q2*w2)*sqrt(tq**2+4d0*amv**2*q2))

c      print *,vmax,vmax_old,vmmm
c      pause
c     -   - 1d-8
      if(cutv.gt.1d-12)vmax=min(vmax,cutv)
       end


**********************************************************
      real*8 function sig_ep(keyrad)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn,bslope
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/phi/phirad,tdif,phidif,tq,vmax,vmmm,cutv,ivec
      common/aa12b/aa1,aa2,bb

       ga2=q2/anu**2

       ipri=1
       call difflt(q2,w2,tdif,sigmal,sigmat)
       ipri=0

       sig_ep=2.d0*an/(xs*ys**2)*(ys**2*sigmat+2d0*(1d0-ys-
     .	.25d0*ys**2*ga2)*(sigmal+sigmat))

      if(keyrad.eq.1)then

      vv1=(aa1+bb*cos(phidif))/2.d0
      vv2=(aa2+bb*cos(phidif))/2.d0

      sum=vacpol(q2)

      ssh=x+q2-vv2
      xxh=s-q2-vv1

       dlm=log(q2/aml2)

      deltavr=(1.5d0*dlm-2.d0-.5d0*log(xxh/ssh)**2
     .		    +fspen(1d0-amp2*q2/ssh/xxh)-pi**2/6.d0)

       delinf=(dlm-1.d0)*log(vmax**2/ssh/xxh)
       extai1=exp(alpha/pi*delinf)
	sig_ep=sig_ep*extai1*(1d0+alpha/pi*(deltavr+sum))

c	icou=icou+1
c	write(*,'(i5)')icou
      endif
      end

**********************************************************
      real*8 function sig_rad(v,ta)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn,bslope
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/phi/phirad,tdif,phidif,tq,vmax,vmmm,cutv,ivec
      dimension tm(2,3),sfm(2),sfm0(2)

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


      sxtm=sx+tdif-v
      aak=((2d0*q2+ta*Sx)*sxtm-(sx-2d0*ta*amp2)*tq)/aly/2d0
c      bbk1_old=q2*Sxtm**2-Sxtm*Sx*tq-amp2*tq**2-amv**2*aly
c	tv=v-tdif
c	tqm=q2+amv**2-tdif
c	bbk1_old2=q2*tv**2-tv*sx*tqm-amp2*tqm**2-tdif*aly
c
	bbk1=q2*(vmax-v)*(vmmm-v)
c	print *,bbk1_old,bbk1_old2,bbk1
c	pause
      if(bbk1.lt.0d0)then
	 write(*,*)' bbk1 =',bbk1
	 bbk1=0d0
      endif
      sqbbk1=sqrt(bbk1)
      bbk2=amp2*(tamax-ta)*(ta-tamin)
      if(bbk2.lt.0d0)then
	 write(*,*)' bbk2 =',bbk2
	 write(*,*)'ta',tamin,ta,tamax
	 bbk2=0d0
      endif
      sqbbk2=sqrt(bbk2)
      bbk=sqbbk1*sqbbk2/aly
      d2kvir=2d0*(aak+bbk*cos(phirad-phidif))
      factor=1d0+ta-d2kvir


      r=v/factor

      tldq2=q2+r*ta
      tldw2=w2-r*(1.d0+ta)
      tldtd=tdif-r*(ta-d2kvir)

       call difflt(q2,w2,tdif,sigmal0,sigmat0)
       call difflt(tldq2,tldw2,tldtd,sigmal,sigmat)

c     write(*,'(5g11.3)')r,ta,v,factor
c     write(*,'(5g11.3)')q2,w2,tdif,sigmal0,sigmat0
c     write(*,'(5g11.3)')tldq2,tldw2,tldtd,sigmal,sigmat

      sfm(1)=(sx-r)*sigmat	     !*fg
      sfm(2)=2.d0*ap2/(sx-r)*tldq2*(sigmat+sigmal) !*fg
      sfm0(1)=sx*sigmat0	   !*fg
      sfm0(2)=2.d0*ap2/sx*q2*(sigmat0+sigmal0) !*fg

      podinl=0.
      do 11 isf=1,2
      do 1 irr=1,3
	pp=sfm(isf)
      if(irr.eq.1)pp=pp-sfm0(isf)*(1.d0+r*ta/q2)**2
      pres=pp*r**(irr-2)/(q2+r*ta)**2/2.
      podinl=podinl-tm(isf,irr)*pres
c     write(*,'(2i4,5g11.3)')isf,irr,pp,pres,tm(isf,irr),podinl
    1 continue
   11 continue

	sig_rad=.5d0*an*alpha/pi*podinl/factor

cc    write(*,'(5g11.3)')sfm,sfm0
cc    write(*,'(5g11.3)')sfm,sfm0
c     write(*,'(5g11.3)')sig_rad,alpha,an,pi,podinl

      end


****************** difflt *************************************

      subroutine difflt(q2,w2,t,sigl,sigt)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn,bslope
      common/phi/phirad,tdif,phidif,tq,vmax,vmmm,cutv,ivec
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
c      ff2=1d0
      ff2=exp(bslope*t)
c      ff2=1d0/(1d0-t/0.71d0)**4

      pt20=eh**2-amv2-(t+q2-amv2+2d0*anu*eh)**2/4d0/(anu**2+q2)
      tqt=t+q2-amv2
      PT2=(-(4d0*(ANU**2+Q2)*AMV2+4d0*ANU*EH*TQt-4d0*EH**2*Q2+TQt**2)
     . )/(4d0*(ANU**2+Q2))

      if(pt2.lt.0d0)then
	write(*,'(a5,7g11.4)')' pt2 ',q2,w2,t,pt2,pt20
	if(pt2.lt.0d0)pt2=0d0
c	stop
      endif
      if(w2.lt.(amp+amv)**2)then
	write(*,'(a5,7g11.4)')' w2  ',w2,(amp+amv)**2,pt2
	stop
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
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn,bslope
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

c$nodebug
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C *						      *
      FUNCTION URAND(IY)
C *						      *
C *   This is a standard pseudo-random generator      *
C *   that work on IBM-370 and IBM-PC. We don't       *
C *   know does it work on SUN? 		      *
C *						      *
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
      INTEGER*4 IY,M2,IA,IC
      DATA S,M2,IA,IC/.46566128E-9,1073741824,843314861,453816693/
      IY=IY*IA+IC
      IF(IY.LT.0)IY=(IY+M2)+M2
      URAND=FLOAT(IY)*S
cc	urand=0.1
      END

*C * * * * * * * * * * * * * * * * * * * * * * * * * * *
*C *						       *
*      FUNCTION EXPVAL(AM,S)
*C *						       *
*C *   This is a simple generator of the pseudo        *
*C *   random numbers with Normal distribution.        *
*C *   AM - mean				       *
*C *   S  - SQRT(variance)			       *
*C *						       *
*C * * * * * * * * * * * * * * * * * * * * * * * * * * *
*c	IMPLICIT REAL*8 (A-H,O-Z)
*      INTEGER*4 IU
*      COMMON /RAN/ IU
*
*      PARAMETER (ZERO = 0D0)
*      PARAMETER (SIX  = 6D0)
*      REAL*4 URAND
*
*      IF(S .GT. ZERO)THEN
*	  A=ZERO
*	  DO I=1,12
*	     A=A+URAND(IU)
*	  END DO
*	  EXPVAL=(A-SIX)*S+AM
*      ELSE
*	  EXPVAL=AM
*      END IF
*
*      END

