
       PARAMETER (NWPAWC = 40000)
       PARAMETER (Nch = 15)
        dimension v(nch),ra(nch),su(nch)
       common/pawc/h(nwpawc)

         call hlimit(nwpawc)

         call hbook1(11,' ',100,-2.,2.,0.)
         call hbook1(12,' ',100,-2.,2.,0.)


       open(8,file='allu.dat')

c       sigm=0.3
       nev=100000

       read(8,*)nnn,cutv,step,rat
       sigm=-cutv
       if(nnn.ne.nch)stop 'nnn.ne.nch'

       su(1)=1.
       v(1)=0.
       do i=2,nch
       read(8,*)xs,q2,t,v(i),ra(i)
       su(i)=su(i-1)+step*ra(i)
       enddo

       print *,su       

       do iev=1,nev

       rsu=rndm(-1)*su(nch)

       

       do i=1,nch
         if(rsu.lt.su(i))goto 10
       enddo

10     continue

       call rannor(agau,bgau)

       vgen=v(i)       
       vsme=vgen+sigm*agau
       v000=     sigm*bgau
 
c       vsme=expval(vgen,sigm)
c       v000=expval(0.,sigm)
c       print *,vgen,vsme

          call hfill(11,vsme,0.,1.) 
          call hfill(12,v000,0.,1.) 
        enddo

        call hrput(0,'evbyev.rz','n')  

        end


* * * * * * * * * * * * * * * * * * * * * * * * * * * 
*                                                   * 
*
      FUNCTION URAND(IY)
*                                                   * 
C *   This is a standard pseudo-random generator      * 
C *   that work on IBM-370 and IBM-PC. We don't       * 
C *   know does it work on SUN?                       * 
C *                                                   *
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
      INTEGER*4 IY,M2,IA,IC
      DATA S,M2,IA,IC/.46566128E-9,1073741824,843314861,453816693/
      IY=IY*IA+IC
      IF(IY.LT.0)IY=(IY+M2)+M2
      URAND=FLOAT(IY)*S
      END


*C * * * * * * * * * * * * * * * * * * * * * * * * * * *
*C *                                                   *
      FUNCTION EXPVAL(AM,S)
*C *                                                   *
*C *   This is a simple generator of the pseudo        *
*C *   random numbers with Normal distribution.        *
*C *   AM - mean                                       *
*C *   S  - SQRT(variance)                             *
*C *                                                   *
*C * * * * * * * * * * * * * * * * * * * * * * * * * * *
*c      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 IU
      COMMON /RAN/ IU

      PARAMETER (ZERO = 0D0)
      PARAMETER (SIX  = 6D0)
      REAL*4 URAND

      IF(S .GT. ZERO)THEN
         A=ZERO
         DO I=1,12
            A=A+URAND(IU)
         END DO
         EXPVAL=(A-SIX)*S+AM
      ELSE
         EXPVAL=AM
      END IF

      END



