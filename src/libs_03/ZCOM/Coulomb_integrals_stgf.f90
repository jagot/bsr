!===============================================================================
      Real(8) Function Fdip(ek1,l1,ek2,l2,ifail,tzed)
!===============================================================================
!     Alan Burgess:
!     Calculates the function I(kappa1,l1,kappa2,l2,1) defined in
!     Phil. Trans. Roy. Soc. A226, 255, (1970),
!     where ek1=kappa1**2 and ek2=kappa2**2.
!     It is suitable for use in equations (8),(9),(10) or (11) of
!     J. Phys. B. 7, L364 (1974).
!
!     Calls:  Fdip0   -  for neutral case, tzed = 0.d0
!             Fdip1, Fdip2, Fdipa
!-------------------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: l1, l2
      Real(8), intent(in) :: ek1, ek2, tzed
      Integer :: iprint, ifail
      Real(8) :: emin , emax , t , fa , fd1 , fd2 , fd , rat
      Real(8) :: eps=1.d-4
      Real(8), external :: Fdip0 , Fdip1 , Fdip2 , Fdipa

      if (tzed.eq.0.d0) then
        Fdip=Fdip0(ek1,l1,ek2,l2,eps,ifail)
        Return
      end if

      Fdip=0.d0
      iprint=ifail

      if (ek1+ek2.le.1.d-40) then   ! Arguments too small
        Fdip=0.d0
        ifail=1
        if (iprint.gt.0) Write(iprint,100) ifail
        Return
      else
        if      (ek1-ek2.le.0.D0) then
          emin=ek1
          emax=ek2
        else if (ek1-ek2.gt.0.D0) then
          emin=ek2
          emax=ek1
        endif
        t=emin/emax

! Compare t

        if      (t.le.0.02944D0) then
          fd1=Fdip1(ek1,l1,ek2,l2)
          if (fd1*fd1.le.1.D-40) then
            ifail=2
            if (iprint.gt.0) Write(6,100) ifail
            Fdip=fd1
            Return
          endif
          fd=fd1

        else if (t.gt.0.02944D0) then
          if      (t.ge.0.16667D0) then
            fd2=Fdip2(ek1,l1,ek2,l2)
            if (fd2*fd2.le.1.D-40) then
              ifail=2
              if (iprint.gt.0) Write(*,100) ifail
              Fdip=fd2
              Return
            endif
            fd=fd2
          else if (t.lt.0.16667D0) then
            fd1=Fdip1(ek1,l1,ek2,l2)
            fd=fd1
            if (fd1*fd1.le.1.D-40) then
              fd2=Fdip2(ek1,l1,ek2,l2)
              if (fd2*fd2.le.1.D-40) then
                ifail=2
                if (iprint.gt.0) write(*,100) ifail
                Fdip=fd2
                Return
              endif
              fd=fd2
            endif
          endif
        endif

! ... get Fdip

        if (fd.lt.0.D0 .or. fd.gt.1.D0) then
          ifail=3
          if (iprint.gt.0) write(iprint,100) ifail
          Fdip=0.D0
          Return
        endif

        fa=Fdipa(ek1,l1,ek2,l2)

        ifail=0
        if (fa.eq.0.D0) then
          fa=Fdip0(ek1,l1,ek2,l2,eps,ifail)
          ifail=-ifail
          if (fa.eq.0.D0) then
            Fdip=fd
            Return
          endif
        endif
        rat=fd/fa
        if (rat.gt.10.D0) then
          ifail=4
          if (iprint.gt.0) write(iprint,100) ifail
          Fdip=0.D0
        endif
        Fdip=fd

      endif

  100 Format ("***FDIP FAILURE: IFAIL=",I2)

      End Function Fdip

!===============================================================================
      Real(8) Function Fdip0(ek1,l1,ek2,l2,eps,ifail)
!===============================================================================
! Alan Burgess:
! Calculates the function I0(k1,l1,k2,l2,1) defined in
! Phil. Trans. Roy. Soc. A266, 255 (1970),
! where ek1=k1*k1, ek2=k2*k2, and the relative accuracy is approximately eps.
! It is suitable for use in equations (13) etc. of
! J. Phys. B. 7, L364 (1974)
!--------------------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: l1, l2
      Real(8), intent(in) :: ek1, ek2, eps
      Integer :: l , i , i0 , iprint, ifail
      Real(8) :: e , el , p , p1 , h0 , t , t1 , ti , f , fd , &
                 x , h , s , a , b , c , d , ai , aii , c1
      Real(8), external :: F21

      iprint=ifail
      ifail=0

      If      (l1.eq.l2) Then
        ifail=1
        Fdip0=0.d0
        Return
      Else If (l1.lt.l2) Then
        l=l1
      Else If (l1.gt.l2) Then
        l=l2
      End If

      el=Dfloat(l)
      fd=0.5d0/(el+1.d0)
      If      (ek1.eq.ek2) Then
        Fdip0=fd
        Return
      Else If (ek1.lt.ek2) Then
        e=ek1/ek2
        p=Dfloat(l1-l)
      Else If (ek1.gt.ek2) Then
        e=ek2/ek1
        p=Dfloat(l2-l)
      End If
      fd = fd * e**((el+p+0.5d0)/2.d0)

! To obtain the function ek1 of
! M. J. Seaton, Proc. Phys. Soc. A68, 457 (1955),
! remove the "c" on the next line.

      If      (e.lt.0.5d0) Then      !  was bug: -0.5

        a=el+1.d0
        b=p-0.5d0
        c=el+p+1.5d0
        f=F21(a,b,c,e,eps,ifail)
        l=l+1
        el=Dfloat(l)
        If      (p.le.0.5d0) Then
          c1=el+el+1.d0
        Else If (p.gt.0.5d0) Then
          c1=1.d0
        End If
        Do i = 1 , l
          ai=Dfloat(i)
          aii=ai+ai
          c1=c1*ai*ai*4.d0/(aii*(aii+1.d0))
        End Do
        fd=fd*f*c1

      Else

        p1=p-0.5d0
        t=p1*(el+1.d0)*(e-1.d0)
        i0=l+1
        h0=0.d0
        Do i = 1 , i0
          ti=Dfloat(i)
          h0=h0+1.d0/ti
        End Do
        x=1.d0-e
        h=1.d0-(p+p+h0+Dlog(x/4.d0))
        s=1.d0+t*h
        a=el+1.d0
        b=p1
        c=1.d0
        d=0.d0
        Do While (Dabs(t1).ge.eps*Dabs(s) .and. c.lt.300.d0)
          a=a+1.d0
          b=b+1.d0
          c=c+1.d0
          d=d+1.d0
          t=t*a*b*x/(c*d)
          h=h+p1/(d*b)+el/(c*a)
          t1=t*h
          s=s+t1
        End Do
        If (c.ge.300) Then
          If (iprint.gt.0) Write(6,'(a)') " FAILED TO CONVERGE IN FDIP0"
          ifail=2
        EndIf
        fd=fd*s

      End If

       FDIP0 = fd

      End Function Fdip0

!======================================================================
      Real(8) Function Fdip1 (ek1,l1,ek2,l2)
!======================================================================
!     Dipole Coulomb integral from monopole Coulomb integrals,
!     Burgess et al (1970),  Equation A3:
!
!     (l+1) I(k1,l,k2,l+1;1) = [1+(l+1)^2 k2^2]^1/2 I(k1,l,k2,l;0) -
!                              [1+(l+1)^2 k1^2]^1/2 I(k1,l+1,k2,l+1;0)
!
!     On entry:  ek1 = k1^2;  ek2 = k2^2
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: l1,l2
      Real(8), intent(in) :: ek1,ek2
      Integer :: l,lp
      Real(8) :: a1,a2,b1,b2
      Real(8), external :: Fmon1   ! I(k1,l,k2,l,0)

      if(abs(l1-l2).ne.1) then
       write(*,*) 'WARNING: |l1-l2| \= 1  in Fdip1'
       Fdip1 = 0.d0
       Return
      end if

      if (l1.lt.l2) then
        a1=ek1;  a2=ek2;   l=l1
      else
        a1=ek2;  a2=ek1;   l=l2
      end If

      lp = l+1
      b1 = sqrt(1.d0 + lp*lp*a2)*Fmon1(ek1,ek2,l )
      b2 = sqrt(1.d0 + lp*lp*a1)*Fmon1(ek1,ek2,lp)
      Fdip1=(b1-b2)/lp

      End Function Fdip1


!=====================================================================
      Real(8) Function F21 (a,b,c,d,eps,ifail)
!=====================================================================
!     generalized hypergeometric function:
!     ...
!---------------------------------------------------------------------
      Implicit none
      Real(8), intent(in) :: a,b,c,d,eps
      Integer :: ifail, i,iprint
      Real(8) :: t,dd,tn1,tn2,sum,ai,at,as

      iprint=ifail
      ifail =0

      t   = (a*b*d)/c
      sum = 1.d0 + t

      dd  = 1.d0/(1.d0-d)
      tn1 = 0.d0
      i   = 1
      at  = 2.d0
      as  = 1.d0
      Do while (at.gt.as .and. i.le.300)   !  300 ? in arguments?
        ai  = dfloat(i)
        t   = t*(a+ai)*(b+ai)*d/((c+ai)*(1.d0+ai))
        tn2 = t*dd
        F21 = sum + tn2
        sum = sum + t
        at  = abs(t+tn2-tn1)
        as  = abs(F21)*eps
        tn1 = tn2
        i   = i + 1
      End Do

      if (i.gt.300) then
        If (iprint.gt.0) Write(iprint,'(a)') " FAILED TO CONVERGE IN F21"
        ifail=3
      end if

      End Function F21


!===============================================================
       Real(8) Function Argam(l,a)
!===============================================================
! Calculates argument of Gamma(l+1+i*a)
! where l is an integer not less than zero
!===============================================================
       Implicit none
       Integer, intent(in) :: l
       Real(8), intent(in) :: a
       Integer :: j , j0 , j1
       Real(8) :: b , c , d , z , g0 , g1 , g2 , g3 , u

       b=Dabs(a)
       b=250.D0*b**0.25d0-a*a

!       b = 250.D0 * Dnroot(4,Dabs(a)) - a*a

       j0=l+1
       c=Dfloat(j0)
       d=c*c
       z=0.D0

       If (d.lt.b) Then
         b = sqrt(b)
         j1 = int(b)
         Do j = j0 , j1
           d = a/Dfloat(j)
           If      (abs(d).lt.0.1D0) Then
             g1=d*d
             g2=-35.D0*g1+45.D0
             g2=-g1*g2+63.D0
             g2=-g1*g2+105.D0
             g1=d-d*g1*g2/315.D0
           Else
             g1=Datan(d)
           End If
           z=z+g1
         End Do
         j0=j1+1
       End If

       d=Dfloat(j0)
       g0=d*d
       u=a*a
       g1=1.D0/(g0+u)
       g2=g1*g1
       g3=10.D0*g0*g0-20.D0*g0*u+2.D0*u*u
       g3=g3*g2-21.D0*g0+7.D0*u
       g3=g3*g2+210.D0
       g1=a*g3*g1/2520.D0

       Argam=-z+0.5D0*a*log(g0+u)+(d-0.5D0)*atan(a/d)-a-g1

       End Function Argam

!======================================================================
      Real(8) Function Fmon1(ek1,ek2,l)
!======================================================================
!    Calculates I(ek1,ek2,l;0) as described in Burgess 1970:  Eq. A14
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: l
      Real(8), intent(in) :: ek1, ek2
      Integer :: m, ma1, ma2, mg
      Real(8) :: pi, vmax, x1, x2, x3, x4, x5, x6, x7, t, a1, a2, b, g, &
                 eta, em, emm, u, v, w, s, u0, u1, t0, t1, tm, s0, s1, sm

      If (ek1+ek2.le.1.d-40) Then
        Fmon1=1.D50                        ! ???
        Return
      End If

      vmax=2000000.d0
      x1=sqrt(ek1)
      x2=sqrt(ek2)
      x3=x1+x2
      x4=x3*x3
      x5=x1*x2
      x6=x2-x1
      x7=4.d0/x4
      pi=3.141592653589793d0
      If (ek1-ek2.le.0.d0) Then
        eta=1.d0/x2
      Else
        eta=1.d0/x1
      End If
      g = 0.5d0*pi*exp(-pi*eta)

      If (g.eq.0.d0) Then         ! NRB            ???
        Fmon1=0.d0
        Return
      End If

      a1=1.d0
      a2=1.d0
      mg=0
      ma1=0
      ma2=0
      m=-1

      Do While (m-l.lt.0)
        m=m+1
        em=Dfloat(m)
        t=em+em+1.d0
        g=g*x7/(t*(t+1.d0))
        emm=em*em
        a1=a1*(1.d0+emm*ek1)
        a2=a2*(1.d0+emm*ek2)
        Do While (g.lt.0.015625d0)
          g=64.d0*g
          mg=mg-1
        End Do
        Do While (g.gt.64.d0)
          g=0.015625d0*g
          mg=mg+1
        End Do
        Do While (a1.gt.64.d0)
          a1=0.015625d0*a1
          ma1=ma1+1
        End Do
        Do While (a2.gt.64.d0)
          a2=0.015625d0*a2
          ma2=ma2+1
        End Do
      End Do

      g=g*(t+1.d0)

      If (x1.ge.300.d0) Then
        b=pi/x1
        a1=1.5d0*a1/(b*(3.d0-b*(3.d0-b*(2.d0-b))))
      Else
        If (x1.gt.0.2d0) Then
          b=-pi/x1
          a1=a1/(1.d0-Dexp(2.d0*b))
        End If
      End If

      If (x2.ge.300.d0) Then
        b=pi/x2
        a2=1.5d0*a2/(b*(3.d0-b*(3.d0-b*(2.d0-b))))
      Else
        If (x2.gt.0.2d0) Then
          b=-pi/x2
          a2=a2/(1.d0-Dexp(2.d0*b))
        End If
      End If

      g=g*sqrt(a1*a2)*(8.0d0)**(mg+mg+ma1+ma2)

      s0=1.d0
      s1=0.d0
      u=Dfloat(l)
      v=0.d0
      w=u+u+1.d0
      t0=1.d0
      t1=0.d0
      s=s0*s0+s1*s1
      t=t0*t0+t1*t1

      Do While (s.lt.1.D24*t)
        u=u+1.d0
        v=v+1.d0
        w=w+1.d0
        If (v.gt.vmax) Then
          Fmon1=0.d0
          Return
        End If
        u0=u*u*x5+1.d0
        u1=u*x6
        t=t0*u0-t1*u1
        t1=t0*u1+t1*u0
        t0=t
        t=x7/(v*w)
        t0=t*t0
        t1=t*t1
        s0=s0+t0
        s1=s1+t1
        s=s0*s0+s1*s1
        t=t0*t0+t1*t1
        sm=1.d0/s
        tm=1.d0/t
        If (sm*tm.eq.0.d0) Then      ! NRB
          Fmon1=0.d0
          Return
        End If
      End Do

      Fmon1=g*Dsqrt(s)

      End Function Fmon1

!======================================================================
      Real(8) Function Fdipa(ek1,l1,ek2,l2)
!======================================================================
!     Asymptotic expression for I(k1,l1,k2,l2,1) based on A40,1 of BHT
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: l1, l2
      Real(8), intent(in) :: ek1, ek2
      Integer :: l, icase
      Real(8) :: pi , x1 , x2 , xp , e , ee , t , t0 ,  tl , tl2 , f0

! ... Arguments too small

      if (ek1*ek2.lt.1.D-50) icase=0
      x1 = 1.d0/sqrt(ek1)
      x2 = 1.d0/sqrt(ek2)
      xp = abs(x1-x2)

! ... Arguments too large

      if (xp.gt.100.d0) icase=0

! ... non zero

      pi = acos(-1.0d0)
      xp = exp(0.5d0*pi*xp)

      if      (ek1-ek2.le.0) then
        e=ek1/ek2
        if      (l1-l2.le.0) then
          l=l1;  icase=1
        Else if (l1-l2.gt.0) then
          l=l2;  icase=2
        End if
      else
        e=ek2/ek1
        if      (l1-l2.le.0) then
          l=l1;  icase=2
        Else if (l1-l2.gt.0) then
          l=l2;  icase=1
        End if
      End if

! ... Calculate for different cases

      if      (icase.eq.0) then

! ... Zero
        Fdipa=0.d0

      Else if (icase.eq.1) then

! .. formula A40:

        tl=Dfloat(l)
        t0=1.d0-e
        if (tl*t0.lt.e) then
          Fdipa=0.d0
        Else
          t=pi*tl
          ee=Dsqrt(e)
          f0=Dsqrt(t*t0*ee)*ee**l
          tl=Dfloat(l+l+1)
          Fdipa=f0*xp/tl
        End if

      Else if (icase.eq.2) then
! .. formular A41:
        t0=1.d0-e
        tl=l
        if (tl*t0.lt.e) then
          Fdipa=0.d0
        Else
          t0=1.d0/t0
          t=tl*pi
          ee=Dsqrt(e)
          f0=Dsqrt(t*t0*ee)*ee**(l+1)
          tl=l+l+1
          tl2=l+l+3
          Fdipa=f0*xp/(tl*tl2)
        End if
      End if

      End Function Fdipa

!===========================================================================
       Real(8) Function Fdip2(ek1,l1,ek2,l2)
!===========================================================================
! Calculates I(ek1,l1,ek2,l2;1) when ek1 is close to ek2.
! Eqn A30/A31
!---------------------------------------------------------------------------
       Implicit none
       Integer l1, l2
       Real(8) ek1, ek2
       Integer l, j1, j2
       Real(8) fd, pi, wmax, eta1, eta2, a, b, c, c1, c2, w0, w1, t, t0, t1, &
               u0, u1, v0, v1, x, x0, y0, y1, y2, z0, z1, p, p0, p1, q0, q1
       Real(8), external :: argam          !Added because math.mod is not working

      wmax=200.D0                                         !A30/A31
      eta1=1.D0/Dsqrt(ek1)
      eta2=1.D0/Dsqrt(ek2)
      w1=eta2-eta1
      pi=3.141592653589793d0
      a=abs(w1)
      b=pi*a
      If     (b.le.1.D-2) Then

        c=Dsqrt(3.D0/(3.D0-b*(3.D0-b*(2.D0-b))))
      Else If (b.gt.1.D-2) Then
        If      (b.lt.14.D0) Then
          b=b+b
          c1=1.D0-Dexp(-b)
          c=Dsqrt(b/c1)
        Else If (b.ge.14.D0) Then
          c=Dsqrt(b+b)
        End If
      End If
      c=0.5D0*c/Dsqrt(eta1*eta2)
      c2=eta1+eta2
      c1=4.D0*eta1*eta2/(c2*c2)
      If      (l2.gt.l1) Then
        l=l1
      Else If (l2.le.l1) Then
        l=l2
        t1=eta1
        eta1=eta2
        eta2=t1
        w1=-w1
      End If
      c=c*c1**(l+1)
      u0=Dfloat(l+1)
      u1=eta1
      v0=u0
      v1=-eta2
      w0=1.D0
      x0=w1/(c2*c2)
      y2=-eta2-eta2
      y0=-u0*w1+y2
      y1=eta2*w1
      t1=x0/(1.D0+w1*w1)
      z0=u0*t1
      z1=u1*t1
      t=z0-z1*w1
      z1=z0*w1+z1
      z0=t
      q0=-1.D0+z0*y0-z1*y1
      q1=z0*y1+z1*y0
      x=w1*x0
      t0=t1

      Do While (abs(t0).le.1.D24*abs(t1))   ! was bug for negative t0

        u0=u0+1.0d0
        v0=v0+1.0d0
        w0=w0+1.0d0

        If (w0.gt.wmax) Then
          Fdip2=0.d0;  Return
        End If

        y0=y0+y2
        t=z0*u0-z1*u1
        z1=z0*u1+z1*u0
        z0=t
        t=z0*v0-z1*v1
        z1=z0*v1+z1*v0
        z0=t
        t=z0*w0-z1*w1
        z1=z0*w1+z1*w0
        z0=t
        x0=x/(w0*(w0*w0+w1*w1))
        z0=z0*x0
        z1=z1*x0
        t0=z0*y0-z1*y1
        t1=z0*y1+z1*y0
        q0=q0+t0
        q1=q1+t1

        t1=t0*t0+t1*t1
        t0=q0*q0+q1*q1
      End Do

      j1=0
      j2=l+1
      p=Argam(j1,w1)+Argam(l,eta1)-Argam(j2,eta2)

      If (a.gt.1.d-40) p = p + w1*log(c2/a)
      p0=cos(p)
      p1=sin(p)
      t=p0*q0-p1*q1
      q1=p0*q1+p1*q0

      Fdip2=c*q1

      Return                                              !Missing rturn statement?

      End Function Fdip2
