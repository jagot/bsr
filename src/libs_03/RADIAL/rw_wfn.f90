!====================================================================
      Subroutine RW_w(nuw,m,job)
!====================================================================
!
!     read radial function 'm' from w-file
!
!--------------------------------------------------------------------

      USE radial

      IMPLICIT NONE
      Integer(4) :: nuw,m
      Character(1) :: job
      Real(8) :: ZZ,EI,ZETA,ZR
      Character(3) :: EL3
      Character(3), External :: ELF3
      Character(4), External :: ELF4

      if(job.eq.'w') then
       EL3 = ELF3(nro(m),lro(m),kro(m))
       WRITE(nuw) Atom,Term,EL3,mro(m),Z,0.d0,0.d0,AZ(m)
       WRITE(nuw) P(1:mro(m),m)
       Return
      end if

      if(m.ge.mrf-2) CALL Alloc_radial(mrf+jrf)

      READ(nuw,end=2) Atom,Term,EL3,mro(m),ZZ,EI,ZETA,AZ(m)
      READ(nuw) P(1:mro(m),m)

      if(ZZ.ne.Z) Stop ' ZZ <> Z  in frm-fail'

      Call EL3_nlk(EL3,nro(m),lro(m),kro(m))

      ero(m) = ELF4(nro(m),lro(m),kro(m))

      if(mro(m).lt.NR) P(mro(m)+1:NR,m) = 0.d0
      mexp(m) = lro(m) + 1
      aexp(m) = - D1/mexp(m)
      ZR = Z*R(1)
      bexp(m) = P(1,m)/(AZ(m)*R2(1)*R(1)**lro(m))
      bexp(m) = (bexp(m) - D1 + ZR/(lro(m)+1) )/ZR**2
      Return

    2 m = -1

      End Subroutine RW_w




!====================================================================
      Subroutine RW_frm(nuw,m,job)
!====================================================================
!
!     read radial function 'm' from rfm-file
!
!--------------------------------------------------------------------

      USE radial

      IMPLICIT NONE
      Integer(4) :: nuw,m
      Character(1) :: job
      Real(8) :: ZZ,EI,ZETA,ZR

      if(job.eq.'w') then
       WRITE(nuw,1) Atom,Term,ero(m),mro(m),Z,0.d0,0.d0,AZ(m)
       WRITE(nuw,2) P(1:mro(m),m)
       Return
      end if

      if(m.ge.mrf-2) CALL Alloc_radial(mrf+jrf)

    1 Format(A6,A6,2X,A4,I6,F6.2,3(E18.10))
    2 Format(4(E18.10))

      READ(nuw,1,end=10) Atom,Term,ero(m),mro(m),ZZ,EI,ZETA,AZ(m)
      READ(nuw,2) P(1:mro(m),m)

      if(ZZ.ne.Z) Stop ' ZZ <> Z  in frm-fail'

      Call EL4_nlk(ero(m),nro(m),lro(m),kro(m))

      if(mro(m).lt.NR) P(mro(m)+1:NR,m) = 0.d0
      mexp(m) = lro(m) + 1
      aexp(m) = - D1/mexp(m)
      ZR = Z*R(1)
      bexp(m) = P(1,m)/(AZ(m)*R2(1)*R(1)**lro(m))
      bexp(m) = (bexp(m) - D1 + ZR/(lro(m)+1) )/ZR**2
      Return

   10 m = -1


      End Subroutine RW_frm


!====================================================================
      Subroutine Read_frm(nuw)
!====================================================================
!     read radial functions from frm-file (new list)
!--------------------------------------------------------------------

      USE radial

      IMPLICIT NONE
      Integer(4) :: nuw,m
      Real(8) :: ZZ,EI,ZETA,ZR,AZZ
      Character(3) :: EL

      CALL Alloc_radial(0)
      read(nuw,'(24x,F6.2)') Z; rewind(nuw);  Call INITR
      CALL Alloc_radial(irf)

    1 Format(A6,A6,2X,A4,I6,F6.2,3(E18.10))
    2 Format(4(E18.10))

      nrf = 0
    5 Continue

      READ(nuw,1,end=10) Atom,Term,EL,m,ZZ,EI,ZETA,AZZ
      nrf = nrf + 1
      if(nrf.gt.mrf-2) Call Alloc_radial(mrf+jrf) 

      ero(nrf) = EL
      mro(nrf) = m;   m = nrf
      AZ(nrf) = AZZ

      READ(nuw,2) P(1:mro(nrf),nrf)

      Call EL4_nlk(ero(nrf),nro(nrf),lro(nrf),kro(nrf))

      if(mro(nrf).lt.NR) P(mro(nrf)+1:NR,nrf) = 0.d0
      mexp(m) = lro(m) + 1
      aexp(m) = - D1/mexp(m)
      ZR = Z*R(1)
      bexp(m) = P(1,m)/(AZ(m)*R2(1)*R(1)**lro(m))
      bexp(m) = (bexp(m) - D1 + ZR/(lro(m)+1) )/ZR**2
      go to 5

   10 m = -1


      End Subroutine Read_frm

