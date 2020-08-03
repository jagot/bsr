!   Bloch conditions at r=a
!   ***********************
!=======================================================================
      Subroutine Bloch0(basis,kappa,n,nu,c,nhm,hm)
!=======================================================================
!
!   method 0:
!                      |       0      nu      |                       
!               L = c  |                      |
!                      |  (nu-1)      0       |                       
!    
!   simple Bloch operator that symmetrize the matrix without forcing 
!   any Q(a)/P(a) relation
!-----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss

      Implicit none
      Character(2) :: basis
      Integer :: kappa, nhm
      Real(8) :: n, nu, c, dn,bn, bp,bm, dp,dm, a11,a21,a12,a22,&
                 hm(nhm,nhm)

      ! parameter n is not used
 
      a21 = (nu-1.d0)   * c 
      a12 =  nu         * c 

      if(basis.eq.'bb') then

       hm(ns,ms) = hm(ns,ms) + a12
       hm(ms,ns) = hm(ms,ns) + a21

      elseif(basis.eq.'bm') then

       hm(nsp,ns+nsq) = hm(nsp,ns+nsq) + a12
       hm(ns+nsq,nsp) = hm(ns+nsq,nsp) + a21

      elseif(basis.eq.'dp') then

       dn=dpbs(nv+1,1,ks  )
       bn=dpbs(nv+1,1,ks-1)

       hm(ns  ,ms-1) = hm(ns  ,ms-1) + a12 * bn
       hm(ns  ,ms  ) = hm(ns  ,ms  ) + a12 * dn
       hm(ms-1,ns  ) = hm(ms-1,ns  ) + a21 * bn
       hm(ms  ,ns  ) = hm(ms  ,ns  ) + a21 * dn

       dn=dpbs(nv+1,2,1)
       bn=dpbs(nv+1,2,2)

       hm(   1,ns+2) = hm(   1,ns+2) - a12 * bn
       hm(   1,ns+1) = hm(   1,ns+1) - a12 * dn
       hm(ns+2,   1) = hm(ns+2,   1) - a21 * bn
       hm(ns+1,   1) = hm(ns+1,   1) - a21 * dn

      elseif(basis.eq.'dm') then

       dn=dmbs(nv+1,1,ks  )
       bn=dmbs(nv+1,1,ks-1)

       hm(ns-1,ms  ) = hm(ns-1,ms  ) + a12 * bn
       hm(ns  ,ms  ) = hm(ns  ,ms  ) + a12 * dn
       hm(ms  ,ns-1) = hm(ms  ,ns-1) + a21 * bn
       hm(ms  ,ns  ) = hm(ms  ,ns  ) + a21 * dn

      elseif(basis.eq.'dd') then

       bp =  dpbs(nv+1,1,ks-1) 
       bm =  dmbs(nv+1,1,ks-1) 
       dp =  dpbs(nv+1,1,ks)   
       dm =  dmbs(nv+1,1,ks)   

       hm(ns  ,ns-1)=hm(ns  ,ns-1)+ a12*bp          
       hm(ns-1,ns  )=hm(ns-1,ns  )+          a21*bp   
       hm(ns  ,ns  )=hm(ns  ,ns  )+ a12*dp + a21*dp  

       hm(ns-1,ms-1)=hm(ns-1,ms-1)+          a21*bp*bm           
       hm(ns  ,ms-1)=hm(ns  ,ms-1)+          a21*dp*bm  
       hm(ns-1,ms  )=hm(ns-1,ms  )+          a21*bp*dm  ! ???  
       hm(ns  ,ms  )=hm(ns  ,ms  )+ a12    + a21*dp*dm   

       hm(ms-1,ns-1)=hm(ms-1,ns-1)+          a12*bp*bm           
       hm(ms  ,ns-1)=hm(ms  ,ns-1)+          a12*bp*dm  
       hm(ms-1,ns  )=hm(ms-1,ns  )+          a12*dp*bm  ! ??? 
       hm(ms  ,ns  )=hm(ms  ,ns  )+ a21    + a12*dp*dm   
      
       hm(ms  ,ms-1)=hm(ms  ,ms-1)+ a21*bm         
       hm(ms-1,ms  )=hm(ms-1,ms  )+          a12*bm   
       hm(ms  ,ms  )=hm(ms  ,ms  )+ a21*dm + a12*dm   

      else

!       Stop 'bloch0: unknown basis'

      end if

      END SUBROUTINE Bloch0


!=======================================================================
      SUBROUTINE Bloch1(basis,kappa,n,nu,c,nhm,hm)
!=======================================================================
!
!   method 1:
!                      | -nu * n      nu      |                       Q(a)
!               L = c  |                      |,  n=(b+kappa)/2*a*c = --- 
!                      |  (nu-1)    (1-nu)/n  |                       P(a)
!
!    This formular from Szmytkowski and Hinze, J.Phys.B,29,761,1996
!    For nu=1/2, it agrees with expression from I.Grant (2007)
!    
!-----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss

      Implicit none
      Character(2) :: basis
      Integer(4) :: kappa,nhm, i1,i2
      Real(8) :: n, nu, c, dn,bn, bp,bm, dp,dm, a11,a21,a12,a22,&
                 hm(nhm,nhm)
 
      a11 = -nu*n       * c 
      a21 = (nu-1.d0)   * c 
      a12 =  nu         * c 
      a22 =-(nu-1.d0)/n * c 

      if(basis.eq.'bb') then

       i1 = ns; i2= ms
       hm(i1,i1) = hm(i1,i1) + a11 
       hm(i1,i2) = hm(i1,i2) + a12  
       hm(i2,i1) = hm(i2,i1) + a21
       hm(i2,i2) = hm(i2,i2) + a22 
  
      elseif(basis.eq.'bm') then

       i1 = nsp; i2= ns+nsq
       hm(i1,i1) = hm(i1,i1) + a11 
       hm(i1,i2) = hm(i1,i2) + a12  
       hm(i2,i1) = hm(i2,i1) + a21
       hm(i2,i2) = hm(i2,i2) + a22 
      
      elseif(basis.eq.'dp') then

       dn=dpbs(nv+1,1,ks  )
       bn=dpbs(nv+1,1,ks-1)

       hm(ns  ,ns  ) = hm(ns  ,ns  ) + a11 
       hm(ns  ,ms-1) = hm(ns  ,ms-1) + a12 * bn
       hm(ns  ,ms  ) = hm(ns  ,ms  ) + a12 * dn
       hm(ms-1,ns  ) = hm(ms-1,ns  ) + a21 * bn
       hm(ms  ,ns  ) = hm(ms  ,ns  ) + a21 * dn
       hm(ms-1,ms-1) = hm(ms-1,ms-1) + a22 * bn * bn 
       hm(ms-1,ms  ) = hm(ms-1,ms  ) + a22 * bn * dn 
       hm(ms  ,ms-1) = hm(ms  ,ms-1) + a22 * dn * bn 
       hm(ms  ,ms  ) = hm(ms  ,ms  ) + a22 * dn * dn 

      elseif(basis.eq.'dm') then

       dn=dmbs(nv+1,1,ks  )
       bn=dmbs(nv+1,1,ks-1)

       hm(ns-1,ns-1) = hm(ns-1,ns-1) + a11 * bn * bn 
       hm(ns-1,ns  ) = hm(ns-1,ns  ) + a11 * bn * dn 
       hm(ns  ,ns-1) = hm(ns  ,ns-1) + a11 * dn * bn 
       hm(ns  ,ns  ) = hm(ns  ,ns  ) + a11 * dn * dn 
       hm(ns-1,ms  ) = hm(ns-1,ms  ) + a12 * bn
       hm(ns  ,ms  ) = hm(ns  ,ms  ) + a12 * dn
       hm(ms  ,ns-1) = hm(ms  ,ns-1) + a21 * bn
       hm(ms  ,ns  ) = hm(ms  ,ns  ) + a21 * dn
       hm(ms  ,ms  ) = hm(ms  ,ms  ) + a22 

      elseif(basis.eq.'dd') then

       bp =  dpbs(nv+1,1,ks-1) 
       bm =  dmbs(nv+1,1,ks-1) 
       dp =  dpbs(nv+1,1,ks)   
       dm =  dmbs(nv+1,1,ks)   

       hm(ns-1,ns-1)=hm(ns-1,ns-1)+                       + a22*bp*bp           
       hm(ns  ,ns-1)=hm(ns  ,ns-1)+       a12*bp          + a22*dp*bp  
       hm(ns-1,ns  )=hm(ns-1,ns  )+                a21*bp + a22*bp*dp   
       hm(ns  ,ns  )=hm(ns  ,ns  )+ a11 + a12*dp + a21*dp + a22*dp*dp   

       hm(ns-1,ms-1)=hm(ns-1,ms-1)+                       + a21*bp*bm           
       hm(ns  ,ms-1)=hm(ns  ,ms-1)+       a11*bm          + a21*dp*bm  
       hm(ns-1,ms  )=hm(ns-1,ms  )+                a22*bp + a21*bp*dm  ! ???  
       hm(ns  ,ms  )=hm(ns  ,ms  )+ a12 + a11*dm + a22*dp + a21*dp*dm   

       hm(ms-1,ns-1)=hm(ms-1,ns-1)+                       + a12*bp*bm           
       hm(ms  ,ns-1)=hm(ms  ,ns-1)+       a22*bp          + a12*bp*dm  
       hm(ms-1,ns  )=hm(ms-1,ns  )+                a11*bm + a12*dp*bm  ! ??? 
       hm(ms  ,ns  )=hm(ms  ,ns  )+ a21 + a22*dp + a11*dm + a12*dp*dm   
      
       hm(ms-1,ms-1)=hm(ms-1,ms-1)+                       + a11*bm*bm           
       hm(ms  ,ms-1)=hm(ms  ,ms-1)+       a21*bm          + a11*dm*bm  
       hm(ms-1,ms  )=hm(ms-1,ms  )+                a12*bm + a11*bm*dm   
       hm(ms  ,ms  )=hm(ms  ,ms  )+ a22 + a21*dm + a12*dm + a11*dm*dm   

      else

       Stop 'bloch1: unknown basis'

      end if

      END SUBROUTINE Bloch1
      

!=======================================================================
      SUBROUTINE Bloch2(basis,kappa,n,nu,c,nhm,hm,cm,ip)
!=======================================================================
!
!   method 2: 
!                       | 0   nu |
!               L = c   |        |
!                       |nu-1   0|
!
!   and the condition Q(a)/P(a)=n is satisfied by using 
!   combined B-spline elements related to r=a 
!
!-----------------------------------------------------------------------

      Use DBS_grid
      Use DBS_gauss

      Implicit none

      Character(2) :: basis
      Integer(4) :: kappa, nhm,i,j, ip(nhm)
      Real(8) :: n, nu, c, a21,a12, dn,bn, bp,bm, dp,dm, &
                 hm(nhm,nhm),cm(nhm,nhm)

      a21 = (nu-1.d0) * c 
      a12 =  nu       * c 

      if(basis.eq.'bb') then

       hm(ns,ms) = hm(ns,ms) + a12
       hm(ms,ns) = hm(ms,ns) + a21
      
       hm(ns,:) = hm(ns,:) + hm(ms,:) * n
       cm(ns,:) = cm(ns,:) + cm(ms,:) * n

       hm(:,ns) = hm(:,ns) + hm(:,ms) * n
       cm(:,ns) = cm(:,ns) + cm(:,ms) * n

       ip(ms)=0

      elseif(basis.eq.'bm') then

       hm(ns,nhm) = hm(ns,nhm) + a12
       hm(nhm,ns) = hm(nhm,ns) + a21
      
       hm(ns,:) = hm(ns,:) + hm(nhm,:) * n
       cm(ns,:) = cm(ns,:) + cm(nhm,:) * n

       hm(:,ns) = hm(:,ns) + hm(:,nhm) * n
       cm(:,ns) = cm(:,ns) + cm(:,nhm) * n

       ip(nhm)=0

      elseif(basis.eq.'dp') then

       dn=dpbs(nv+1,1,ks  )
       bn=dpbs(nv+1,1,ks-1)

       hm(ns  ,ms-1) = hm(ns  ,ms-1) + a12 * bn
       hm(ns  ,ms  ) = hm(ns  ,ms  ) + a12 * dn
       hm(ms-1,ns  ) = hm(ms-1,ns  ) + a21 * bn
       hm(ms  ,ns  ) = hm(ms  ,ns  ) + a21 * dn

       dn = dn/n; bn=bn/n 

       i=ms;  j=ns

       hm(i,:) = hm(i,:) + hm(j,:) * dn
       cm(i,:) = cm(i,:) + cm(j,:) * dn

       hm(:,i) = hm(:,i) + hm(:,j) * dn
       cm(:,i) = cm(:,i) + cm(:,j) * dn

       i=ms-1; j=ns

       hm(i,:) = hm(i,:) + hm(j,:) * bn
       cm(i,:) = cm(i,:) + cm(j,:) * bn

       hm(:,i) = hm(:,i) + hm(:,j) * bn
       cm(:,i) = cm(:,i) + cm(:,j) * bn

       i=ms;  j=ms-1

       ip(ns)=0

      elseif(basis.eq.'dm') then

       dn=dmbs(nv+1,1,ks  )
       bn=dmbs(nv+1,1,ks-1)

       hm(ns-1,ms  ) = hm(ns-1,ms  ) + a12 * bn
       hm(ns  ,ms  ) = hm(ns  ,ms  ) + a12 * dn
       hm(ms  ,ns-1) = hm(ms  ,ns-1) + a21 * bn
       hm(ms  ,ns  ) = hm(ms  ,ns  ) + a21 * dn

       dn = dn*n; bn=bn*n 

       i=ns;  j=ms
  
       hm(i,:) = hm(i,:) + hm(j,:) * dn
       cm(i,:) = cm(i,:) + cm(j,:) * dn
  
       hm(:,i) = hm(:,i) + hm(:,j) * dn
       cm(:,i) = cm(:,i) + cm(:,j) * dn
  
       i=ns-1; j=ms

       hm(i,:) = hm(i,:) + hm(j,:) * bn
       cm(i,:) = cm(i,:) + cm(j,:) * bn
  
       hm(:,i) = hm(:,i) + hm(:,j) * bn
       cm(:,i) = cm(:,i) + cm(:,j) * bn
  
       i=ns;  j=ms
  
       ip(ms)=0

      elseif(basis.eq.'dd') then        

       Stop 'Bloch2: method 2 does not work for dd basis'

      else

       Stop 'Bloch2: unknown basis'

      end if

      END SUBROUTINE Bloch2


!=======================================================================
      SUBROUTINE Bloch3(basis,kappa,n,nu,c,nhm,hm,cm,ip)
!=======================================================================
!   method 3: 
!
!   Impliment only the condition Q(a)/P(a)=n, which is forcing by using 
!   combined B-spline elements related to r=a 
!-----------------------------------------------------------------------

      Use DBS_grid
      Use DBS_gauss

      Implicit none

      Character(2) :: basis
      Integer(4) :: kappa, nhm, i,j, ip(nhm)
      Real(8) :: n, nu, c, a21,a12, dn,bn, bp,bm, dp,dm, &
                 hm(nhm,nhm),cm(nhm,nhm)

      a21 = (nu-1.d0) * c 
      a12 =  nu       * c 

      if(basis.eq.'bb') then

!       hm(ns,ms) = hm(ns,ms) + a12
!       hm(ms,ns) = hm(ms,ns) + a21
      
       hm(ns,:) = hm(ns,:) + hm(ms,:) * n
       cm(ns,:) = cm(ns,:) + cm(ms,:) * n

       hm(:,ns) = hm(:,ns) + hm(:,ms) * n
       cm(:,ns) = cm(:,ns) + cm(:,ms) * n

       ip(ms)=0

      elseif(basis.eq.'bm') then


       hm(ns,:) = hm(ns,:) + hm(nhm,:) * n
       cm(ns,:) = cm(ns,:) + cm(nhm,:) * n

       hm(:,ns) = hm(:,ns) + hm(:,nhm) * n
       cm(:,ns) = cm(:,ns) + cm(:,nhm) * n

       ip(nhm)=0

      elseif(basis.eq.'dp') then

       dn=dpbs(nv+1,1,ks  )
       bn=dpbs(nv+1,1,ks-1)

!       hm(ns  ,ms-1) = hm(ns  ,ms-1) + a12 * bn
!       hm(ns  ,ms  ) = hm(ns  ,ms  ) + a12 * dn
!       hm(ms-1,ns  ) = hm(ms-1,ns  ) + a21 * bn
!       hm(ms  ,ns  ) = hm(ms  ,ns  ) + a21 * dn

       dn = dn/n; bn=bn/n 

       i=ms;  j=ns

       hm(i,:) = hm(i,:) + hm(j,:) * dn
       cm(i,:) = cm(i,:) + cm(j,:) * dn

       hm(:,i) = hm(:,i) + hm(:,j) * dn
       cm(:,i) = cm(:,i) + cm(:,j) * dn

       i=ms-1; j=ns

       hm(i,:) = hm(i,:) + hm(j,:) * bn
       cm(i,:) = cm(i,:) + cm(j,:) * bn

       hm(:,i) = hm(:,i) + hm(:,j) * bn
       cm(:,i) = cm(:,i) + cm(:,j) * bn

       i=ms;  j=ms-1

       ip(ns)=0

      elseif(basis.eq.'dm') then

       dn=dmbs(nv+1,1,ks  )
       bn=dmbs(nv+1,1,ks-1)

!       hm(ns-1,ms  ) = hm(ns-1,ms  ) + a12 * bn
!       hm(ns  ,ms  ) = hm(ns  ,ms  ) + a12 * dn
!       hm(ms  ,ns-1) = hm(ms  ,ns-1) + a21 * bn
!       hm(ms  ,ns  ) = hm(ms  ,ns  ) + a21 * dn

       dn = dn*n; bn=bn*n 

       i=ns;  j=ms
  
       hm(i,:) = hm(i,:) + hm(j,:) * dn
       cm(i,:) = cm(i,:) + cm(j,:) * dn
  
       hm(:,i) = hm(:,i) + hm(:,j) * dn
       cm(:,i) = cm(:,i) + cm(:,j) * dn
  
       i=ns-1; j=ms

       hm(i,:) = hm(i,:) + hm(j,:) * bn
       cm(i,:) = cm(i,:) + cm(j,:) * bn
  
       hm(:,i) = hm(:,i) + hm(:,j) * bn
       cm(:,i) = cm(:,i) + cm(:,j) * bn
  
       i=ns;  j=ms
  
       ip(ms)=0

      elseif(basis.eq.'dd') then        

       Stop 'Bloch3: method 3 does not work for dd basis'

      else

       Stop 'Bloch3: unknown basis'

      end if

      END SUBROUTINE Bloch3



!=======================================================================
      SUBROUTINE Bloch(basis,n,nhm,hm,cm,ip)
!=======================================================================
!     Impliment only the condition Q(a)= n P(a), 
!     which is forcing by using combined B-spline elements related to r=a 
!-----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss

      Implicit none
      Character(2) :: basis
      Integer(4) :: kappa, nhm, i,j, ip(nhm)
      Real(8) ::  n, dn,bn, hm(nhm,nhm),cm(nhm,nhm)

      if(basis.eq.'bb') then

       Call Replace_line(ms,hm,cm,ip,ns,ms,n)

      elseif(basis.eq.'bm') then

       Call Replace_line(ms,hm,cm,ip,nsp,ns+nsq,n)

      elseif(basis.eq.'bd') then

       ip(ms)=0
       bn = n /bspd(nv+1,1,ks-1) 
       Call Replace_line(ms,hm,cm,ip,ns,ms-1,bn)

      elseif(basis.eq.'db') then

       ip(ns)=0
       bn = n * bspd(nv+1,1,ks-1) 
       Call Replace_line(ms,hm,cm,ip,ns-1,ms,bn)

      elseif(basis.eq.'aa') then

       dn=bspd(nv+1,1,ks  ) *n
       bn=bspd(nv+1,1,ks-1) *n

       Call Replace_line(ms,hm,cm,ip,ns  ,ms,dn)
       Call Replace_line(ms,hm,cm,ip,ns-1,ms,bn)

      else

       Stop 'Bloch: unknown basis'

      end if

      End Subroutine Bloch


!=======================================================================
      SUBROUTINE Replace_line(ms,hm,cm,ip,i,j,C)
!=======================================================================
!     Bloch-operator modification of Dirac Hamiltonian
!     in B-spline basis
!-----------------------------------------------------------------------
      Implicit none
      Integer ::  ms, i,j, ip(ms)
      Real(8) ::  hm(ms,ms),cm(ms,ms),C

      ip(j)=0

      hm(i,:) = hm(i,:) + hm(j,:) * C
      cm(i,:) = cm(i,:) + cm(j,:) * C

      hm(:,i) = hm(:,i) + hm(:,j) * C
      cm(:,i) = cm(:,i) + cm(:,j) * C

      END SUBROUTINE Replace_line

