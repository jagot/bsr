!====================================================================
      SUBROUTINE Mod_bsp(kappa,m)
!====================================================================
!  ... some new basis ???
!--------------------------------------------------------------------

      Use DBS_grid
      Use DBS_gauss
      Use zconst, only: c_au
      Use DBS_nuclear, Z => atomic_number

      Implicit none
      Integer, Intent(in) :: kappa,m
      Integer :: i,j,nw
      Real(8) :: gamma

      if(m.le.0) Return      

      gamma = sqrt(kappa**2 - (z/c_au)**2)

      Do i=1,nv            ! intervals
       Do j=1,ks           ! B-splines
        if(i+j-1.gt.m) Cycle
        gx(:) = gr(i,:)**gamma
        dpbs(i,:,j) = (dpbs(i,:,j)+gamma*bsp(i,:,j)*grm(i,:)) * gx(:)
        dmbs(i,:,j) = (dmbs(i,:,j)-gamma*bsp(i,:,j)*grm(i,:)) * gx(:)
        bspd(i,:,j) = (bspd(i,:,j)+gamma*bsp(i,:,j)*grm(i,:)) * gx(:)
        bsp (i,:,j) = bsp(i,:,j)* gx(:)
       End do
      End do

       i = nv+1
       gx(1) = t(ks)**gamma
       dpbs(i,2,:) = (dpbs(i,1,:)+gamma*bsp(i,1,:)/t(ns+1)) * gx(1)
       dmbs(i,2,:) = (dmbs(i,1,:)-gamma*bsp(i,1,:)/t(ns+1)) * gx(1)
       bspd(i,2,:) = (bspd(i,1,:)+gamma*bsp(i,1,:)/t(ns+1)) * gx(1)
       bsp (i,2,:) = bsp(i,1,:)* gx(1)

      if(m.ge.ns) then
       i = nv+1
       gx(1) = t(ns+1)**gamma
       dpbs(i,1,:) = (dpbs(i,1,:)+gamma*bsp(i,1,:)/t(ns+1)) * gx(1)
       dmbs(i,1,:) = (dmbs(i,1,:)-gamma*bsp(i,1,:)/t(ns+1)) * gx(1)
       bspd(i,1,:) = (bspd(i,1,:)+gamma*bsp(i,1,:)/t(ns+1)) * gx(1)
       bsp (i,1,:) = bsp(i,1,:)* gx(1)
      end if

      End Subroutine Mod_bsp



