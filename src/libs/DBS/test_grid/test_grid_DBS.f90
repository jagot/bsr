      USE DBS_grid

      Implicit real(8) (A-H,O-Z)

      Call Def_grid

      nu=1; open(nu,file='knot.tab')

      Do i=1,ns+ks
       write(nu,'(i8,D24.16)')  i,t(i)
      End do


      End 

