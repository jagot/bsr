!=======================================================================
       Integer Function Index_TR(ion,i1,i2,np,ni)
!=======================================================================
!      index of matrix element OM(i1,i2) in one-dimention OM-matrix array
!      used UPPER part of matrix !!!
!      ion - neutral atom or ion (0 or > 0)
!      np - number of physical states
!      ni - number of 'ionization' states
!----------------------------------------------------------------------

       Index_TR  = 0
       if(ion.ne.0.and.i1.eq.i2) Return
       if(i2.lt.i1) Return
       if(i1.gt.ni.and.i2.gt.np) Return

       if(i2.le.np) then
        itr = (i2-1)*i2/2+i1
        if(ion.ne.0) itr = itr - i2 + 1
       else
        itr = (np+1)*np/2 + (i2-np-1)*ni + i1
        if(ion.ne.0)  itr = (np-1)*np/2 + (i2-np-1)*ni + i1
       end if

       Index_TR = itr

       End Function Index_TR

