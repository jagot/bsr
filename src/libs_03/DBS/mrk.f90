subroutine mrk_common_gen_array(ksi, ksj, rkd1, rkd2, rkd3, rkd4)
  Use DBS_grid
  Use DBS_moments, only: rkd
  Use DBS_integrals
  Use Timer

  Implicit none
  Integer, intent(in) :: ksi, ksj
  Real(8), intent(in) :: rkd1(1:ks*ks,1:nv), rkd2(1:ks*ks,1:nv), rkd3(1:ks*ks,1:nv), rkd4(1:ks*ks,1:nv)

  Integer :: i,j, ii,jj, iv,jv, ih,jh, ihp,jhp, ip,jp
  Real(8) :: c

  call TimerStart('mrk: generate array')
  rkb=0.d0

  do jv = 1,nv
     jj = 0
     do jh = 1,ksj
        j  = jv  + jh - 1
        do jhp=jh,ksj
           jp = jhp - jh + 1
           jj = jj  + 1

           do iv = 1,nv
              ii = 0
              do ih = 1,ksi
                 i  = iv  + ih - 1
                 do ihp=ih,ksi
                    ip = ihp - ih + 1
                    ii = ii  + 1

                    if     ( iv < jv ) then
                       c = rkd1(ii,iv)*rkd2(jj,jv)
                    else if( iv > jv ) then
                       c = rkd3(jj,jv)*rkd4(ii,iv)
                    else
                       c = rkd(ii,jj,iv)
                    end if

                    rkb(i,j,ip,jp) = rkb(i,j,ip,jp) +  c
                 end do
              end do
           end do
        end do
     end do
  end do
  call TimerStop('mrk: generate array')
end subroutine mrk_common_gen_array
