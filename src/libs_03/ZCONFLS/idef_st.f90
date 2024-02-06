!=======================================================================
      Integer Function Idef_st(in)
!=======================================================================
!     define number of states in l- or j-files (unit 'in')
!-----------------------------------------------------------------------
      Character AS*80

      Idef_st = 0
      rewind(in)

    1 read(in,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(15:20).ne.'NUMBER') go to 1
      read(AS(23:),*) n
      Idef_st=Idef_st+n
      go to 1
    2 Rewind(in)

      End Function Idef_st

