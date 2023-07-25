!======================================================================
      Subroutine tma_om (ilsp,nopen,nom,tmar,tmai,om)
!======================================================================
      Use target
      Use channels

      Implicit none
      Integer, Intent(in)  :: ilsp,nopen
      Real(8), intent(in)  :: tmar(*),tmai(*)
      Integer, Intent(out) :: nom
      Real(8), intent(out) :: om(*)

      Integer :: i,j, ij, itr,itr1,itr2
      Real(8) :: g,s

      g = (2*lpar(ilsp)+1) * iabs(ispar(ilsp)) / 2.d0
      if(ispar(ilsp).eq.0) g = (lpar(ilsp)+1)/2.d0
      
      i=iptar(ilsp,nopen);  nom=(i+1)*i/2; om(1:nom)=0.d0

      Do i=1,nopen; itr1=iptar(ilsp,i)
       Do j=i,nopen; itr2=iptar(ilsp,j)
        ij=(j-1)*j/2+i; itr=(itr2-1)*itr2/2+itr1
        s = (tmar(ij)*tmar(ij)+tmai(ij)*tmai(ij))*g
        if(itr1.eq.itr2.and.i.ne.j) s = s + s
        om(itr) = om(itr) + s
       End do
      End do

      End Subroutine tma_om
