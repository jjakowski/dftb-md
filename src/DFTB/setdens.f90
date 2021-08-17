module   density 
 use slkomod
 use datamod


interface  setdens 
   module procedure setdens_mo       !   nbasis,    Co,focc,U,Ui,    Pn
   module procedure setdens_mo_v2    !   nbasis,    Co,focc,U,Ui,    Pn, Po
   module procedure setdens_focc     !   nbasis,       focc,         Pn 
   module procedure setdens_qatom    !   nbasis,nat,qatom,iat,       Pn
end interface setdens 


 
contains 
!----------
   subroutine setdens_mo(nbasis,Co,focc,U,Ui,Pn)
     implicit none
     integer      ,intent (in)  :: nbasis
     real (kind=8),intent (in)  :: focc(nbasis),U(nbasis,nbasis),Ui(nbasis,Nbasis),Co(nbasis,nbasis)
     real (kind=8),intent (out) :: Pn(nbasis,nbasis)
     real (kind=8),allocatable  :: Po(:,:)
     integer                    ::   i,j,k ,Nsq
  ! set o/n transformation flags
    integer, parameter          :: iPo2Pn=1, iPn2Po=2, iFo2Fn=3, iFn2Fo=4
   
     allocate (Po(nbasis,nbasis))
 
     call  dcopy(Nbasis*Nbasis, Co,1, Pn,1) 
     do i=1,Nbasis 
        call  dscal(Nbasis,focc(i),Pn(1,i),1)
        enddo
     call  dgemm('N','T',Nbasis,Nbasis,Nbasis,1.0d0,Co,Nbasis,Pn, Nbasis,0.0d0,Po,Nbasis)
     call  ontrans(U,Ui,Po(1,1), Pn(1,1),nbasis,iPo2Pn) 

   end subroutine setdens_mo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine setdens_mo_v2(nbasis,Co,focc,U,Ui,Pn,Po)
     implicit none
     integer      ,intent (in)  :: nbasis
     real (kind=8),intent (in)  :: focc(nbasis),U(nbasis,nbasis),Ui(nbasis,Nbasis),Co(nbasis,nbasis)
     real (kind=8),intent (out) :: Pn(nbasis,nbasis), Po(nbasis, nbasis)
     integer                    ::   i,j,k ,Nsq
  ! set o/n transformation flags
    integer, parameter          :: iPo2Pn=1, iPn2Po=2, iFo2Fn=3, iFn2Fo=4
   
 
     call  dcopy(Nbasis*Nbasis, Co,1, Pn,1) 
     do i=1,Nbasis 
        call  dscal(Nbasis,focc(i),Pn(1,i),1)
        enddo
     call  dgemm('N','T',Nbasis,Nbasis,Nbasis,1.0d0,Co,Nbasis,Pn, Nbasis,0.0d0,Po,Nbasis)
     call  ontrans(U,Ui,Po(1,1), Pn(1,1),nbasis,iPo2Pn) 

   end subroutine setdens_mo_v2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine setdens_focc(nbasis,focc,Pn)
   !here  focc  has all Nel  electrons  
     implicit none
     integer,       intent (in)  ::  nbasis
     real (kind=8), intent (in)  ::  focc(nbasis)
     real (kind=8), intent (out) ::  Pn(nbasis,nbasis)
      Pn=0.0d0 ;
      call  dcopy(nbasis,focc,1,Pn,nbasis+1)
   end subroutine setdens_focc  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine setdens_qatom(nbasis,nat,qatom,iatt,Pn)
    !qatom only has  Delta-Q charges on each atom (no shell resolution)
     implicit none
     integer,       intent (in)  ::  nbasis, Nat,iatt(nat)
     real (kind=8), intent (in)  ::  qatom(Nat)
     real (kind=8), intent (out) ::  Pn(nbasis,nbasis)
     !real (kind=8)   ::  Pn(nbasis,nbasis)
     integer                     ::  i,j,k,jtype
     real (kind=8)               ::  charge
     real (kind=8), parameter    :: one=1.0d0, three=3.0d0, five=5.0d0

      j=1
      Pn=0.0d0 ;
      do i=1,nat
         jtype =iatt(i)
         if (lmax(jtype)==1)    charge =  qatom(i)/1.0d0             
         if (lmax(jtype)==2)    charge =  qatom(i)/4.0d0             
         if (lmax(jtype)==3)    charge =  qatom(i)/9.0d0             
          
         Pn(j,j)     = atoms(jtype)%fs / one   - charge          
         j=j+1
      if (lmax(jtype)>1)   then 
         Pn(j  ,j  ) = atoms(jtype)%fp / three - charge
         Pn(j+1,j+1) = atoms(jtype)%fp / three - charge
         Pn(j+2,j+2) = atoms(jtype)%fp / three - charge 
         j= j+3
         endif
      if (lmax(jtype)>2)   then
         Pn(j  ,j  ) = atoms(jtype)%fd / five  - charge 
         Pn(j+1,j+1) = atoms(jtype)%fd / five  - charge 
         Pn(j+2,j+2) = atoms(jtype)%fd / five  - charge 
         Pn(j+3,j+3) = atoms(jtype)%fd / five  - charge 
         Pn(j+4,j+4) = atoms(jtype)%fd / five  - charge 
         j= j+5
         endif
         

      enddo 
   end subroutine setdens_qatom 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  fermi(nbasis,Tel, qel, eigs, efermi,focc)
  !-----
  ! input:
  !    nbasis
  !    Tel <-- temperature
  !    qel <--  No. of electrons (real, can be fractional)
  !    eig <--  sorted eigenvalues
  !    
  ! output: 
  !   efermi   <--  fermi level at Tel
  !   focc     <--  Fermi occupations a t Tel
  !----------
    implicit none
    integer,        intent (in)  ::  nbasis
    real (kind=8),  intent (in)  ::  Tel, qel,eigs(nbasis) 
    real (kind=8),  intent (out) ::  efermi,  focc(nbasis)

    integer                  :: i,l
    real (kind=8)            :: efu, efl,  charge, beta  
    real (kind=8), parameter :: kT2Eh=3.16683d-6, qtol =1.0d-6, zero =0.0d0

    beta = 1.0d0/(Tel*kT2Eh) 
    charge= zero  
 
    if (Qel< zero .or. Qel> real(nbasis))  then 
        write(*,*)'Fermi: wrong Qel =',  Qel
        stop
    endif

    efl =  eigs(1)       
    efu =  eigs(nbasis)  

    if (efermi <  efl)  then 
       efermi = efl
    else if (efermi >  efu ) then 
       efermi = efu
    endif 
    do i = 1, nbasis 
       focc(i) = 1.0d0/(1.0d0 +exp(beta*(eigs(i)-efermi)))
       charge = charge+focc(i)
    enddo
    if (Abs(charge-Qel) < qtol)  return    ! 


    !    -- update fermi energy  and focc
     do  while (Abs(charge-Qel) > qtol)       
       if (charge  <  Qel )   then  
          efl   =  efermi  
       else if (charge  >  Qel) then 
          efu   =  efermi
       endif   
       efermi =  0.5d0*(efl + efu)
       charge= zero  
       do i = 1, nbasis 
          focc(i) = 1.0d0/(1.0d0 +exp(beta*(eigs(i)-efermi)))
          charge = charge+focc(i)
       enddo
     enddo 
    
   ! normalize 
   !  focc = focc* (Qel/charge) 

  end subroutine fermi
end module density
