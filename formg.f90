
module hamilsccmod 
use datamod
use utilmod
!----------------------------------------------------------------------------------------------
!  shell resolved versions of the formG
interface formG_shrs
   module procedure formG_shrs_SP           !  <-- input: S,P 
   module procedure formG_shrs_Q            !  <-- input: dQ
   module procedure formG_shrs_SPGamm       !  <-- input: S,P, Gamma-matrix
   module procedure formG_shrs_QGamm        !  <---input:  dQ, Gamma-matrix
end interface formG_shrs                     

!---------------------
! atom resolved versions of the FormG
interface formG_atrs
   module procedure formG_atrs_SP           !  <-- input: S,P   
   module procedure formG_atrs_Q            !  <-- input: dQ                  (G,nbasis,C1,iat1,Nat1,C2,iat2,Nat2,dQ)
   module procedure formG_atrs_SPGamm       !  <-- input: S,P, Gamma-matrix  
   module procedure formG_atrs_QGamm        !  <-- input:  dQ, Gamma-matrix   (G,nbasis,C1,iat1,Nat1,C2,iat2,Nat2,dQ,gammat2)
end interface formG_atrs                              


!----
!interface  formgammat 
!    module  procedure  formgammat_shrs      ! (gammat,C1,iat1,Nat1,nbasis1,C2,iat2,Nat2,nbasis2) 
!    module  procedure  formgammat_atrs      ! (gammat,C1,iat1,Nat1,Nbasis1,C2,iat2,Nat2,Nbasis2)
!end interface  formgammat 


!----------------------------------------------------------------------------------------------
contains 
  subroutine  formg_shrs_SP(G,nbasis,C1,iat1,Nat1,C2,iat2,Nat2,S,P) 
    implicit none
    real (kind=8)       :: C1(3,Nat1), C2(3,Nat2), S(Nbasis,Nbasis), G(nbasis,nbasis) , P(nbasis,nbasis) 
    integer             ::  i,j,Nat1,Nat2, nbasis,iat1(:),iat2(:)
    real (kind=8)       ::  R,x1,y1,z1,x2,y2,z2,dx,dy,dz , ddot, qtot, gtmp
    real (kind=8),allocatable  ::  qatom(:),q0atom(:), dqatom(:),qshell(:),q0shell(:),dqshell(:)       ! Mulliken
    real (kind=8),allocatable  ::  gammat(:,:),  gammxq(:)                 ! JJ-version shell -resolution gamma and G ,(kloppman-ohno)
    integer                    ::  nbasis1,nbasis2,iact ,ii,jj
    
    nbasis1=nbasis  ;   nbasis2=nbasis  
    allocate (gammat(nbasis1,nbasis2), gammxq(nbasis))                 ! full shelli-resolution  versions
    gammat=0.d0   
    call formgammat_shrs(gammat,C1,iat1,Nat1,nbasis1,C2,iat2,Nat2,nbasis2)      ! gamma_JJ

    !
    !  get mulliken charges:
    !  reference charges  for neutral atoms/shells  and calculate  difference charges 
    ! q0shell, q0atom,  dqatom, dqshell
     allocate(qatom(Nat1),    q0atom(Nat1),    dqatom(Nat1))
     allocate(qshell(nbasis), q0shell(nbasis), dqshell(nbasis))

     do i=1,Nat1   
       q0atom(i) = atoms(iat1(i))%fs  + atoms(iat1(i))%fp  +atoms(iat1(i))%fd
       enddo 
    
    !--------------
    !   Mulliken
     qtot=0.0d0
     do i =1,Nbasis
       qshell(i) = ddot(Nbasis,S(1,i),1, P(1,i),1)
       qtot= qtot +qshell(i)
     enddo
     
    !----------
    !    set q0shell and qatom
    qatom=0.0d0    

    iact=0
    do i=1,Nat1
       do j=1,lmax(iat1(i))* lmax(iat1(i))                            !  loop over shells
           if (j==1)                then                              !  s-shell
               q0shell(iact +j)   = atoms(iat1(i))%fs 
           elseif (j>1 .and. j<5)   then                              !  p-shell
               q0shell(iact +j)   = atoms(iat1(i))%fp/3.0d0
           elseif (j>4 .and. j<10)  then                              !  d-shell
               q0shell(iact +j)   = atoms(iat1(i))%fd/5.0d0
           endif 
           qatom(i) = qatom(i) + qshell(iact + j)
       enddo
       iact =iact + lmax(iat1(i))* lmax(iat1(i))
    enddo

    !--------
    !    difference charges
        dqatom  = qatom  - q0atom  
        dqshell = qshell - q0shell 

    !
    ! now  the   G -matrix, shell-resolution  (JJ version),  Kloppman-Ohno
     gammxq=0.0d0 
     call  dgemv('N',nbasis,nbasis, 1.0d0, gammat,nbasis,dqshell,1, 0.0d0, gammxq,1)
     do i=1,nbasis 
       do j=1,nbasis 
         G(i,j) =  0.5d0*S(i,j)*(gammxq(i) +gammxq(j))
       enddo
     enddo
   
  end subroutine    formg_shrs_SP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  formg_shrs_SPGamm(G,nbasis,C1,iat1,Nat1,C2,iat2,Nat2,S,P,gammat) 
    implicit none
    real (kind=8)       :: C1(3,Nat1), C2(3,Nat2), S(Nbasis,Nbasis), G(nbasis,nbasis) , P(nbasis,nbasis) 
    integer             ::  i,j,Nat1,Nat2, nbasis,iat1(:),iat2(:)
    real (kind=8)       ::  R,x1,y1,z1,x2,y2,z2,dx,dy,dz , ddot, qtot, gtmp , gammat(nbasis,nbasis)
    real (kind=8),allocatable  ::  qatom(:),q0atom(:), dqatom(:),qshell(:),q0shell(:),dqshell(:)       ! Mulliken
    real (kind=8),allocatable  ::  gammxq(:)                 ! JJ-version shell -resolution gamma and G ,(kloppman-ohno)
    integer                    ::  nbasis1,nbasis2,iact ,ii,jj
    
    nbasis1=nbasis  ;   nbasis2=nbasis  
    allocate (gammxq(nbasis))                 ! full shelli-resolution  versions

    !
    !  get mulliken charges:
    !  reference charges  for neutral atoms/shells  and calculate  difference charges 
    ! q0shell, q0atom,  dqatom, dqshell
     allocate(qatom(Nat1),    q0atom(Nat1),    dqatom(Nat1))
     allocate(qshell(nbasis), q0shell(nbasis), dqshell(nbasis))

     do i=1,Nat1   
       q0atom(i) = atoms(iat1(i))%fs  + atoms(iat1(i))%fp  +atoms(iat1(i))%fd
       enddo 
    
    !--------------
    !   Mulliken
     qtot=0.0d0
     do i =1,Nbasis
       qshell(i) = ddot(Nbasis,S(1,i),1, P(1,i),1)
       qtot= qtot +qshell(i)
     enddo
     
    !----------
    !    set q0shell and qatom
    qatom=0.0d0    

    iact=0
    do i=1,Nat1
       do j=1,lmax(iat1(i))* lmax(iat1(i))                            !  loop over shells
           if (j==1)                then                              !  s-shell
               q0shell(iact +j)   = atoms(iat1(i))%fs 
           elseif (j>1 .and. j<5)   then                              !  p-shell
               q0shell(iact +j)   = atoms(iat1(i))%fp/3.0d0
           elseif (j>4 .and. j<10)  then                              !  d-shell
               q0shell(iact +j)   = atoms(iat1(i))%fd/5.0d0
           endif 
           qatom(i) = qatom(i) + qshell(iact + j)
       enddo
       iact =iact + lmax(iat1(i))* lmax(iat1(i))
    enddo

    !--------
    !    difference charges
        dqatom  = qatom  - q0atom  
        dqshell = qshell - q0shell 

    !
    ! now  the   G -matrix, shell-resolution  (JJ version),  Kloppman-Ohno
     gammxq=0.0d0 
     call  dgemv('N',nbasis,nbasis, 1.0d0, gammat,nbasis,dqshell,1, 0.0d0, gammxq,1)

     do i=1,nbasis 
       do j=1,nbasis 
         G(i,j) =  0.5d0*S(i,j)*(gammxq(i) +gammxq(j))
       enddo
     enddo
   
   
  end subroutine  formg_shrs_SPGamm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  formg_shrs_Q(G,nbasis,C1,iat1,Nat1,C2,iat2,Nat2,S,dQ) 
    implicit none
    real (kind=8)       :: C1(3,Nat1), C2(3,Nat2), S(Nbasis,Nbasis), G(nbasis,nbasis) , P(nbasis,nbasis) 
    integer             ::  i,j,Nat1,Nat2, nbasis,iat1(:),iat2(:)
    real (kind=8)       ::  R,x1,y1,z1,x2,y2,z2,dx,dy,dz , ddot, qtot, gtmp,dQ(:)
    real (kind=8),allocatable  ::  qatom(:),q0atom(:), dqatom(:),qshell(:),q0shell(:)       ! Mulliken
    real (kind=8),allocatable  ::  gammat(:,:),  gammxq(:)                 ! JJ-version shell -resolution gamma and G ,(kloppman-ohno)
    integer                    ::  nbasis1,nbasis2,iact ,ii,jj
    
    nbasis1=nbasis  ;   nbasis2=nbasis  
    allocate (gammat(nbasis1,nbasis2), gammxq(nbasis))                 ! full shelli-resolution  versions

    gammat=0.d0   
    call formgammat_shrs(gammat,C1,iat1,Nat1,nbasis1,C2,iat2,Nat2,nbasis2)      ! gamma_JJ

    !
    ! now  the   G -matrix, shell-resolution  (JJ version),  Kloppman-Ohno
     gammxq=0.0d0 
     call  dgemv('N',nbasis,nbasis, 1.0d0, gammat,nbasis,dq,1, 0.0d0, gammxq,1)
     do i=1,nbasis 
       do j=1,nbasis 
         G(i,j) =  0.5d0*S(i,j)*(gammxq(i) +gammxq(j))
       enddo
     enddo
   
  end subroutine  formg_shrs_Q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  formg_shrs_QGamm(G,nbasis,C1,iat1,Nat1,C2,iat2,Nat2,S,dQ,gammat) 
    implicit none
    real (kind=8)       :: C1(3,Nat1), C2(3,Nat2), S(Nbasis,Nbasis), G(nbasis,nbasis) , P(nbasis,nbasis) 
    integer             ::  i,j,Nat1,Nat2, nbasis,iat1(:),iat2(:)
    real (kind=8)       ::  R,x1,y1,z1,x2,y2,z2,dx,dy,dz , ddot, qtot, gtmp,dQ(:),gammat(nbasis,nbasis)
    real (kind=8),allocatable  ::  gammxq(:)                 ! JJ-version shell -resolution gamma and G ,(kloppman-ohno)
    integer                    ::  nbasis1,nbasis2,iact ,ii,jj
    
    nbasis1=nbasis  ;   nbasis2=nbasis  
    allocate (gammxq(nbasis))                 ! full shelli-resolution  versions

    !
    ! now  the   G -matrix, shell-resolution  (JJ version),  Kloppman-Ohno
     gammxq=0.0d0 
     call  dgemv('N',nbasis,nbasis, 1.0d0, gammat,nbasis,dq,1, 0.0d0, gammxq,1)
     do i=1,nbasis 
       do j=1,nbasis 
         G(i,j) =  0.5d0*S(i,j)*(gammxq(i) +gammxq(j))
       enddo
     enddo
   
   
  end subroutine  formg_shrs_QGamm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! atomic resolution  charges only  (M.Elst. version)
!
  subroutine  formG_atrs_SP(G,nbasis,C1,iat1,Nat1,C2,iat2,Nat2,S,P) 
    implicit none
    real (kind=8)       :: C1(3,Nat1), C2(3,Nat2), S(Nbasis,Nbasis), G(nbasis,nbasis) , P(nbasis,nbasis) 
    integer             ::  i,j,Nat1,Nat2, nbasis,iat1(:),iat2(:)
    real (kind=8)       ::  R,x1,y1,z1,x2,y2,z2,dx,dy,dz , ddot, qtot, gtmp
    real (kind=8),allocatable  ::  qatom(:),q0atom(:), dqatom(:),qshell(:),q0shell(:),dqshell(:)       ! Mulliken
    real (kind=8),allocatable  ::  gammat2(:,:), gammxq2(:),scr(:)         ! MElst-version   atomic resolution (Kloppman -Ohno)
    integer                    ::  nbasis1,nbasis2,iact ,ii,jj,M
    
    
    nbasis1=nbasis  ; nbasis2=nbasis  
    allocate (gammat2(nat1,Nat2)     , gammxq2(nbasis), scr(nbasis))   ! atomic resolution version  (standards Marcus)
    call formgammat_atrs(gammat2,C1,iat1,Nat1,nbasis1,C2,iat2,Nat2,nbasis2)      ! gamma_ME

    !
    !  get mulliken charges:
    !  reference charges  for neutral atoms/shells  and calculate  difference charges 
     allocate(qatom(Nat1),    q0atom(Nat1),    dqatom(Nat1))
     allocate(qshell(nbasis), q0shell(nbasis), dqshell(nbasis))

     do i=1,Nat1   
       q0atom(i) = atoms(iat1(i))%fs  + atoms(iat1(i))%fp  +atoms(iat1(i))%fd
       enddo 

    ! set qshell
    
     qtot=0.0d0
     do i =1,Nbasis
       qshell(i) = ddot(Nbasis,S(1,i),1, P(1,i),1)
       qtot= qtot +qshell(i)
     enddo
     
    !
    ! set q0shell and qatom
    
    qatom=0.0d0    

    iact=0
    do i=1,Nat1
       do j=1,lmax(iat1(i))* lmax(iat1(i))                            !  loop over shells
           if (j==1)                then                              !  s-shell
               q0shell(iact +j)   = atoms(iat1(i))%fs 
           elseif (j>1 .and. j<5)   then                              !  p-shell
               q0shell(iact +j)   = atoms(iat1(i))%fp/3.0d0
           elseif (j>4 .and. j<10)  then                              !  d-shell
               q0shell(iact +j)   = atoms(iat1(i))%fd/5.0d0
           endif 
           qatom(i) = qatom(i) + qshell(iact + j)
       enddo
       iact =iact + lmax(iat1(i))* lmax(iat1(i))
    enddo

    !
    ! difference charges
        dqatom  = qatom  - q0atom  
        dqshell = qshell - q0shell 

  ! now  the   G -matrix, atomic-resolution  (M.Elst. version),  Kloppman-Ohno
     gammxq2=0.0d0
  !shift of hamiltonian
     call  dgemv('N',nat1,nat1, 1.0d0, gammat2,nat1,dqatom,1, 0.0d0, gammxq2,1)

     iact=1
     do i=1,nat1
        M =lmax(iat1(i))* lmax(iat1(i))
        scr(iact:iact+M-1) =  gammxq2(i)
        iact = iact + M
     enddo

     do i=1,nbasis
       do j=1,nbasis
         G(i,j) =  0.5d0*S(i,j)*(scr(i)+scr(j))
       enddo
     enddo

     end subroutine formG_atrs_SP
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  formG_atrs_SPGamm(G,nbasis,C1,iat1,Nat1,C2,iat2,Nat2,S,P,gammat2) 
    implicit none
    real (kind=8)       :: C1(3,Nat1), C2(3,Nat2), S(Nbasis,Nbasis), G(nbasis,nbasis), P(nbasis,nbasis), gammat2(nat1,nat2)
    integer             ::  i,j,Nat1,Nat2, nbasis,iat1(:),iat2(:)
    real (kind=8)       ::  R,x1,y1,z1,x2,y2,z2,dx,dy,dz , ddot, qtot, gtmp , ecoul
    real (kind=8),allocatable  ::  qatom(:),q0atom(:), dqatom(:),qshell(:),q0shell(:),dqshell(:)       ! Mulliken
    !real (kind=8),allocatable  ::  gammat(:,:),  gammxq(:)                 ! JJ-version shell -resolution gamma and G ,(kloppman-ohno)
    real (kind=8),allocatable  ::  gammxq2(:),scr(:)         ! MElst-version   atomic resolution (Kloppman -Ohno)
    integer                    ::  nbasis1,nbasis2,iact ,ii,jj,M
    
    nbasis1=nbasis   ; nbasis2=nbasis  
    allocate (gammxq2(nbasis), scr(nbasis))   ! atomic resolution version  (standards Marcus)

    !
    !  get mulliken charges:
    !  reference charges  for neutral atoms/shells  and calculate  difference charges 
     allocate(qatom(Nat1),    q0atom(Nat1),    dqatom(Nat1))
     allocate(qshell(nbasis), q0shell(nbasis), dqshell(nbasis))

     do i=1,Nat1   
       q0atom(i) = atoms(iat1(i))%fs  + atoms(iat1(i))%fp  +atoms(iat1(i))%fd
       enddo 

    ! set qshell
    
     qtot=0.0d0
     do i =1,Nbasis
       qshell(i) = ddot(Nbasis,S(1,i),1, P(1,i),1)
       qtot= qtot +qshell(i)
     enddo
     
    !
    ! set q0shell and qatom
    
    qatom=0.0d0    

    iact=0
    do i=1,Nat1
       do j=1,lmax(iat1(i))* lmax(iat1(i))                            !  loop over shells
           if (j==1)                then                              !  s-shell
               q0shell(iact +j)   = atoms(iat1(i))%fs 
           elseif (j>1 .and. j<5)   then                              !  p-shell
               q0shell(iact +j)   = atoms(iat1(i))%fp/3.0d0
           elseif (j>4 .and. j<10)  then                              !  d-shell
               q0shell(iact +j)   = atoms(iat1(i))%fd/5.0d0
           endif 
           qatom(i) = qatom(i) + qshell(iact + j)
       enddo
       iact =iact + lmax(iat1(i))* lmax(iat1(i))
    enddo

    !
    ! difference charges
        dqatom  = qatom  - q0atom  
        dqshell = qshell - q0shell 

  ! now  the   G -matrix, atomic-resolution  (M.Elst. version),  Kloppman-Ohno
     gammxq2=0.0d0
  !shift of hamiltonian
     call  dgemv('N',nat1,nat1, 1.0d0, gammat2,nat1,dqatom,1, 0.0d0, gammxq2,1)
     ecoul = ddot(nat1*nat1,gammxq2,1,dqatom,1) 
     write(*,*)'current Ecoul =',ecoul
     write(*,*)'dqatom=',dqatom

     iact=1
     do i=1,nat1
        M =lmax(iat1(i))* lmax(iat1(i))
        scr(iact:iact+M-1) =  gammxq2(i)
        iact = iact + M
     enddo

     do i=1,nbasis
       do j=1,nbasis
         G(i,j) =  0.5d0*S(i,j)*(scr(i)+scr(j))
       enddo
     enddo

     end subroutine formG_atrs_SPGamm
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  formG_atrs_Q(G,nbasis,C1,iat1,Nat1,C2,iat2,Nat2,S,dQ)
    implicit none
    real (kind=8)       :: C1(3,Nat1), C2(3,Nat2), S(Nbasis,Nbasis), G(nbasis,nbasis) , P(nbasis,nbasis) ,dQ(nat1)
    integer             ::  i,j,Nat1,Nat2, nbasis,iat1(:),iat2(:)
    real (kind=8)       ::  R,x1,y1,z1,x2,y2,z2,dx,dy,dz , ddot, qtot, gtmp
    real (kind=8),allocatable  ::  qatom(:),q0atom(:), dqatom(:),qshell(:),q0shell(:),dqshell(:)       ! Mulliken
    !real (kind=8),allocatable  ::  gammat(:,:),  gammxq(:)                 ! JJ-version shell -resolution gamma and G ,(kloppman-ohno)
    real (kind=8),allocatable  ::  gammat2(:,:), gammxq2(:),scr(:)         ! MElst-version   atomic resolution (Kloppman -Ohno)
    integer                    ::  nbasis1,nbasis2,iact ,ii,jj,M
    
    !------
    !   calculate gamma matrix, Kloppman -Ohno
    nbasis1=nbasis   ; nbasis2=nbasis  
    allocate (gammat2(nat1,Nat2)     , gammxq2(nbasis), scr(nbasis))   ! atomic resolution version  (standards Marcus)
    call formgammat_atrs(gammat2,C1,iat1,Nat1,nbasis1,C2,iat2,Nat2,nbasis2)      ! gamma_ME

  ! now  the   G -matrix, atomic-resolution  (M.Elst. version),  Kloppman-Ohno
     gammxq2=0.0d0
     call  dgemv('N',nat1,nat1, 1.0d0, gammat2,nat1,dQ,1, 0.0d0, gammxq2,1)

     iact=1
     do i=1,nat1
        M =lmax(iat1(i))* lmax(iat1(i))
        scr(iact:iact+M-1) =  gammxq2(i)
        iact = iact + M
     enddo

     do i=1,nbasis
       do j=1,nbasis
         G(i,j) =  0.5d0*S(i,j)*(scr(i)+scr(j))
       enddo
     enddo

     end subroutine formG_atrs_Q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  formG_atrs_QGamm(G,nbasis,C1,iat1,Nat1,C2,iat2,Nat2,S,dQ,gammat2)
    implicit none
    real (kind=8)       :: C1(3,Nat1), C2(3,Nat2), S(Nbasis,Nbasis),G(nbasis,nbasis),P(nbasis,nbasis),dQ(nat1),gammat2(nat1,nat2)
    integer             ::  i,j,Nat1,Nat2, nbasis,iat1(:),iat2(:)
    real (kind=8)       ::  R,x1,y1,z1,x2,y2,z2,dx,dy,dz , ddot, qtot, gtmp
    real (kind=8),allocatable  ::  qatom(:),q0atom(:), dqatom(:),qshell(:),q0shell(:),dqshell(:)       ! Mulliken
    !real (kind=8),allocatable  ::  gammat(:,:),  gammxq(:)                 ! JJ-version shell -resolution gamma and G ,(kloppman-ohno)
    real (kind=8),allocatable  ::  gammxq2(:),scr(:)         ! MElst-version   atomic resolution (Kloppman -Ohno)
    integer                    ::  nbasis1,nbasis2,iact ,ii,jj,M
    
    !------
    !   calculate  hamiltonian shift from gamma matrix, Kloppman -Ohno
     !write(*,*)'entering formG_atrs_QGamm:'
     !write(*,*)'*** dQ='; call  dumpm(dQ)
     nbasis1=nbasis   ; nbasis2=nbasis  
     allocate (gammxq2(nbasis), scr(nbasis))   ! atomic resolution version  (standards Marcus)
     gammxq2=0.0d0
     call  dgemv('N',nat1,nat1, 1.0d0, gammat2,nat1,dQ,1, 0.0d0, gammxq2,1)
     iact=1
     do i=1,nat1
        M =lmax(iat1(i))* lmax(iat1(i))
        scr(iact:iact+M-1) =  gammxq2(i)
        iact = iact + M
     enddo

     do i=1,nbasis
       do j=1,nbasis
         G(i,j) =  0.5d0*S(i,j)*(scr(i)+scr(j))
       enddo
     enddo
     !stop

     end subroutine formG_atrs_QGamm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  subroutine formgammat_shrs(gammat,C1,iat1,Nat1,Nbasis1,C2,iat2,Nat2,Nbasis2)
    implicit none 
    real (kind=8)        ::  C1(3,Nat1), C2(3,Nat2)
    real (kind=8),intent(out)         ::  gammat(nbasis1,nbasis2)
    integer              ::  i,j,Nat1,Nat2, nbasis1,nbasis2,iat1(:),iat2(:)
    real (kind=8)        ::  R,x1,y1,z1,x2,y2,z2,dx,dy,dz 
    real (kind=8)        ::  Uas,Uap,Uad,  Ubs,Ubp,Ubd , one

    integer              :: mapat1,mapat2,  ntmp1, ntmp2,iact ,iact1,iact2 ,ii, jj
    integer, allocatable      :: lmvec1(:), lmvec2(:) 
    real (kind=8),allocatable :: Uvec1i(:),  Uvec2i(:)

    
    one =1.0d0
!
!  check the input size
!   
   ntmp1 = 0   ; ntmp2=0
   do i=1,Nat1
     ntmp1 = ntmp1 + lmax(iat1(i))* lmax(iat1(i))
     enddo
   do i=1,Nat2
     ntmp2 = ntmp2 + lmax(iat1(i))* lmax(iat1(i))
     enddo
   if (ntmp1.ne.nbasis1 .or.  ntmp2.ne.nbasis2) then 
     write(*,*)'Problem with gammat  size. Stop'
     write(*,*)'nbas1, nsize1 =', nbasis1, ntmp1
     write(*,*)'nbas2, nsize2 =', nbasis2, ntmp2
     stop
   endif 


!@
!@ now set Uvec1,lmvec1
!@  Uvec  <-- vector of   U chemical hardness for C
!@  lmvec <-- vector that holds  basis set info (where is  the last  basis for the C1) 
!@
    allocate(lmvec1(nat1),lmvec2(nat2), Uvec1i(nbasis1), Uvec2i(nbasis2))
    Uvec1i =0.0d0   ; Uvec2i =0.d0  ; lmvec1= 0 ; lmvec2 = 0
    iact = 0
    do i=1,Nat1
      ntmp1= lmax(iat1(i))*lmax(iat1(i))
      lmvec1(i) = iact+ ntmp1 
        Uvec1i(iact+1) = one/atoms(iat1(i))%Us
        if  (ntmp1>1)  then 
          Uvec1i(iact+2) = one/atoms(iat1(i))%Up
          Uvec1i(iact+3) = one/atoms(iat1(i))%Up
          Uvec1i(iact+4) = one/atoms(iat1(i))%Up
        endif 
        if  (ntmp1>4)  then 
          Uvec1i(iact+5) = one/atoms(iat1(i))%Ud
          Uvec1i(iact+6) = one/atoms(iat1(i))%Ud
          Uvec1i(iact+7) = one/atoms(iat1(i))%Ud
          Uvec1i(iact+8) = one/atoms(iat1(i))%Ud
          Uvec1i(iact+9) = one/atoms(iat1(i))%Ud
        endif 
      iact =iact +ntmp1
    enddo

!@
!@ same as above but now set Uvec2,lmvec2
!@  Uvec  <-- vector of   U chemical hardness for C
!@  lmvec <-- vector that holds  basis set info (where is  the last  basis for the C1) 
!@
    iact = 0
    do i=1,Nat2
      ntmp2= lmax(iat2(i))*lmax(iat2(i))
      lmvec2(i) = iact+ ntmp2 
        Uvec2i(iact+1) = one/atoms(iat2(i))%Us
        if  (ntmp2>1)  then 
          Uvec2i(iact+2) = one/atoms(iat2(i))%Up
          Uvec2i(iact+3) = one/atoms(iat2(i))%Up
          Uvec2i(iact+4) = one/atoms(iat2(i))%Up
        endif 
        if  (ntmp2>4)  then 
          Uvec2i(iact+5) = one/atoms(iat2(i))%Ud
          Uvec2i(iact+6) = one/atoms(iat2(i))%Ud
          Uvec2i(iact+7) = one/atoms(iat2(i))%Ud
          Uvec2i(iact+8) = one/atoms(iat2(i))%Ud
          Uvec2i(iact+9) = one/atoms(iat2(i))%Ud
        endif 
      iact =iact +ntmp2
    enddo

!@
!@  now do the gamma matrix
!@
    iact1 =1
    iact2 =1
    do i=1,Nat1
      x1 = C1(1,i)  
      y1 = C1(2,i)  
      z1 = C1(3,i)  
      do j=1,Nat2
         dx=x1-C2(1,j)
         dy=y1-C2(2,j)
         dz=z1-C2(3,j)
         R=sqrt(dx*dx +dy*dy+dz*dz)      !  R is expected to be already in  bohrs
         
         do ii=iact1,   lmvec1(i)  
           do jj=iact2, lmvec2(j)  
             gammat(ii,jj)   =1.0d0 /sqrt( R*R +0.25d0*(Uvec1i(ii)  + Uvec2i(jj))**2)
           enddo
         enddo
         iact2 = lmvec2(j) +1
      enddo  
      iact2= 1
      iact1 = lmvec1(i) +1
        
    enddo  
    
   100  format(2(' ',I6),6(' ',E12.6))
  end subroutine formgammat_shrs



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this  version uses total atomic charges instead  of  shells
!  
  subroutine formgammat_atrs(gammat,C1,iat1,Nat1,Nbasis1,C2,iat2,Nat2,Nbasis2)
    implicit none 
    real (kind=8)        ::  C1(3,Nat1), C2(3,Nat2)
    real (kind=8), intent (out)        ::  gammat(nat1,Nat2)
    integer              ::  i,j,Nat1,Nat2, nbasis1,nbasis2,iat1(:),iat2(:)
    real (kind=8)        ::  R,x1,y1,z1,x2,y2,z2,dx,dy,dz 
    real (kind=8)        ::  Uas,Uap,Uad,  Ubs,Ubp,Ubd ,one

    integer              :: mapat1,mapat2,  ntmp1, ntmp2,iact ,iact1,iact2 
    real (kind=8),allocatable :: Uvec1i(:),  Uvec2i(:)

   one =1.0d0
!@
!@ now set Uvec1,lmvec1
!@  Uvec  <-- vector of   U chemical hardness for C
!@
    allocate(Uvec1i(nbasis1), Uvec2i(nbasis2))
   
    Uvec1i =0.0d0   ; Uvec2i =0.d0  
    iact = 0
    do i=1,Nat1
      Uvec1i(i) = one/atoms(iat1(i))%Us
    enddo
    do i=1,Nat2
      Uvec2i(i) = one/atoms(iat2(i))%Us
    enddo

!@
!@  now do the gamma matrix
!@
    do i=1,Nat1
      x1 = C1(1,i)  
      y1 = C1(2,i)  
      z1 = C1(3,i)  
      do j=1,Nat2
         dx=x1-C2(1,j)
         dy=y1-C2(2,j)
         dz=z1-C2(3,j)
         R=sqrt(dx*dx +dy*dy+dz*dz)      !  R is expected to be already in  bohrs
         gammat(i,j)   =1.0d0 /sqrt( R*R +0.25d0*(Uvec1i(i)  + Uvec2i(j))**2)
      enddo 
    enddo
         
   100  format(2(' ',I6),6(' ',E12.6))
  end subroutine  formgammat_atrs


end module   hamilsccmod




