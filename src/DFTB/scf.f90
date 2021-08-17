

subroutine scf_diis_shrs(Tel,qeltot,MdMax,Nat,C,iiat,nnbasis,Pn,U,Ui,S,Hcore,Gn,Fn,iprint)
  use  utilmod
  use hamilsccmod
  use density
  use datamod

  !  MdMax <---  max size of DIIS space
  implicit none
!  integer    ::  nnbasis
  integer,       intent (in) :: nnbasis,  Nat,  iiat(Nat) , iprint
  real (kind=8), intent (in) :: Tel, C(3,Nat), qeltot
  real (kind=8)              :: Pn(nnbasis,nnbasis),Gn(nnbasis,nnbasis), Fn(nnbasis,nnbasis)
  real (kind=8)              :: S(nnbasis,nnbasis), U(nnbasis,nnbasis), Ui(nnbasis,nnbasis) 
  real (kind=8), intent (in) :: Hcore(nnbasis,nnbasis)
  

  integer                    ::  i,j,k , iact, Mact,  iter , MaxScf, MdMax
  real (kind=8)              ::  qtot , qel, efermi, En, En0
  real (kind=8), allocatable ::  Ptab(:,:,:),  dPtab(:,:,:), B(:,:)     ! save  density form diagon
  real (kind=8), allocatable ::  PdTab(:,:,:), dPdTab(:,:,:)            !  save DIIS denisties   
  real (kind=8), allocatable ::  x(:), scr(:,:),Pdiis(:,:)    ! diis vector
  real (kind=8), allocatable ::  gammat(:,:),  Fo(:,:),Po(:,:), Co(:,:), eigs(:), focc(:) , work(:) 

! set o/n transformation flags
  integer, parameter  :: iPo2Pn=1, iPn2Po=2, iFo2Fn=3, iFn2Fo=4
  integer             :: lwork, info       
  integer             ::  Nel

  real (kind=8), external    :: ddot
  real (kind=8)              ::  scferr, scftol ,dEtol, dE
  logical                    ::  converged

  
  MdMax=1

  write(*,*)'------- scf_diis_shrs  module -------'

    lwork = 5*nnbasis*nnbasis
    allocate( gammat(nnbasis,nnbasis), Fo(nnbasis,nnbasis), Po(nnbasis,nnbasis),Co(nnbasis,nnbasis))
    allocate(scr(nnbasis,nnbasis),B(MdMax,MdMax),x(MdMax))
    allocate(dPtab(nnbasis,nnbasis,MdMax),Ptab(nnbasis,nnbasis,MdMax))
    allocate(eigs(nnbasis), Work(lwork), focc(nnbasis), Pdiis(nnbasis,nnbasis))
    call  formgammat_shrs(gammat, C,iiat,Nat,nnbasis, C,iiat,Nat,nnbasis)
    call  formg_shrs (Gn,nnbasis,C,iiat,Nat,C,iiat,Nat ,S,Pn,gammat)

    dPtab =0.d00 
    Ptab  =0.0d0

    Fn =  Hcore +  Gn 
    iact=1 ; iter =0 
   ! current energy 
    scr = Hcore +0.5d0*Gn
    En =  ddot(nnbasis*nnbasis, scr,1, Pn,1)    
    En0= 0.0d0

1001 format('Iter,  En, dEn,  err: ',I6, ' ',F16.8,2('  ',E12.6))
   
!---  now  diag ------ check here!!!!!!!!!!
   !Hn =  Hcore + Gn
    call  ontrans(U,Ui,Pn     , Po     ,nnbasis,iPn2Po)
    call  ontrans(U,Ui,Fn     , Fo     ,nnbasis,iFn2Fo)

!   error for !  u  pn  fn ui  -   li fn  pn l !-commutator
!DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
    !scr= matmul(Po,Fo) 
    call dgemm('N','N',nnbasis,nnbasis,nnbasis,1.0d0,Po,nnbasis,Fo,nnbasis,0.0d0,scr,nnbasis)
    scr = scr -transpose(scr)
    !Ptab(:,:,iact) = Po 
    !dPtab(:,:,iact) = scr
    ScfErr =  ddot(nnbasis*nnbasis, scr,1, scr,1)
    B=0.0d0
    !B(iact,iact)  = ScfErr
    scftol=1.0d-8 
    MaxSCF =200  
    !MaxSCF =20  
    iter =0 
    write(*,*)'-------------------------------------------'
    !write(*,1001)  iter, En, En-En0, scferr
    
    iter =1 
    SCFtol=1.0d-8
    dEtol =1.0d-8
    de =En-En0
    Pdiis = Po 
    write(*,*) 'scferr, scfto =', scferr,scftol
    !write(*,*) 'iter, MaxSCF =', iter, Maxscf 
     !Result1 =  .NOT.  Something .AND. Another
    converged =   (SCFerr < SCFtol) .and. (dE <  dEtol)
    !do while  (iter < MaxSCF  .and.  SCFerr > SCFtol) 
    do while  (iter < MaxSCF  .and. .not. converged )
         
        iact =  mod(iter,MdMax) ; if (iact.eq.0) iact=MdMax
        Mact =  min(MdMax,iter)  
        Ptab(:,:,iact) = Po 
        !Ptab(:,:,iact) = Pdiis 
        dPtab(:,:,iact) = scr
        B(iact,iact)  = ScfErr
        do i=1, Mact
          scferr =  ddot(nnbasis*nnbasis, dPtab(1,1,i), 1, dPtab(1,1,iact),1)   
          B(i,iact) = scferr 
          B(iact,i) = scferr 
          enddo
        write(*,*)'-----------------------------------------------------------'
        write(*,*) converged
        !write(*,*)'DIIS: Iact, Mact, MdMax ',  iact, Mact,MdMax
        !write(*,*) '**** B:' ;  call dumpm(B)
        !write(*,*) '**** dPtab :' ;  call dumpm(dPtab )
        !write(*,*) '****  Ptab :' ;  call dumpm(Ptab )
       
       
        !call solve_diis(Mact,B,x) 
        call solve_diis(MdMax,B,Mact,x) 
        Pdiis=0.0d0
        do i=1,Mact
           !Pdiis= Pdiis +x(i)* Ptab(:,:,iact)     ! orthogonal
           Pdiis= Pdiis +x(i)* Ptab(:,:,i)     ! orthogonal
           enddo
        write(*,*)'**** Diis: x ='; call dumpm(x)
        call   ontrans(U,Ui,Pdiis   , scr    ,nnbasis,iPo2Pn) 
        call   formg_shrs (Gn,nnbasis,C,iiat,Nat,C,iiat,Nat ,S,scr,gammat)
        Fn = Hcore +  Gn          
        call   ontrans(U,Ui,Fn, Fo  ,nnbasis, iFn2Fo)
        !rite(*,*)'***** Pdiis:'; call  dumpm(Pdiis)

        ! prepare error vectors for next DIIS iteration
        !scr = matmul(Po,Fo)
        call dgemm('N','N',nnbasis,nnbasis,nnbasis,1.0d0,Po,nnbasis,Fo,nnbasis,0.0d0,scr,nnbasis)
        scr = scr -transpose(scr)
        scferr  =  ddot(nnbasis*nnbasis, scr,1, scr,1)
        !write(*,*) '--- scferr =',scferr 

        Co = Fo
        call   DSYEV( 'Vectors', 'Upper', Nnbasis, Co, Nbasis, eigs, Work, LWORK, INFO )
        !write(*,*)' Eigs:';  if (iprint>0) write(*,1000)(eigs(i),i=1,Nbasis)
        focc=0.d0             ! clear focc
        qel = 0.5d0 *qeltot   ; ! for fermi 
        call fermi(nnbasis,Tel, qel , eigs, efermi,focc)
        focc= 2.0d0*focc      ! since we have restricted-scc
        call setdens(nnbasis,Co,focc,U,Ui,Pn,Po)

        ! current energy 
        En0 =En
       
        scr = Hcore +0.5d0*Gn
        En =  ddot(nnbasis*nnbasis, scr,1, Pn,1)    
        !write(*,*)' Current En, scferr = ', En, scferr
        de =En -En0 
        write(*,1001)  iter, En, dE, scferr

        iter = iter + 1
        
        !---
        ! Pdiis=P0; 
        ! do i=1,Mact ;  Pdiis= Pdiis + x(i)* Ptab(:,:,iact) ;enddo
 

       !diis -Po
       !Po ->Pn
       ! G(P)
       !  F = H +G
       !  Fo 
       ! diag 
       !  fermi 
       !   Po -> Pn
       !   comutt,  err, B

    enddo
 ! no do while


    stop


    Co = Fo
    call   DSYEV( 'Vectors', 'Upper', Nbasis, Co, Nbasis, eigs, Work, LWORK, INFO )
    
 ! write(*,*)  'Eig-max =', eig(idamax(Nbasis,eig,1))

   !write(*,*)' Eigs:';  if (iprint>0) write(*,1000)(eigs(i),i=1,Nbasis)
    focc=0.d0             ! clear focc
    qel = 0.5d0 *qeltot   ; ! for fermi 
    call fermi(nnbasis,Tel, qel , eigs, efermi,focc)
    focc= 2.0d0*focc      ! since we have restricted-scc
    call setdens(nnbasis,Co,focc,U,Ui,Pn,Po)

    write(*,*)'*** Po ='; call dumpm(Po)
    write(*,*)'*** Fo ='; call dumpm(Fo)

    stop

!


    
    




    write(*,*) '**  gammat:';  call  dumpm(gammat)  
    write(*,*) '**  Pn :';  call  dumpm(Pn)  
    write(*,*) '**  Gn :';  call  dumpm(Gn)  
    write(*,*) '**  Fn :';  call  dumpm(Fn)
    write(*,*) '**  Fo :';  call  dumpm(Fo)
    write(*,*) '**  focc:'; call  dumpm(focc)  
    write(*,*) '**  eigs:'; call  dumpm(eigs)
    write(*,*)' *** S ';call  dumpm(S)
    write(*,*)' *** Po (Cmo):';call  dumpm(Po)
    write(*,*)' *** Pn (Cmo):';call  dumpm(Pn)
   

    
  

1000 format(100(' ',E16.8))


end subroutine  scf_diis_shrs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine  solve_diis(N,B,CC) 
subroutine  solve_diis(Nb,B,N,CC) 
  !---------------------------------------------------
  ! this subroutine solves DIIS problem:
  !   http://vergil.chemistry.gatech.edu/notes/diis/node2.html
  ! CPL, Pullay 1982
  ! scusria, JCP, 2012 .....   etc.
  !----
  !   N <-- how many vectors  used for DIIS
  !   [B,  -1 ]     [c       ]      [0 ]
  !   [-1   0 ]  *  [-lambda ]   =  [-1]
  !   Pdiis =  summ_i{ C_i*P_i }
  ! where  B(i,j) = <DP_i| DP_j > 
  !  and   error vectors  DP_i  are  commutators:
  !     DP_i =  [F(Po), Po]
  !  Then  
  !   Fdiis =  sum_i{  CC_i* F_i  } 
  !---
  !  also condition number is estimated:    acond = |A|*|A^-1|
  !   
  !   This is the same as   function      cond(A,1) from  octave
  ! The smaller 'acond' the better !
  !---------------------------------------------------
 use utilmod
  implicit none 
    integer,       intent (in)  :: N, Nb
    !real (kind=8), intent (in)  :: B(N,N)
    real (kind=8), intent (in)  :: B(Nb,Nb)
    real (kind=8), intent (out) :: CC(N)
    integer                     :: nrhs, M, lda, ldx,  info 
    character,  parameter       ::  trans='N'

    real (kind=8), allocatable  :: A(:,:) ,x(:), work(:)
    integer,       allocatable  :: ipiv(:),  iwork(:) 
    real (kind=8),  external    ::  dlange                 ! lapack  natrix norm routine 'norm'=/o, 1, I/ , infinity
    real (kind=8)               ::  anorm, rcond, acond    ! norm of A, rcond=/reciprocal of condition number/, acond=/condition number for A/ 
    character, parameter        ::  norm='1'

    

    
    M= N +1 ; 
    lda = M ; ldx =M  ; nrhs =1     ! nrhs  :  how many vectors on the right-hand-side of    A*x=y

    allocate (A(M,M), x(M), ipiv(M))
    allocate (iwork(M), work(4*M))     ! for condition  number analysis
        
    A(1:N  ,1:N  )   =     B(1:N,1:N) 
    A(1:N  ,  N+1)   =  -1.0d0
    A(  N+1,1:N  )   =  -1.0d0
    A(  N+1 , N+1)   =   0.0d0

    x(1:N)           =   0.0d0
    x(N+1 )          =  -1.0d0

    !write(*,*)' *** DIIS:  A=' ; call dumpm(A) 
  
    !  get the norm of A before LU decomposition
    anorm  =  dlange('1',lda, lda, A,lda,work)

    call dgetrf(M,M,A,lda, ipiv, info) 
    CALL DGECON( '1', lda, A , lda , ANORM, RCOND, WORK, IWORK, INFO ) 
    acond  =  1.0d0/rcond
    !write(*,*)'diis, condition number:', acond

    if  (INFO.EQ.0) THEN
      call dgetrs(trans,M,nrhs,   A,lda,ipiv,x,ldx,info)
    else 
      write(*,*)' problem with DIIS'
      stop
    endif 
    !write(*,*)'x*** x =' ; call dumpm(x) 

   CC = x(1:N)

end 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  solve_diis2(N,B,CC) 
  !---------------------------------------------------
  ! this subroutine solves DIIS problem:
  ! aith constraint  sum_i{c_i*c_i} = 1
  ! which leads to diagonalization of matrix error
  !  B*c =  lambda*c
  !--------
  !   Pdiis =  summ_i{ C_i*P_i }
  ! where  B(i,j) = <DP_i| DP_j > 
  !  Then  
  !   Pdiis =  sum_i{  CC_i* P_i  } 
  !---------------------------------------------------
 use utilmod
  implicit none 
    integer,       intent (in)  :: N
    real (kind=8), intent (in)  :: B(N,N)
    real (kind=8), intent (out) :: CC(N)
    integer                     :: Lwork, info


!    integer                     :: nrhs, M, lda, ldx,  info 
!    character,  parameter       ::  trans='N'


    
    real (kind=8), allocatable  :: A(:,:),  Work(:) , eigs(:)


    lwork = 5*n*n

    allocate (Work(lwork), A(n,n), eigs(n))
        
    A = B ;

   call dsyev( 'Vectors', 'Upper', n,A,n, eigs, Work,lwork,info)
 
    if (info.ne.0) then 
      write(*,*)' problem with DIIS, info=',info
     write(*,*)'*** B  =' ; call dumpm(B)
      stop
    endif 


    write(*,*)'*** errs =' ; call dumpm(eigs) 
   CC= A(:,1)


end 


!####################

 


subroutine scf_diis_atrs(Tel,qeltot,MdMax,Nat,C,iiat,nnbasis,Pn,U,Ui,S,Hcore,Gn,Fn,iprint)
  use timing 
  use  utilmod
  use hamilsccmod
  use density
  use datamod

  !  MdMax <---  max size of DIIS space
  implicit none
!  integer    ::  nnbasis
  integer,       intent (in) :: nnbasis, Nat, iiat(Nat) , iprint
  real (kind=8), intent (in) :: Tel, C(3,Nat), qeltot
  real (kind=8)              :: Pn(nnbasis,nnbasis),Gn(nnbasis,nnbasis), Fn(nnbasis,nnbasis)
  real (kind=8)              :: S(nnbasis,nnbasis), U(nnbasis,nnbasis), Ui(nnbasis,nnbasis) 
  real (kind=8), intent (in) :: Hcore(nnbasis,nnbasis)
  

  integer                    ::  i,j,k , iact, Mact,  iter , MaxScf , MdMax
  real (kind=8)              ::  qtot , qel, efermi, En, En0
  !real (kind=8), allocatable ::  Ptab(:,:,:),  dPtab(:,:,:), B(:,:)     ! save  density form diagon
  real (kind=8), allocatable ::  Ftab(:,:,:),  dFtab(:,:,:), B(:,:)     ! save  density form diagon
  real (kind=8), allocatable ::  PdTab(:,:,:), dPdTab(:,:,:)            !  save DIIS denisties   
  real (kind=8), allocatable ::  x(:), scr(:,:),Pdiis(:,:)              ! diis vector
  real (kind=8), allocatable ::  gammat(:,:),  Fo(:,:),Po(:,:), Co(:,:), eigs(:), focc(:),  work(:) 
  real (kind=8), allocatable ::  dQatom(:), Qatom(:), Q0atom(:), dQshell(:), Qshell(:), Q0shell(:) ,  dQ(:)

! set o/n transformation flags
  integer, parameter  :: iPo2Pn=1, iPn2Po=2, iFo2Fn=3, iFn2Fo=4
  integer             :: lwork, info       
  integer             :: Nel

  real (kind=8), external    :: ddot
  real (kind=8)              ::  scferr, scftol ,dEtol, dE, Encc,Encc0, Ecoul,Erep ,Eel
  logical                    ::  converged

  integer                    :: iread_gamma    !!!!! debug 

  MdMax=10

  write(*,*)'------- scf_diis_atrs  module -------'
   call time_stamp('..entering  scf_diis_atrs: ')
   call start_time(13)

    lwork = 5*nnbasis*nnbasis
    allocate( gammat(Nat,Nat), Fo(nnbasis,nnbasis), Po(nnbasis,nnbasis),Co(nnbasis,nnbasis))
    allocate(scr(nnbasis,nnbasis),B(MdMax,MdMax),x(MdMax))
    !allocate(dPtab(nnbasis,nnbasis,MdMax),Ptab(nnbasis,nnbasis,MdMax))
    allocate(dFtab(nnbasis,nnbasis,MdMax),Ftab(nnbasis,nnbasis,MdMax))
    allocate(eigs(nnbasis), Work(lwork), focc(nnbasis), Pdiis(nnbasis,nnbasis)) 

    allocate(dQatom(nat),Qatom(nat),Q0atom(nat),dQshell(nnbasis),Qshell(nnbasis),Q0shell(nnbasis))
    allocate(dQ(nbasis))

   !
   !  set q0-shell and q0-atom
   !
      iact = 0
      do i=1,Nat
        q0atom(i) = atoms(iiat(i))%fs  + atoms(iiat(i))%fp  +atoms(iiat(i))%fd
        do j=1,lmax(iiat(i))* lmax(iiat(i))                           !  loop over shells
           if (j==1)                then                              !  s-shell
               q0shell(iact +j)   = atoms(iiat(i))%fs
           elseif (j>1 .and. j<5)   then                              !  p-shell
               q0shell(iact +j)   = atoms(iiat(i))%fp/3.0d0
           elseif (j>4 .and. j<10)  then                              !  d-shell
               q0shell(iact +j)   = atoms(iiat(i))%fd/5.0d0
           endif
           !qatom(i) = qatom(i) + qshell(iact + j)
        enddo
        iact =iact + lmax(iiat(i))* lmax(iiat(i))
      enddo

   !
   !  generate gammamatrix,  Mulliken chares
   !  
    call  formgammat_atrs(gammat, C,iiat,Nat,nnbasis, C,iiat,Nat,nnbasis)
    call  Mulliken(nnbasis, nat,Qatom,Qshell,iiat,S,Pn)
    dQatom  = Qatom - Q0atom 
    dQshell = Qshell- Q0shell

    if (iprint > 1)  then 
      write(*,*)'*** Qatom=';  call dumpm(Qatom)
      write(*,*)'*** Q0atom='; call dumpm(Q0atom)
      write(*,*)'*** dQatom='; call dumpm(dQatom)
      write(*,*)'*** S  ='; call dumpm(S) 
      write(*,*)'*** Pn ='; call dumpm(Pn)
    endif 
    call  elenergy(nbasis,nat,Encc0,Ecoul,en,dQatom,gammat,Pn,Hcore) 

    write(*,*)'-------------------------------------------'
    write(*,*)' NCC energy = ',  Encc0

    !stop 
!----------debug start ----------------
!        write(*,*)'Entering  scf_diis_atrs....'
!        read(*,*)IRead_gamma
!        write(*,*)'IRead_gamma = ', iread_gamma
!        if (IRead_gamma==1) then 
!          write(*,*)' reading gamma from input file'  
!          do i=1,Nat
!             read(*,*)( gammat(i,j),j=1,i)
!            do j=1,i
!              gammat(j,i)= gammat(i,j)
!            enddo  
!          enddo
!          write(*,*)' reading charges  from input file' 
!              read(*,*)( focc(j),j=1,nat)
!        else 
!          write(*,*) 'gamma   calculated'
!        endif 
!
!       write(*,*) '***  gammat(read in)= ';   call  dumpm(gammat)
!       write(*,*) '*** charges(read in)= ',(focc(j),j=1,nat)    ! here focc stores dQ 
!       do i=1,Nat
!         focc(i) = focc(i) - (atoms(iiat(i))%fs  + atoms(iiat(i))%fp  +atoms(iiat(i))%fd )
!       enddo
!
!       call  formg_atrs (Gn,nnbasis,C,iiat,Nat,C,iiat,Nat,  S, focc,gammat)
!       write(*,*) '*** dQ: charges: ',(focc(j),j=1,nat)    ! here focc stores dQ 
!       write(*,*) '*** Gn (gammat)= ';   call  dumpm(Gn)
!       !stop
!----------- debug ends -----------------
    !dQatom =  focc   ! here focc was read in

!    call  time_start(15)
!    call  formg_atrs (Gn,nnbasis,C,iiat,Nat,C,iiat,Nat,  S, dQatom, gammat)
!    call  time_stop(15)
!    call  elenergy(nbasis,nat,Encc,Ecoul,en,dQatom,gammat,Pn,Hcore) 
 



  !
  !  current energy 
  !     Pn -->  Qatom -->  Gn -->  ontrans -->  ScfErr=[F,P]  --> Energy  
    call time_start(19)
    call  Mulliken(nnbasis, nat,Qatom,Qshell,iiat,S,Pn)
    call time_stop(19)
    dQatom  = Qatom - Q0atom  ;  dQshell = Qshell- Q0shell
    call time_start(15)
    call  formg_atrs (Gn,nnbasis,C,iiat,Nat,C,iiat,Nat,  S, dQatom, gammat)
    call time_stop(15)

    call  elenergy(nbasis,nat,Encc,Ecoul,en,dQatom,gammat,Pn,Hcore) 
    
  
1001 format('Iter,  En, dEn,  err: ',I6, ' ',F16.8,2('  ',E12.6))


  !
  !  initialize DIIS
  !
    En0= 0.0d0
    !dPtab =0.d00   ;   Ptab  =0.0d0
    B=0.0d0
    scftol=1.0d-8 
    dEtol  =1.e-6
    !MaxSCF =200  
    MaxSCF =20  
    iact =1  
    iter =0 
    Fn =  Hcore +  Gn 

   
!---  now  diag ------ check here!!!!!!!!!!
    call  time_start(18)
    call  ontrans(U,Ui,Pn     , Po     ,nnbasis,iPn2Po)
    call  ontrans(U,Ui,Fn     , Fo     ,nnbasis,iFn2Fo)
    if (iprint>1) then ;  write(*,*)'*** Fo='; call dumpm(Fo,10) ;endif
    call  time_stop(18)


    call  time_start(16)
    !scr= matmul(Po,Fo) 
    call dgemm('N','N',nnbasis,nnbasis,nnbasis,1.0d0,Po,nnbasis,Fo,nnbasis,0.0d0,scr,nnbasis)
    scr = scr -transpose(scr)

!---  here iact=1
    !dFtab(:,:,1) =  scr
    dFtab(:,:,1) =  scr
    ScfErr =  ddot(nnbasis*nnbasis, scr,1, scr,1)
    B(1,1)    = ScfErr
    Ftab(:,:,1) =  Fo 
    call time_stop(16)

    iter =1 
    de =En-En0
    !Pdiis = Po 
    converged =   (SCFerr < SCFtol) .and. (dE <  dEtol)

    ! 1)  diag  Fo ,  get new Po
    ! 2)  form  new Fo and Fo   ,  <--- save Fo
    ! 3)  test the error  scferr=  [Fo,Po}     <--- converged? 
    ! 4)  converged=no   if not then diis:  B  matrix, get Fdiis, <------- save scferr=  [Fo,Po}   
    iter =1  

    write(*,*)'-------------------------------------------'
    call  time_stamp('.. time before  SCC:')
    write(*,1102)
    write(*,*)'Iter,          E(1el),              E(2el),               Escf                dE '    
    write(*,1101)
    write(*,1100)iter,Encc, Ecoul, En,dE  

1100  format(I6,10(' ',F18.10))
1101  format(90('-')) 
1102  format(90('=')) 

    do while  (iter < MaxSCF  .and. .not. converged )
    
    ! update iteration counters
        iter = iter+1 
        Mact = min(MdMax,iter)
        iact = mod(iter,Mact) ;  if (iact.eq.0)  iact=Mact
        !iact = mod(iter,MdMax) 
        !write(*,*)'------------------------------------------'
        !write(*,*)'iter,iact,Mact,i=',iter,iact,mact,i

     
    ! diag  Fo ,  get new Po
        Co = Fo
        !write(*,*)' ... diagonalize: Fo ';  call  dumpm(Fo,10) 
        call   start_time(14)
        call   DSYEV( 'Vectors', 'Upper', Nnbasis, Co, Nbasis, eigs, Work, LWORK, INFO )
        call   stop_time(14)
        !write(*,*)' ... eigs: ';  call  dumpm(eigs) 
        !write(*,*)' Eigs:';  if (iprint>0) write(*,1000)(eigs(i),i=1,Nbasis)
        focc=0.d0             ! clear focc
        qel = 0.5d0 *qeltot   ; ! for fermi 
        call fermi(nnbasis,Tel, qel , eigs, efermi,focc)
        focc= 2.0d0*focc      ! since we have restricted-scc
        call start_time(17)
        call setdens(nnbasis,Co,focc,U,Ui,Pn,Po)
        !write(*,*)' ... form: Po ';  call  dumpm(Po,10) 
        call stop_time(17)

    ! form  new Gn and Fo   ,  <--- save Fo
        call time_start(19)
        call  Mulliken(nnbasis, nat,Qatom,Qshell,iiat,S,Pn)
        call time_stop(19)
        dQatom  = Qatom - Q0atom  ;  dQshell = Qshell- Q0shell
        call  start_time(15)
        call  formg_atrs (Gn,nnbasis,C,iiat,Nat,C,iiat,Nat,  S, dQatom, gammat)
        call  stop_time(15)
        Fn = Hcore + Gn   
        call  start_time(18)
        call   ontrans(U,Ui,Fn, Fo  ,nnbasis, iFn2Fo)
        call  stop_time(18)
        Ftab(:,:,iact)= Fo         ! <--- save  F

    ! get current energy (Mulliken charges are needed)
        En0 = En
        call  elenergy(nbasis,nat,Encc,Ecoul,en,dQatom,gammat,Pn,Hcore) 
        de =En -En0 

    ! test the error  scferr=  [Fo,Po}     <--- save dFtab !converged? 
        call  start_time(16)
        !scr = matmul(Po,Fo)
        call dgemm('N','N',nnbasis,nnbasis,nnbasis,1.0d0,Po,nnbasis,Fo,nnbasis,0.0d0,scr,nnbasis)
        scr = scr -transpose(scr)
        !write(*,*)'*** Po'; call dumpm(Po,10)
        !write(*,*)'*** Fo'; call dumpm(Fo,10)
        !write(*,*)'*** [Po,Fo]'; call dumpm(scr,10)
        dFtab(:,:,iact)=  scr
        scferr =  ddot(nnbasis*nnbasis, dFtab(1,1,iact), 1, dFtab(1,1,iact),1)   
        !write(*,*)'iact, scferr :',iact, scferr
        converged =   (abs(SCFerr) < SCFtol) .and. (abs(dE) <  dEtol) 
        !write(*,*)'iter,scferr,scftol=',iter,scferr,scftol
        !write(*,*)'iter,de,     detol=',iter,de,detol 


    ! converged=no ?  if not then diis:  B  matrix, get Fdiis, <------- save scferr=  [Fo,Po}   
            B(iact,iact) = scferr
       if  (.not.converged)  then 
        do i=1, Mact
          if (i.eq.iact) then 
            !B(iact,iact) = scferr
          else
            scferr =  ddot(nnbasis*nnbasis, dFtab(1,1,i), 1, dFtab(1,1,iact),1)   
            !scr= dFtab(:,:,iact) ;  write(*,*)'***   dFtab(iact)',iact ;  call dumpm(scr,11)
            !scr= dFtab(:,:,i   ) ;  write(*,*)'***   dFtab(i   )',i    ;  call dumpm(scr,11)
            !scferr= iter
            B(i,iact) = scferr 
            B(iact,i) = scferr 
          endif
        enddo
        !write(*,*)' *** B:';  call dumpm(B,Mact)
        !call solve_diis(Mact,B,x) 
        call solve_diis(MdMax,B,Mact,x) 
        !write(*,*)' *** x:';  call dumpm(x,Mact)
   
        !write(*,*)'  ok -1'; call flush(6)
        Fo=0.0d0
        do i=1,Mact
           Fo= Fo + x(i)* Ftab(:,:,i)     ! orthogonal
           !write(*,*)'*** Ftab(i)=',i; scr=Ftab(:,:,i); call dumpm(scr,10)
           enddo
       endif 
       !write(*,*)' *** Fo,diis =';  call dumpm(Fo,10) 
       call  stop_time(16)

    ! print results
      write(*,1100)iter,Encc, Ecoul, En,dE  
       

    enddo   !while

    Nel= int(qel+1.0d-8) 
    write(*,1102)
    write(*,*)'*** Orbital energies(Occupied)  a.u.: ' ;  write(*,1200) (eigs(i),i=1,Nel)
    write(*,*)'*** Orbital energies(Virtual)   a.u.: ' ;  write(*,1200) (eigs(i),i=Nel+1,Nnbasis)

    call  stop_time(13)   
    call  stop_time(10)   

1200 format(5(' ',F16.8))
1201 format(A30,F16.4)


    write(*,1101)
    write(*,*)' NCC   energy = ', Encc0
    write(*,*)' Band  energy = ', Encc
    write(*,*)' El-El energy = ', Ecoul
    write(*,*)' SCC   energy = ', En 
     

    write(*,*)
    write(*,*) '--------------- Timing information --------------'  
    write(*,1201) ' Total time        (sec) :  ',  timetab(10) 
    write(*,1201) '    Total NCC time (sec) :  ',  timetab(2)
    write(*,1201) '        ncc Hamilt (sec) :  ',  timetab(1)
    write(*,1201) '        cholesky   (sec) :  ',  timetab(3)+timetab(4)
    write(*,1201) '            S^(1)  (sec) :  ',  timetab(3)
    write(*,1201) '            U,Ui   (sec) :  ',  timetab(4)
    write(*,1201) '    Total SCC time (sec) :  ',  timetab(13)
    write(*,1201) '        scc hamilt (sec) :  ',  timetab(15)
    write(*,1201) '        diag  time (sec) :  ',  timetab(14)
    write(*,1201) '        set densty (sec) :  ',  timetab(17)
    write(*,1201) '        ontrans    (sec) :  ',  timetab(18)
    write(*,1201) '        Mulliken   (sec) :  ',  timetab(19)
    write(*,1201) '        diis  time (sec) :  ',  timetab(16)
    write(*,*) 
    call  time_stamp('.. time after  SCC:')
    
    stop
    write(*,1102)


    !   write(*,*) '*** gammat= ';   call  dumpm(gammat)



    stop

    
    




    write(*,*) '**  gammat:';  call  dumpm(gammat)  
    write(*,*) '**  Pn :';  call  dumpm(Pn)  
    write(*,*) '**  Gn :';  call  dumpm(Gn)  
    write(*,*) '**  Fn :';  call  dumpm(Fn)
    write(*,*) '**  Fo :';  call  dumpm(Fo)
    write(*,*) '**  focc:'; call  dumpm(focc)  
    write(*,*) '**  eigs:'; call  dumpm(eigs)
    write(*,*)' *** S ';call  dumpm(S)
    write(*,*)' *** Po (Cmo):';call  dumpm(Po)
    write(*,*)' *** Pn (Cmo):';call  dumpm(Pn)
   

    
  

1000 format(100(' ',E16.8))


end subroutine  scf_diis_atrs 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine elenergy(nbasis,nat,iat1,energy,S,dQatom P,Hcore,G,F,gammat)
!subroutine elenergy(nbasis,nat,energy,Encc,Ecoul,dQatom,gammat,P,Hcore)
subroutine elenergy(nbasis,nat,Encc,Ecoul,energy,dQatom,gammat,P,Hcore)
    integer  :: nbasis,  nat
    integer  :: iat1(nat)
    real (kind=8) :: S(nbasis,nbasis), P(nbasis,nbasis), Hcore(nbasis,nbasis), F(nbasis,nbasis)
    real (kind=8) :: gammat(nat,nat),dQatom(nat)
    real (kind=8) ::   Encc,  Ecoul, energy
    real (kind=8),allocatable :: scr(:)
    integer       :: i,j,k
    real (kind=8), external :: ddot

   

    allocate(scr(nat))
    Encc =  ddot(nbasis*nbasis, Hcore,1, P,1) 
    call  dgemv('N',nat,nat, 1.0d0, gammat,nat,dqatom,1, 0.0d0, scr,1)
    Ecoul = 0.50d0  *   ddot(nat,dqatom,1,scr, 1 ) 
    
    energy =  Encc +  Ecoul

!    write(*,100) Encc,  Ecoul, Energy
100 format ('Encc, Ecoul, Etot=',  3(' ',F18.10))


end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  mulliken(Nnbasis,Nnat,Qatom,Qshell,iat1,S,P)
  use datamod

  implicit none
  integer       :: nnbasis, Nnat, i,j,k,ll, iact
  integer       ::  iat1(nnat)
  !real (kind=8) :: Qatom(nnat), Qshell(nnat), S(nnbasis,nnbasis), P(nnbasis,nnbasis)
  real (kind=8) :: Qatom(nnat), Qshell(nnbasis), S(nnbasis,nnbasis), P(nnbasis,nnbasis)
  integer, allocatable ::  MapBs2At(:)
  real (kind=8), external ::  ddot

  !
  !--   mark  which  basis  belong to which atom
  !
     allocate (MapBs2At(Nnbasis))  
     iact=1
     do i=1,Nnat
        ll= lmax(iat1(i))* lmax(iat1(i))
        MapBs2At(iact:iact+ll-1) = i 
        iact = iact +  ll
     enddo
  
  !
  ! set set shell and atomic charges: qshell, qatom
  !
     Qatom= 0.d0
     do i =1,Nnbasis
        qshell(i) = ddot(Nnbasis,S(1,i),1, P(1,i),1)
        j= MapBs2At(i)
        qatom(j) = qatom(j)  +  qshell(i)
     enddo

end  subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
