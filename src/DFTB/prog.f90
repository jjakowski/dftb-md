program  prog
use timing
use slkomod
use splinemod 
use datamod
use hamilmod 
use hamilsccmod 
use density

!use temp
  implicit none
  character ::  infile*80,  outfile*80, slkofile*128
  integer   ::  Nat ,Nel , MdMax
  logical   ::  ex
  integer   ::  ntype,   istatus, iHomo,iHetero 

  real (kind=8), allocatable, dimension(:,:)  :: C0, C,CC, Frc, Rij
  real (kind=8), allocatable  ::  Hcore(:,:),S(:,:), G(:,:) , Stemp(:,:), Htemp(:,:), gammat(:,:)
  integer,  allocatable  :: iat0(:)

!  integer , allocatable  ::  iat(:), lmax(:)
  real (kind=8), allocatable, dimension(:,:) :: Cn,Co,Pn,Po,Hn,Ho,Gn,Go,U,Ui,Si, Fn,Fo

! set o/n transformation flags
  integer :: iPo2Pn=1, iPn2Po=2, iFo2Fn=3, iFn2Fo=4 
  
!-diagon:
  real (kind=8), allocatable, dimension(:) :: work,eig,focc
  integer ::  lwork, info
  integer :: idamax, iprint=1

!- second order SCC-term
  real (kind=8),  allocatable,  dimension(:) :: qshell, qatom
  real (kind=8) :: ddot , qtot, qel, charge, Tel, efermi=0.0d0 ,  energy

  !integer ::  infile_fd=11
 integer :: istdin=5, istdout=6 
!  integer :: istdin=101, istdout=6 

  integer :: i,j,k  ,IRead_gamma

  call time_stamp(' Starting   dftb. ') 

  infile='IN_dftb'
  inquire(file=infile,exist=ex)
  if(ex)then
    open(istdin,file=infile,status='old')
  else
    write(*,*)'STOP: Cannot find file: ',infile
    stop
  endif

 

! istdout =7 
 write(istdout,*) 'How many different atoms: Ntype, iprint =?'
 read(istdin,*) ntype, iprint
 write(istdout,*) ntype, iprint
 allocate (lmax(ntype),STAT=istatus) ; if (istatus/=0) stop "cannot allocate lmax"
 read(istdin,*)(lmax(i),i = 1,ntype)

!--- read   slkofiles
 iHomo=1
 iHetero=0
 do i=1,ntype
  do j=1,ntype 
    read(istdin,*)slkofile
    call  rdslko(ntype,i,j,slkofile)    !<--- uncomment here!
  enddo
  enddo
  write(*,*)'  Done with reading SLKO files'
!
!-- prepare second dervis table for  Cubic splines (setsktab2) 
  call setsktab2(ntype)                 !<--- uncomment here!

!
!----  read charge
  read(istdin,*)  charge
!
!--- read geometry 

  read(istdin,*) Nat
  allocate (iat0(Nat),C0(3,Nat), STAT = iStatus)
  allocate (iat(Nat) ,C(3,Nat) , Frc(3,Nat), STAT = iStatus)
  if (iStatus /= 0) stop "****** Exiting: Cannot allocate iat,C,Frc ******"  
  read(istdin,*)        !   blank/comment line
  do i=1,Nat
    read(istdin,*) iat0(i), C0(1,i),C0(2,i),C0(3,i)
    enddo
   C0 = C0/ bohr2A
  call start_time(10)

  nbasis = 0
  do i=1,Nat 
    Nbasis =nbasis +  lmax(iat0(i))**2 
    enddo
    write(*,*)'nbasis =', nbasis

   
   !----------------
   ! form S, Hcore

   write(*,*)'iat0 =',iat0
   call  sortatoms(Nat,ntype,C0,iat0,C,iat)
   write(*,*)'iat1 =',iat

   allocate(Hcore(nbasis,nbasis),S(nbasis,nbasis))
   allocate(Htemp(nbasis,nbasis),Stemp(nbasis,nbasis))
   allocate(Fn(nbasis,nbasis))
   Stemp =0.0d0;  Htemp=0.0d0  

   call  start_time(1)
   call  start_time(2)
   call formsh(Htemp,Stemp,C,iat,Nat,Nbasis)
   call  trnsl2i(Stemp,S,mapbas(1,1),nbasis)
   call  trnsl2i(Htemp,Hcore,mapbas(1,1),nbasis)
   call  stop_time(1)

   if (iprint > 1) then 
    write(*,*)'*** trnsl: main S:';  call  dumpf(S,nbasis)
    write(*,*)'*** trnsl: main H:';  call  dumpf(Hcore,nbasis)
   endif 

   write(*,*)'dRGrids =',  dRGrids
   write(*,*)' NGrids =', NGrids

!-  get the number of electrons
   qtot= -charge 
   do i=1,Nat
     j =iat0(i)
     qtot=qtot +  atoms(j)%fs + atoms(j)%fp + atoms(j)%fd
   enddo
   Nel=qtot  +1.0d-8      !  add dQ=1.0d-8  to ensure correct float->integer   converion
   write(*,*) '********* charge, Nel , qtot =',charge,Nel,qtot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   allocate (Si(nbasis,nbasis), U(nbasis,nbasis),  Ui(nbasis,nbasis), Co(nbasis,nbasis), Cn(nbasis,nbasis))
   allocate (Hn(Nbasis,Nbasis), Ho(nbasis,nbasis), Gn(nbasis,nbasis), Po(nbasis,nbasis), Pn(nbasis,nbasis))
   allocate( eig(Nbasis),focc(nbasis), work(5*Nbasis*Nbasis) )
   !allocate (Si(nbasis,nbasis), Hn(Nbasis,Nbasis), Ho(nbasis,nbasis), Po(nbasis,nbasis), Pn(nbasis,nbasis), U(nbasis,nbasis),Ui(nbasis,nbasis),Co(nbasis,nbasis),Cn(nbasis,nbasis),Gn(nbasis,nbasis))
   lwork =5*Nbasis*Nbasis
 
   !
   !   set intial  guess density
   !  cholesky decomposition and  O/N transformation matrices U,Ui
   call dochol(S,Si,U,Ui,Nbasis) ! includes timing inside dochol (split on Sinverse and U,Ui
   Hn=Hcore
   call  ontrans(U,Ui,Hn     , Ho     ,nbasis,iFn2Fo) !;  write(*,*)'***  main Ho:';  call  dumpm(Ho,nbasis)  
   Co = Ho
   call   DSYEV( 'Vectors', 'Upper', Nbasis, Co, Nbasis, eig, Work, LWORK, INFO )
   ! write(*,*)  'Eig-max =', eig(idamax(Nbasis,eig,1))
   
   !write(*,*)' Eigs:';  if (iprint>0) write(*,1000)(eig(i),i=1,Nbasis)

    Tel=100.0d0  ; efermi= -2.0d0  ; focc=0.0d0 
    qel = 0.5d0 *qtot 
    call fermi(nbasis,Tel, qel , eig, efermi,focc)

    write(*,*)'----iter 1------'
    call setdens(nbasis,Co,focc,U,Ui,Pn)
    Pn=  2.0d0*Pn                   ! Ptot for  restricted  SCF  
  
    if (iprint > 1) then 
      write(*,*)' *** Pn (Cmo):';call  dumpm(Pn)
      write(*,*)' efermi =', efermi
      write(*,*)'focc =', focc
      write(*,*)'eig =', eig


      write(*,*)'****  Hn ';  call  dumpm(Hn)
      write(*,*)'****  Ho ';  call  dumpm(Ho)
    endif

    !allocate( gammat(nat,nat))
    !call formgammat_atrs(gammat, C0,iat0,Nat,nbasis, C0,iat0,Nat,nbasis)
   call stop_time(2)
!    write(*,*) ' Total NCC time (sec) :  ',  timetab(2)
!    write(*,111) ' Total NCC time (sec) :  ',  timetab(2)
!    write(*,112) ' Total NCC time (sec) :  ',  timetab(2)
!111 format(A30,F20.6)
!112 format(A,F20.6)

!-----scf here

     !qtot  <-- total number of electrons 
     MdMax=4 
     !call  scf_diis_shrs(Tel,qtot,MdMax,Nat,C0,iat0,nbasis,Pn,U,Ui,S,Hn, Gn, Fn, iprint)
     tel=1.0e5
     tel=1.0e6
     tel=1.0d0
     call  scf_diis_atrs(Tel,qtot,MdMax,Nat,C0,iat0,nbasis,Pn,U,Ui,S,Hn, Gn, Fn, iprint)

    
     !call  ElEnergy(nbasis,nat,energy,Pn,Sn,Hn,Gn,gammat) 
     


    
    stop



    allocate( gammat(nbasis,nbasis))
    call formgammat_shrs(gammat, C0,iat0,Nat,nbasis, C0,iat0,Nat,nbasis)

      read(*,*)IRead_gamma
       if (IRead_gamma==1) then 
         write(*,*)' reading gamma from input file'  
         do i=1,Nat
             read(*,*)( gammat(i,j),j=1,i)
           do j=1,i
              gammat(j,i)= gammat(i,j)
           enddo  
         enddo
       else 
         write(*,*) 'gamma   calculated'
        endif 
      
    !write(*,*) call dumpm(gammat)      
    stop
    !-- iter-1
    write(*,*)'----iter 1------'
    write(*,*)'focc =', focc
    call setdens(nbasis,Co,focc,U,Ui,Pn)
    write(*,*)' *** Pn (Cmo):';call  dumpm(Pn)
    call  formg_shrs   (Gn,nbasis,C0,iat0,Nat,C0,iat0,Nat ,S,Pn,gammat)
    Hn =  Hcore + Gn 
    call  ontrans(U,Ui,Hn     , Ho     ,nbasis,iFn2Fo)
    Co = Ho
    call   DSYEV( 'Vectors', 'Upper', Nbasis, Co, Nbasis, eig, Work, LWORK, INFO )
   

    !--- iter-2
    write(*,*)'----iter 2------'

    focc=0.d0             ! clear focc
    qel = 0.5d0 *qtot   ; ! for fermi 
    call fermi(nbasis,Tel, qel , eig, efermi,focc)
    focc= 2.0d0*focc      ! since we have restricted-scc
    call setdens(nbasis,Co,focc,U,Ui,Pn)
    write(*,*)' *** Pn (Cmo):';call  dumpm(Pn)
  
 
   


    write(*,*)'****  Gammat(nat)';  call  dumpm(gammat)  
    write(*,*)'****  Hcore';  call  dumpm(Hcore)
    write(*,*)'****  Gn ';  call  dumpm(Gn)
    write(*,*)'****  Hn ';  call  dumpm(Hn)
    write(*,*)'****  Ho ';  call  dumpm(Ho)


    




    stop 
    deallocate(gammat)

    write(*,*)'****  Gammat(nbasis)';  call  dumpm(gammat)  





   









   stop

 write(*,*)' *** Pn (Cmo):';call  dumpm(Pn)
    
    call setdens(nbasis,Co,focc,U,Ui,Pn)
!!!!!  different Pn:
   write(*,*)'****** Pn from focc (Q-shrs)' 
   read(*,*)(focc(i),i=1,nbasis)
   write(*,*)'in; focc=',focc
   call  setdens(nbasis,focc,Pn)
    write(*,*)' *** Pn(focc):';call  dumpm(Pn)

!!---
   write(*,*)'****** Pn from atoms (Q-atrs)'
   allocate(qatom(nat))
   read(*,*)(qatom(i),i=1,nat)
   write(*,*)'in; qatom=',qatom 
  
   call  setdens_qatom(nbasis,nat,qatom,iat0,Pn)
    write(*,*)' *** Pn(qatoms):';call  dumpm(Pn)

    


 
  stop
!
!  O/N transformation Pn-Po
!      syntax:   ontrans (U,Ui, Ain, Aout,Nbasis,ifon)
!
   call start_time(5)
   call  ontrans(U,Ui,S      , Ho     ,nbasis,iFn2Fo) !;  write(*,*)'***  main So=Li*S*Ui:';  call  dumpm(Ho,nbasis) ; Ho=0.0d0
   call  ontrans(U,Ui,Hn     , Ho     ,nbasis,iFn2Fo) !;  write(*,*)'***  main Ho:';  call  dumpm(Ho,nbasis)  
   call  ontrans(U,Ui,Pn(1,1), Po(1,1),nbasis,iPn2Po)
   call  ontrans(U,Ui,Pn(1,1), Po(1,1),nbasis,iPn2Po)
   call stop_time(5)

!--- diagonalize Ho
   Hn=Hcore
   lwork =5*Nbasis*Nbasis

  allocate( eig(Nbasis),focc(nbasis), work(5*Nbasis*Nbasis) )
   Co = Ho
  write(*,*)'***  main Ho:';  call  dumpm(Co,nbasis) 

write(*,*)' now diag'
   call   DSYEV( 'Vectors', 'Upper', Nbasis, Co, Nbasis, eig, Work, LWORK, INFO )
    write(*,*)  'Eig-max =', eig(idamax(Nbasis,eig,1))
    if (iprint>0) write(*,1000)(eig(i),i=1,Nbasis)

  
!  prepare density matrix
   focc=0.0d0
   !do i=1,Nel 
   do i=1,6
       focc(i) = 1.0d0
       enddo
   Pn=0.0d0  ; Po=0.0d0
   call  dcopy(Nel*Nbasis, Co,1, Pn,1) ! here Pn is scratch
   do i=1,Nel
     call  dscal(Nbasis,focc(i),Pn,1)  
     enddo
   write(*,*)'***  main focc*Co:';  call  dumpm(Pn,nbasis)
   call dgemm('N','T',Nbasis,Nbasis,Nel,1.0d0,Co,Nbasis,Pn, Nbasis,0.0d0,Po,Nbasis) 
   write(*,*)'***  main Po=Co*trsp(focc*Co):';  call  dumpm(Po,nbasis)
!-- check trace
qtot=0.0d0
do i=1,nbasis
   qtot=qtot +Po(i,i)
   enddo
write(*,*)'***trace Po=',qtot   

  
!  --non-orthogonal   density
    call  ontrans(U,Ui,Po(1,1), Pn(1,1),nbasis,iPo2Pn)
 write(*,*)'**** main  Pn:';  call dumpm(Pn,nbasis)

! Mulliken population
   allocate(  qatom(Nat)  , qshell(nbasis))
   qtot=0.0d0
   do i =1,Nbasis
     qshell(i) = ddot(Nbasis,S(1,i),1, Pn(1,i),1) 
     qtot= qtot +qshell(i)
     enddo
    write(*,*)'*****  charges:' ; call dumpm(qshell,nbasis)
    write(*,*)'*****  Q-tot:  ', qtot
!---    
qtot=0.0d0
do i=1,nbasis
   qtot=qtot +Pn(i,i)
   enddo
write(*,*)'***trace Pn=',qtot   

!!!!!  different Po:
   read(*,*)(focc(i),i=1,nbasis)
   write(*,*)'******different version Po:'
   write(*,*)'in; focc=',focc
   Pn=0.0d0  ; Po=0.0d0 
   call  dcopy(nbasis,focc,1,Pn,nbasis+1)
   call  ontrans(U,Ui,Pn(1,1), Po(1,1),nbasis,iPn2Po)
  write(*,*)'**** main-new  Pn:';  call dumpm(Pn,nbasis)
  write(*,*)'**** main-new  Po:';  call dumpm(Po,nbasis)


do i=1,nat
   j = iat0(i)
   !fs= atoms(j)%fs
   !fp= atoms(j)%fp
   !fd= atoms(j)%fd
write(*,*)i, j, atoms(j)%fs, atoms(j)%fp, atoms(j)%fd 
enddo
write(*,*)'chargei,Nel =',charge,Nel
!-------

 allocate(gammat(nbasis,nbasis))
 gammat=0.0d0
!call  setdens(Nat,nbasis,nel,1,iat0,Pn,Po,focc,Hn,U,Ui)

!call  formg(G,nbasis,C,iat0,Nat,C,iat,Nat,S,gammat) 
!write(*,*) '**** gammat'; call dumpm(gammat)

!-- shel resolution G (JJ versionthat utilizes all hubbards)
call  formg_shrs(Gn,nbasis,C0,iat0,Nat,C0,iat0,Nat,S,Pn) 

! --- atomic resolution G
allocate(Go(nbasis,nbasis))
call  formg_atrs(Go,nbasis,C0,iat0,Nat,C0,iat0,Nat,S,Pn) 

write(*,*)'****  J.J. version G:'; call dumpm(Gn)
write(*,*)'**** M.Els version G:'; call dumpm(Go)
write(*,*)'**** JJ-M.Els version dG:'; call dumpm(Gn-Go)


write(*,*)' *** Uhubb:'
j=1; write(*,1002)j,j , atoms(j)%Us, atoms(j)%Up, atoms(j)%Ud 
j=2; write(*,1002)j,j , atoms(j)%Us, atoms(j)%Up, atoms(j)%Ud 
write(*,*)' *** occ:'
j=1; write(*,1002)j,j , atoms(j)%fs, atoms(j)%fp, atoms(j)%fd 
j=2; write(*,1002)j,j , atoms(j)%fs, atoms(j)%fp, atoms(j)%fd 


1000 format(100(' ',E16.8))
1001 format(100(' ',E12.6))
1002 format(2(' ',I6),100(' ',E16.8))  

  
end

