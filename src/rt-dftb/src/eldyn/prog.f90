program prog
!-------------------------------------------------------------------------------
! 
! N          -(input)  number of basis functions
! S(N,N)     -(input) overlap matrix dim =N*N
! F(N,N)     -(input) Hamiltonian *dt , nonorthogonal basis  
! Pn0(N,N,2) -(input) initial density matrix (real matrix followed by imaginary)
! Pn1(N,N,2) -(output) final density matrix (real  matrix followed by imaginary)
! Iprint     -(input) debug print (=0 means no debug print, otherwise will print parts of matrices)
! ieldyn     -(input) select the algorithm for dynamics:
!            =1: BCH commutator expanson,
!            =2: diagonalization followed by explicit formation of  time-evolution operator
!            =0:  run  both (BCH followed by explicit propagator), this is good for  timing and debugging.
!-------------------------------------------------------------------------------
use  timing
use  utilmod
use  initmd
implicit none
  integer ::  N, iprint,ieldyn
  real*8, allocatable , dimension(:,:,:) ::  Pn0,Pn1 
  real*8, allocatable , dimension(:,:)   ::  S,F


!
!   get basic info
!
 write(*,*)'N ,ieldyn/0-2/, iprint =?'
 read(*,*) N, ieldyn, iprint
 write(*,*) N,ieldyn,  iprint

!
! generate random S, Pn, Fdt 
!   
  allocate(S(n,n), F(n,n), Pn0(n,n,2), Pn1(n,n,2))
  call start_time(2)
  call init(S,Pn0,F,N)   !    F   matrix  is expected to have (F*dt)
  call stop_time(2)

!
! run dynamics
! 
   call eldyn(N,S,F,Pn0,Pn1,Ieldyn,iprint)
 
!
!  dump time statistics
!
  write(*,*)'====================================================='
  write(*,*)' Total     time   (sec) :', timetab(1)
  write(*,*)' Init-SPF  time   (sec) :', timetab(2)
  write(*,*)' chol:Sinv time   (sec) :', timetab(3)
  write(*,*)' chol:U,Ui time   (sec) :', timetab(4)
  write(*,*)' n/o trans time   (sec) :', timetab(5)
  write(*,*)' bch-prp   time   (sec) :', timetab(6)
  write(*,*)' diag-prp  time   (sec) :', timetab(7)
  write(*,*)
  write(*,*) ' Note: total time does not include  Init-SPF'
  write(*,*)'====================================================='

end program 

