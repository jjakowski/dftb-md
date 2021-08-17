module  datamod
 
    integer , allocatable  ::  iat(:), lmax(:),MapAt0(:,:), MapAt1(:,:),MapBas(:,:)
    public ::  MapAt0, MapAt1,MapBas

!-- skdata *spl file contents:
! line 1: nGridPoints ,  gridDist
! line 2:
!   Ed, Ep, Es  : on site energies
!   SPE         : spin polarization error
!   Ud, Up, Us  : Hubbard values 
!   fd, fp, fs  : occupations for neutral (guess)
! line 3:
!     mass :  in amu (ignored for hetero pair)
!     c2-c9,rcut : polynomial coeff. and cutoff  radius for Erep. Can be  zero when splines are used
!                   polynomial is    	f(r) = SUM_k  [ C_i* (rcut-r)^k ] 
!     d1-d10   :  dummy variables (not used) 
!      real (kind=8)  ::  mass,c2,c3,c4,c5,c6,c7,c8,c9,rcut,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10
!line  4++: 
!     Hdd0 Hdd1 Hdd2 Hpd0 Hpd1 Hpp0 Hpp1 Hsd0 Hsp0 Hss0 Sdd0 Sdd1 Sdd2 Spd0 Spd1 Spp0 Spp1 Ssd0 Ssp0 Sss0      
!     0=sigma, 1=pi, 2=delta
    integer, save  ::    lensktab  = 20
    integer, save  ::   lenatomtab = 9    
     
    integer, allocatable, dimension(:,:)       :: NGrids 
    real (kind=8), allocatable, dimension(:,:) :: dRGrids
    real (kind=8), allocatable, dimension(:,:,:,:) :: sktab,  sktab2    !dim(20,ngrids,ntype,ntype)
    !real (kind=8), allocatable, dimension(:,:,:)   :: atomtab           !order (Ed, Ep, Es,  Ud, Up, Us  ,  fd, fp, fs)
    !!! http://courses.physics.illinois.edu/phys466/comp_info/derived.html
    type  onsite 
      sequence
      real (kind=8) ::  Ed,Ep,Es,Ud,Up,Us,fd,fp,fs,mass
    end type onsite 
    type (onsite), allocatable, dimension(:)     ::  atoms


    real (kind=8),save  :: bohr2A =  0.529177249d0
    real (kind=8),save  :: A2bohr=  1.0d0/0.529177249d0
    integer ::  nbasis
!real (kind=8),save  :: A2bohr=  1.0d0/A2bohr
!public   :: sktab,sktab2 



end module
