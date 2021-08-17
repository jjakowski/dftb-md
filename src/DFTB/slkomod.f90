module  slkomod
use splinemod

use datamod

!sktab, sktab2:  parameters  and second derivs from  cubic splines

!integer, allocatable, dimension(:,:)       :: NGrids 
!real (kind=8), allocatable, dimension(:,:) :: tab,  dRGrids
!real (kind=8), allocatable, dimension(:,:,:,:) :: sktab,  sktab2
!public   ::tab, sktab,sktab2 

!20,ngrid,ntype,ntype

contains 
  subroutine rdslko(ntype,itype,jtype,slkofile) 
    implicit none 
    integer, intent(in)   :: ntype ,itype,jtype
    character*128, intent(in)  :: slkofile
    logical :: ex
    integer :: i,j ,k ,iHomo !Ihomo=1:  Homonuclear case, else heteronuclear!
 
    !-- skdata
    ! line 1:
          integer ::  nGridPoints
          real (kind=8)  ::  gridDist
    ! line 2:
    !   Ed, Ep, Es  : on site energies
    !   SPE         : spin polarization error
    !   Ud, Up, Us  : Hubbard values 
    !   fd, fp, fs  : occupations for neutral (guess)
          real (kind=8)  ::Ed, Ep, Es, SPE, Ud, Up, Us, fd, fp, fs
    ! line 3:
    !     mass :  in amu (ignored for hetero pair)
    !     c2-c9,rcut : polynomial coeff. and cutoff  radius for Erep. Can be  zero when splines are used
    !                   polynomial is    	f(r) = SUM_k  [ C_i* (rcut-r)^k ] 
    !     d1-d10   :  dummy variables (not used) 
          real (kind=8)  ::  mass,c2,c3,c4,c5,c6,c7,c8,c9,rcut,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10
    !line  4++: 
    !     Hdd0 Hdd1 Hdd2 Hpd0 Hpd1 Hpp0 Hpp1 Hsd0 Hsp0 Hss0 Sdd0 Sdd1 Sdd2 Spd0 Spd1 Spp0 Spp1 Ssd0 Ssp0 Sss0      
    !     0=sigma, 1=pi, 2=delta
     real (kind=8), allocatable, dimension(:,:) :: skbuf  !sktab
   

     write(*,*)' Reading file:  ',slkofile
     inquire(file=slkofile,exist=ex)
     if(ex)then
     open(12,file=slkofile,status='old')
     else
       write(*,*)'STOP: Cannot find file: ',slkofile  ; stop
     endif
     rewind 12
 
     read(12,*) gridDist,nGridPoints
     if (itype==jtype) then   
        read(12,*) Ed, Ep, Es, SPE, Ud, Up, Us, fd, fp, fs
        if (.not.allocated(atoms)) allocate(atoms(ntype))
        atoms(itype)%Ed= Ed 
        atoms(itype)%Ep= Ep 
        atoms(itype)%Es= Es 
        atoms(itype)%Ud= Ud 
        atoms(itype)%Up= Up 
        atoms(itype)%Us= Us 
        atoms(itype)%fd= fd 
        atoms(itype)%fp= fp 
        atoms(itype)%fs= fs 
     endif  
     read(12,*) mass,c2,c3,c4,c5,c6,c7,c8,c9,rcut
 
     allocate(skbuf(lensktab,nGridPoints))
     do i=1,nGridpoints-1 
       read(12,*) (skbuf(j,i),j=1,20)   
       enddo
      close(12)

     if  (.not.allocated(sktab)) allocate(sktab(20,0:nGridPoints,ntype,ntype))
     call  resizeSK(NGridPoints+2,ntype,sktab)
     sktab(:,1:NGridPoints,itype,jtype) =  skbuf
     deallocate(skbuf)
     ! Hdd0 Hdd1 Hdd2 Hpd0 Hpd1 Hpp0 Hpp1 Hsd0 Hsp0 Hss0 Sdd0 Sdd1 Sdd2 Spd0 Spd1 Spp0 Spp1 Ssd0 Ssp0 Sss0 
     if (itype==jtype)  then 
        sktab( 1,0,itype,itype) = Ed    ; sktab(11,0,itype,itype) = 1.0d0 ;  
        sktab( 2,0,itype,itype) = Ed    ; sktab(12,0,itype,itype) = 1.0d0 ;
        sktab( 3,0,itype,itype) = Ed    ; sktab(13,0,itype,itype) = 1.0d0 ;
        sktab( 4,0,itype,itype) = 0.0d0 ; sktab(14,0,itype,itype) = 0.0d0 ;
        sktab( 5,0,itype,itype) = 0.0d0 ; sktab(15,0,itype,itype) = 0.0d0 ;
        sktab( 6,0,itype,itype) = Ep    ; sktab(16,0,itype,itype) = 1.0d0 ;
        sktab( 7,0,itype,itype) = Ep    ; sktab(17,0,itype,itype) = 1.0d0 ;
        sktab( 8,0,itype,itype) = 0.0d0 ; sktab(18,0,itype,itype) = 0.0d0 ;
        sktab( 9,0,itype,itype) = 0.0d0 ; sktab(19,0,itype,itype) = 0.0d0 ;
        sktab(10,0,itype,itype) = Es    ; sktab(20,0,itype,itype) = 1.0d0 ;
       endif
   
     if  (.not.allocated(NGrids)) allocate(NGrids(ntype,ntype)) 
     if  (.not.allocated(dRGrids)) allocate(dRGrids(ntype,ntype)) 
     NGrids(itype,jtype)  =  NGridPoints
     dRGrids(itype,jtype) =   gridDist 
    
  end subroutine 

  !@
  !@c===========================================================================
  !@  zero-th line in sktab is always  =0
  !@
  subroutine resizeSK(NGridPoints,ntype,sktab)
    implicit none 
    real (kind=8) , allocatable , dimension(:,:,:,:) :: temp,sktab
    integer  ::  NGridOld  , low, up
    integer, intent(in) :: NGridPoints,ntype
 
    if  (.not.allocated(sktab))  then 
        allocate(sktab(20,0:nGridPoints,ntype,ntype))
        sktab=0.0d0
        endif
   
    if  (NGridPoints+1.gt.size(sktab,2))  then   ! here assumed  that  ntype  is OK! and only NgridPoints need to be increased 
        allocate (temp(20,0:NGridPoints,ntype,ntype))
        temp =0.0d0
        low =lbound(sktab,2)
        up  =ubound(sktab,2)
        temp(:,low:up,:,:)  =sktab          ! might not  work if ntype mismatch 
        deallocate(sktab)  
        allocate (sktab(20,0:NGridPoints,ntype,ntype))       
        sktab =temp 
        deallocate (temp)
        endif
  end subroutine 

  !@
  !@c===========================================================================
  !@ set second derivatives for cubic interpolation of SK parameters
  !@ 
  subroutine  setsktab2(ntype) 
    implicit none 
    integer, intent(in) ::ntype
    integer  :: i,j,k ,Ngridmax,NGridPoints ,ipar
    ! global:  sktab,sktab2, Ngrids, dRGrids
    real (kind=8), allocatable, dimension(:)  ::  y2tab, ytab
    real (kind=8) ::  dR    

    Ngridmax= size(sktab,2)
    if  (.not.allocated(sktab2 )) allocate(sktab2(20,0:Ngridmax,ntype,ntype) )
    sktab2=0.0d0

    do i=1,Ntype 
      do j=1,Ntype 
        NGridPoints = NGrids(i,j)
        dR = drgrids(i,j) 

        allocate(ytab(NGridPoints), y2tab(NGridPoints) ) 
        do ipar=1,20                          ! loop over H,S parameters
          ytab = sktab(ipar,:,i,j)            ! copy column
          call gety2tab(NGridPoints, dR, ytab,y2tab)
          sktab2(ipar,0            ,i,j)  = 0.0d0    !  line #0 contains onsite  energies Ed,Ep,Es 
          sktab2(ipar,1:NGridPoints,i,j)  = y2tab 
        enddo
        deallocate (ytab,y2tab)
      enddo 
    enddo
  end subroutine

   

end module slkomod
