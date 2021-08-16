module hamilmod 
use  datamod
use  utilmod
use  splinemod

contains 
  subroutine formsh(H,S,C,iat,Nat,Nbasis) 
    implicit none 
    integer, intent(in) :: Nat, Nbasis
    integer, intent(in) ::  iat(Nat)
    real (kind=8)       :: C(3,Nat), S(Nbasis,Nbasis), H(Nbasis,Nbasis)
    real (kind=8),allocatable       :: Stemp(:,:), Htemp(:,:) 
    integer  :: NatS, NatP, NatD, NatF, i,j,k 
    integer  :: NbasS, NbasP, NbasD, NbasF
    integer  :: IshS, IshP, IshD, IshF, Iend
    integer  :: IatS, IatP, IatD, IatF 
    !--
    integer  :: Nat1, Nat2
    integer, allocatable :: iat1(:) , iat2(:)
    real (kind=8), allocatable ::  C1(:,:), C2(:,:)
    !---
    Nat1= Nat ;  Nat2= Nat
    allocate(iat1(Nat1), iat2(Nat2), C1(3,Nat1), C2(3,Nat2))
    iat1 = iat ; C1= C 
    iat2 = iat ; C2= C 
    
    
    !---
    write(*,*)'........... entering formsh:'    
    ! 
    !write(*,*)'nbasis,Nat,C = ',nbasis,C
    !write(*,*)  '*** S = '; call dumpm(S,6)

    NatS=0 ; NatP=0 ; NatD =0  ;  NatF =0
    do i=1,Nat
      if     (lmax(iat(i))==1) then
         NatS = NatS+1
      elseif (lmax(iat(i))==2) then      
         NatP = NatP +1  
      elseif (lmax(iat(i))==3) then 
         NatD = NatD +1
      elseif (lmax(iat(i))==4) then 
         NatF = NatF +1
      endif
    enddo
    !write(*,*)'Nat:S,P,D,F,tot:',NatS,NatP,NatD,NatF, Nat 
    NbasS =    NatS+NatP+NatD +NatF
    NbasP = 3*(NatP+NatD+NatF)
    NbasD = 5*(NatD+NatF) 
    NbasF = 7*NatF

! matrix ranges:
    ishS   = 1  
    ishP   = ishS + NbasS 
    ishD   = ishP + NbasP
    ishF   = ishD + NbasD 
    iend   = ishF + NbasF -1
!  atoms range
    iatS   =  1
    iatP   =  iatS + NatS 
    iatD   =  iatP + NatP 
    iatF   =  iatD + NatD 
    
  

!
! set block   ( S x S )  of  S and Hcore
!
    Nat1= Nat  ;   Nat2= Nat
    call  set_block_SS(H,S,nbasis,C1,iat1,Nat1,C2,iat2,Nat2) 
    !write(*,*)' ***  Hcore: ss1-block =' ; call dumpf(H,10)
    !write(*,*)' ***  S:     ss1-block =' ; call dumpf(S,10)


!
!  check if block P exists and form  the S,H       
!  set block  ( S x P ) 
!
    if (NatP+NatD+NatF >0 )  then 
      allocate(Stemp(NbasS,NbasP) ,Htemp(NbasS,NbasP))
      Stemp =0.0d0    ;    Htemp =0.0d0;   
      Nat1 = Nat                   ! all atoms have S shell
      Nat2 = NatP +NatD +NatF      ! all atoms with P shell
      
      call  set_block_SP(Htemp,Stemp,NbasS,C1,iat1,Nat1,C2(1,iatP),iat2(iatP),Nat2)

      H(IshP:IshD-1,  IshS:IshP-1)  = transpose(Htemp)   
      H(IshS:IshP-1,  IshP:IshD-1)  =           Htemp    
      S(IshP:IshD-1,  IshS:IshP-1)  = transpose(Stemp)   
      S(IshS:IshP-1,  IshP:IshD-1)  =           Stemp    
      !  write(*,*)' ***  Hcore: sp1-block =' ; call dumpf(Htemp,12)
      !  write(*,*)' ***  S:     sp1-block =' ; call dumpf(Stemp,12)
      deallocate(Stemp, Htemp)
      endif



!
!  check if block P exists and form  the S,H       
!  set block  ( P x P ) 
!
    if (NatP+NatD+NatF >0 )  then 
      allocate(Stemp(NbasP,NbasP) ,Htemp(NbasP,NbasP))
      Stemp =0.0d0    ;    Htemp =0.0d0
      Nat1 = NatP +NatD +NatF 
      Nat2 = NatP +NatD +NatF 

      call  set_block_PP(Htemp,Stemp,NbasP,C1(1,iatP),iat1(iatP),Nat1,C2(1,iatP),iat2(iatP),Nat2)     ! (P x P)

      H(IshP:IshD-1,  IshP:IshD-1)  =           Htemp    
      S(IshP:IshD-1,  IshP:IshD-1)  =           Stemp    
      !write(*,*)' ***  Hcore: pp-block =' ; call dumpf(Htemp,6)
      !write(*,*)' ***  S:     pp-block =' ; call dumpf(Stemp,6)
      deallocate(Stemp, Htemp)
      endif    


!
!   checj if D block exists and fomr  S,H
!   set block (S x D)
!
    if (NatD+NatF > 0 ) then 
      allocate(Stemp(NbasS,NbasD), Htemp(NbasS, NbasD))
      !allocate(Stemp(Nbasis,Nbasis), Htemp(Nbasis, Nbasis))
      Stemp =0.0d0   ;   Htemp = 0.0d0
      Nat1 =     Nat           ! all atoms
      Nat2 =     NatD +NatF    ! all atoms with D shell

      call  set_block_SD(Htemp,Stemp,NbasS,C1(1,iatS),iat1(iatS),Nat1,C2(1,iatD),iat2(iatD),Nat2)     ! (P x P)

      H(IshS:IshP-1,  IshD:IshF-1) = Htemp  
      S(IshS:IshP-1,  IshD:IshF-1) = Stemp  

      H(IshD:IshF-1,  IshS:IshP-1) = transpose(Htemp)  
      S(IshD:IshF-1,  IshS:IshP-1) = transpose(Stemp)  
     
      !write(*,*)' ***  Hcore: sd-block =' ; call dumpf(Htemp,5)
      !write(*,*)' ***  S:     sd-block =' ; call dumpf(Stemp,5)
      deallocate(Stemp, Htemp)
      endif
     
!
!   checj if D block exists and fomr  S,H
!   set block (P x D)
!
   if  (NatD+NatF > 0 ) then 
      allocate(Stemp(NbasP,NbasD), Htemp(NbasP, NbasD))
      !allocate(Stemp(Nbasis,Nbasis), Htemp(Nbasis, Nbasis))
      Stemp =0.0d0   ;   Htemp = 0.0d0
      Nat1 =     NatP+NatD+NatF          ! all atoms with P shell 
      Nat2 =     NatD +NatF              ! all atoms with D shell

      call  set_block_PD(Htemp,Stemp,NbasP,C1(1,iatP),iat1(iatP),Nat1,C2(1,iatD),iat2(iatD),Nat2)     ! (P x P)

      H(IshP:IshD-1,  IshD:IshF-1) = Htemp  
      S(IshP:IshD-1,  IshD:IshF-1) = Stemp  

      H(IshD:IshF-1,  IshP:IshD-1) = transpose(Htemp)  
      S(IshD:IshF-1,  IshP:IshD-1) = transpose(Stemp)  
     
      !write(*,*)' ***  Hcore: pd-block =' ; call dumpf(Htemp,5)
      !write(*,*)' ***  S:     pd-block =' ; call dumpf(Stemp,5)
      deallocate(Stemp, Htemp)
      endif

!
!   checj if D block exists and fomr  S,H
!   set block (D x D)
!
   if  (NatD+NatF > 0 ) then 
      allocate(Stemp(NbasD,NbasD), Htemp(NbasD, NbasD))
      !allocate(Stemp(Nbasis,Nbasis), Htemp(Nbasis, Nbasis))
      Stemp =0.0d0   ;   Htemp = 0.0d0
      Nat1 =     NatD +NatF              ! all atoms with D shell 
      Nat2 =     NatD +NatF              ! all atoms with D shell

      call  set_block_DD(Htemp,Stemp,NbasD,C1(1,iatD),iat1(iatD),Nat1,C2(1,iatD),iat2(iatD),Nat2)     ! (D x D)

      H(IshD:IshF-1,  IshD:IshF-1) = Htemp  
      S(IshD:IshF-1,  IshD:IshF-1) = Stemp  

     
      !write(*,*)' ***  Hcore: dd-block =' ; call dumpf(Htemp,10)
      !write(*,*)' ***  S:     dd-block =' ; call dumpf(Stemp,10)
      deallocate(Stemp,Htemp)
      endif
     
     

  end subroutine

 !@
 !@
 !@
  subroutine set_block_SS(H,S,LDA,C1,iat1,Nat1,C2,iat2,Nat2)   !sktab, sktab2
    implicit none 
    integer        :: LDA 
    real (kind=8)  :: H(LDA,*), S(LDA,*)  
    real (kind=8)  :: C1(3,*) ,  C2(3,*)
    integer        :: Iat1(*), Iat2(*), Nat1, Nat2
    integer             :: i,j,k,j0,itype,jtype, nGridPoints
    real (kind=8)       ::  x,y,z,hh,R, R2,dR,  cosl,cosm,cosn,xx
    real (kind=8)       :: VssSig0, VssSig1, d2VssSig0, d2VssSig1, Vss
    real (kind=8)       :: SssSig0, SssSig1, d2SssSig0, d2SssSig1, Sss
    real (kind=8) :: x1,y1,z1, x2,y2,z2 

    hh=1.0d-20   ! infinitesimal number for onsite  data (added to avoid dividing by R=zero, needed for cosines) 
     
    do i=1,Nat1
      do j=1,Nat2 
        
        itype =iat1(i) 
        jtype =iat2(j) 
        NGridPoints = NGrids(itype,jtype)  
        dR = dRGrids(itype,jtype)
        x1 =C1(1,i)  ;   y1 =C1(2,i)  ;   z1 =C1(3,i)  ;   
        x2 =C2(1,j)  ;   y2 =C2(2,j)  ;   z2 =C2(3,j)  ;   
        x = x2 -x1 + hh 
        y = y2 -y1  
        z = z2 -z1  
        R2   = x*x +y*y +z*z 
        R    =  sqrt(R2)
        j0= R/dR 
        xx=  R -j0*dR
        j0= min(j0,NgridPoints+1) 
        !write(*,*) ' set_block_SS: Ngridpoints, j0 =',  NgridPoints, j0
        VssSig0   =  sktab(10,j0  ,itype,jtype)
        VssSig1   =  sktab(10,j0+1,itype,jtype)
        d2VssSig0 = sktab2(10,j0  ,itype,jtype)
        d2VssSig1 = sktab2(10,j0+1,itype,jtype)
        call cubics(xx,Vss,dR,VssSig0,d2VssSig0,VssSig1,d2VssSig1)

        SssSig0   =  sktab(20,j0  ,itype,jtype)
        SssSig1   =  sktab(20,j0+1,itype,jtype)
        d2SssSig0 = sktab2(20,j0  ,itype,jtype)
        d2SssSig1 = sktab2(20,j0+1,itype,jtype)
        call cubics(xx,Sss,dR,SssSig0,d2SssSig0,SssSig1,d2SssSig1)
        H(i,j) = Vss 
        S(i,j) = Sss 
      enddo 
    enddo 


  end subroutine set_block_ss
 !@
 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !@
  subroutine set_block_SP(H,S,LDA,C1,iat1,Nat1,C2,iat2,Nat2) 
    implicit none
    integer        :: LDA 
    real (kind=8)  :: H(LDA,*), S(LDA,*)  
    real (kind=8)  :: C1(3,*) ,  C2(3,*)
    real (kind=8)  :: x1,y1,z1, x2,y2,z2 
    integer        :: Iat1(*), Iat2(*), Nat1, Nat2
    integer             :: i,j,k,j0,itype,jtype
    !integer             :: NatS, NatP, NatD
    real (kind=8)       :: x,y,z,hh,R, Ri, R2,dR,  cosl,cosm,cosn,xx
    real (kind=8)       :: VspSig0, VspSig1, d2VspSig0, d2VspSig1, Vsp
    real (kind=8)       :: SspSig0, SspSig1, d2SspSig0, d2SspSig1, Ssp
    integer  :: iatS, iat_s0,iat_s1, iatP, iat_P0,iat_P1,iblS,iblP
    integer  ::  nGridPoints

    hh=1.0d-20   ! infinitesimal number for onsite  data (added to avoid dividing by R=zero, needed for cosines) 
     
    do i=1,Nat1
      do j=1,Nat2 
        itype =iat1(i) 
        jtype =iat2(j) 
        dR = dRGrids(itype,jtype)
        NGridPoints = NGrids(itype,jtype)  
        x1 =C1(1,i)  ;   y1 =C1(2,i)  ;   z1 =C1(3,i)  ;   
        x2 =C2(1,j)  ;   y2 =C2(2,j)  ;   z2 =C2(3,j)  ;   
        x = x2 -x1 + hh 
        y = y2 -y1  
        z = z2 -z1  
        R2   = x*x +y*y +z*z 
        R    =  sqrt(R2)
        Ri   = 1.0d0/R
        cosl = x*Ri
        cosm = y*Ri 
        cosn = z*Ri 
        j0= R/dR 
        xx=  R -j0*dR
        j0= min(j0,NgridPoints+1) 
        !write(*,*)'SP:  i,j, j0=', i,j,j0
        VspSig0   =  sktab( 9, j0  ,itype,jtype)
        VspSig1   =  sktab( 9, j0+1,itype,jtype) 
        d2VspSig0 = sktab2( 9, j0  ,itype,jtype)
        d2VspSig1 = sktab2( 9, j0+1,itype,jtype)
        call cubics(xx,Vsp,dR,VspSig0,d2VspSig0,VspSig1,d2VspSig1)
        !write(*,*)' VspSig0 ,  VspSig1, d2VspSig0,d2VspSig1,Vsp='
        !write(*,*) VspSig0 ,  VspSig1, d2VspSig0,d2VspSig1,Vsp

        SspSig0   =  sktab(19, j0  ,itype,jtype)
        SspSig1   =  sktab(19, j0+1,itype,jtype) 
        d2SspSig0 = sktab2(19, j0  ,itype,jtype)
        d2SspSig1 = sktab2(19, j0+1,itype,jtype)
        call cubics(xx,Ssp,dR,SspSig0,d2SspSig0,SspSig1,d2SspSig1)
        !write(*,*)' SspSig0 , SVspSig1, d2SspSig0,d2SspSig1,Ssp='
        !write(*,*) SspSig0 ,  SspSig1, d2SspSig0,d2SspSig1,Ssp
        !write(*,*)d2SspSig1

        iblS =  i  
        iblP =  (j-1)*3  +1
        H(iblS,iblP  ) =  Vsp*cosl 
        H(iblS,iblP+1) =  Vsp*cosm 
        H(iblS,iblP+2) =  Vsp*cosn 
        S(iblS,iblP  ) =  Ssp*cosl 
        S(iblS,iblP+1) =  Ssp*cosm 
        S(iblS,iblP+2) =  Ssp*cosn 
        !---- this is transpose: 
        !H(iblP  ,iblS) =  Vsp*cosl 
        !H(iblP+1,iblS) =  Vsp*cosm 
        !H(iblP+2,iblS) =  Vsp*cosn 
        !S(iblP  ,iblS) =  Ssp*cosl 
        !S(iblP+1,iblS) =  Ssp*cosm 
        !S(iblP+2,iblS) =  Ssp*cosn 

      enddo
    enddo 

  end subroutine set_block_sp
 !@
 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !@
  subroutine set_block_PP(H,S,LDA,C1,iat1,Nat1,C2,iat2,Nat2) 
    implicit none 
    integer        :: LDA
    real (kind=8)  :: H(LDA,*), S(LDA,*)
    real (kind=8)  :: C1(3,*) ,  C2(3,*)
    real (kind=8)  :: x1,y1,z1, x2,y2,z2
    integer        :: Iat1(*), Iat2(*), Nat1, Nat2
    integer             :: i,j,k,j0,itype,jtype
    real (kind=8)       :: x,y,z,hh,R, Ri, R2,dR,  cosl,cosm,cosn,xx
    real (kind=8)       :: cosll,coslm, cosln, cosmm, cosmn, cosnn
    real (kind=8)       :: VppSig0, VppSig1, d2VppSig0, d2VppSig1, VppSig
    real (kind=8)       :: SppSig0, SppSig1, d2SppSig0, d2SppSig1, SppSig
    real (kind=8)       :: VppPi0 , VppPi1 , d2VppPi0 , d2VppPi1 , VppPi
    real (kind=8)       :: SppPi0 , SppPi1 , d2SppPi0 , d2SppPi1 , SppPi
    integer             :: iatP, iat_P0, iat_P1,  iblP
    integer             :: jatP, jat_P0, jat_P1,  jblP
    integer  ::  nGridPoints


    hh =1.0d-20
    ! set the indices for  P, P atoms

    do i=1,Nat1
       do j=1,Nat2
        itype = iat1(i)
        jtype = iat2(j)
        dR = dRGrids(itype,jtype)
        NGridPoints = NGrids(itype,jtype)  
        x1 =C1(1,i)  ;   y1 =C1(2,i)  ;   z1 =C1(3,i)  ;
        x2 =C2(1,j)  ;   y2 =C2(2,j)  ;   z2 =C2(3,j)  ;
        x = x2 -x1 + hh
        y = y2 -y1
        z = z2 -z1
        R  = sqrt(x*x + y*y +z*z)
        Ri = 1.0d0/ R
        cosl = x*Ri
        cosm = y*Ri
        cosn = z*Ri

        cosll =cosl*cosl
        coslm =cosl*cosm   ;  cosmm = cosm*cosm
        cosln =cosl*cosn   ;  cosmn = cosm*cosn  ; cosnn =cosn*cosn

        j0 = R/dR
        xx = R-j0*dR
        j0= min(j0,NgridPoints+1) 

        VppSig0   =  sktab( 6, j0  ,itype,jtype)
        VppSig1   =  sktab( 6, j0+1,itype,jtype)
        d2VppSig0 = sktab2( 6, j0  ,itype,jtype)
        d2VppSig1 = sktab2( 6, j0+1,itype,jtype)
        call cubics(xx,VppSig,dR,VppSig0,d2VppSig0,VppSig1,d2VppSig1)

        VppPi0    =  sktab( 7, j0  ,itype,jtype) 
        VppPi1    =  sktab( 7, j0+1,itype,jtype) 
        d2VppPi0  = sktab2( 7, j0  ,itype,jtype) 
        d2VppPi1  = sktab2( 7, j0+1,itype,jtype) 
        call cubics(xx,VppPi,dR,VppPi0,d2VppPi0,VppPi1,d2VppPi1)

        SppSig0   =  sktab(16, j0  ,itype,jtype) 
        SppSig1   =  sktab(16, j0+1,itype,jtype) 
        d2SppSig0 = sktab2(16, j0  ,itype,jtype) 
        d2SppSig1 = sktab2(16, j0+1,itype,jtype) 
        call cubics(xx,SppSig,dR,SppSig0,d2SppSig0,SppSig1,d2SppSig1)

        SppPi0    =  sktab(17, j0  ,itype,jtype) 
        SppPi1    =  sktab(17, j0+1,itype,jtype) 
        d2SppPi0  = sktab2(17, j0  ,itype,jtype) 
        d2SppPi1  = sktab2(17, j0+1,itype,jtype) 
        call cubics(xx,SppPi,dR,SppPi0,d2SppPi0,SppPi1,d2SppPi1)

        iblP  =  (i -1)*3 + 1
        jblP  =  (j -1)*3 + 1
        
        H(iblP  ,jblP)   =  cosll*VppSig + (1.0d0-cosll)*VppPi
        H(iblP+1,jblP)   =  coslm*VppSig -         coslm*VppPi 
        H(iblP+2,jblP)   =  cosln*VppSig -         cosln*VppPi 

        H(iblP  ,jblP+1) =  coslm*VppSig -         coslm*VppPi 
        H(iblP+1,jblP+1) =  cosmm*VppSig + (1.0d0-cosmm)*VppPi
        H(iblP+2,jblP+1) =  cosmn*VppSig -         cosmn*VppPi 

        H(iblP  ,jblP+2) =  cosln*VppSig -         cosln*VppPi 
        H(iblP+1,jblP+2) =  cosmn*VppSig -         cosmn*VppPi
        H(iblP+2,jblP+2) =  cosnn*VppSig + (1.0d0-cosnn)*VppPi 


        S(iblP  ,jblP)   =  cosll*SppSig + (1.0d0-cosll)*SppPi
        S(iblP+1,jblP)   =  coslm*SppSig -         coslm*SppPi 
        S(iblP+2,jblP)   =  cosln*SppSig -         cosln*SppPi 

        S(iblP  ,jblP+1) =  coslm*SppSig -         coslm*SppPi 
        S(iblP+1,jblP+1) =  cosmm*SppSig + (1.0d0-cosmm)*SppPi
        S(iblP+2,jblP+1) =  cosmn*SppSig -         cosmn*SppPi 

        S(iblP  ,jblP+2) =  cosln*SppSig -         cosln*SppPi 
        S(iblP+1,jblP+2) =  cosmn*SppSig -         cosmn*SppPi
        S(iblP+2,jblP+2) =  cosnn*SppSig + (1.0d0-cosnn)*SppPi 

      enddo
    enddo 


     
  end subroutine set_block_PP 
 !@
 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !@ d-orbitals:  xy, yz, zx, x^2-y^2, 3*z^2-r^2
 !@
  subroutine set_block_SD(H,S,LDA,C1,iat1,Nat1,C2,iat2,Nat2) 
    implicit none
    integer        :: LDA 
    real (kind=8)  :: H(LDA,*), S(LDA,*)  
    real (kind=8)  :: C1(3,*) ,  C2(3,*)
    real (kind=8)  :: x1,y1,z1, x2,y2,z2 
    integer        :: Iat1(*), Iat2(*), Nat1, Nat2
    integer             :: i,j,k,j0,itype,jtype
    real (kind=8)       :: x,y,z,hh,R, Ri, R2,dR,  cosl,cosm,cosn,cosll, cosmm, cosnn, xx
    real (kind=8)       :: VsdSig0, VsdSig1, d2VsdSig0, d2VsdSig1, Vsd
    real (kind=8)       :: SsdSig0, SsdSig1, d2SsdSig0, d2SsdSig1, Ssd
    
    integer  :: iblS, iblD 
    real (kind=8) :: sqrt3
    integer  ::  nGridPoints

    hh=1.0d-20   ! infinitesimal number for onsite  data (added to avoid dividing by R=zero, needed for cosines) 
     
    sqrt3=  sqrt(3.0d0)
    
    do i=1,Nat1
      do j=1,Nat2 
        itype =iat1(i) 
        jtype =iat2(j) 
        dR = dRGrids(itype,jtype)
        NGridPoints = NGrids(itype,jtype)  
        x1 =C1(1,i)  ;   y1 =C1(2,i)  ;   z1 =C1(3,i)  ;   
        x2 =C2(1,j)  ;   y2 =C2(2,j)  ;   z2 =C2(3,j)  ;   
        x = x2 -x1 + hh 
        y = y2 -y1  
        z = z2 -z1  
        R2   = x*x +y*y +z*z 
        R    =  sqrt(R2)
        Ri   = 1.0d0/R
        cosl = x*Ri
        cosm = y*Ri 
        cosn = z*Ri 
        cosll =cosl*cosl 
        cosmm =cosm*cosm 
        cosnn =cosn*cosn 

        j0= R/dR 
        j0= min(j0,NgridPoints+1) 
        xx=  R -j0*dR
        VsdSig0   =  sktab( 8, j0  ,itype,jtype)
        VsdSig1   =  sktab( 8, j0+1,itype,jtype) 
        d2VsdSig0 = sktab2( 8, j0  ,itype,jtype)
        d2VsdSig1 = sktab2( 8, j0+1,itype,jtype)
        call cubics(xx,Vsd,dR,VsdSig0,d2VsdSig0,VsdSig1,d2VsdSig1)

        SsdSig0   =  sktab(18, j0  ,itype,jtype)
        SsdSig1   =  sktab(18, j0+1,itype,jtype) 
        d2SsdSig0 = sktab2(18, j0  ,itype,jtype)
        d2SsdSig1 = sktab2(18, j0+1,itype,jtype)
        call cubics(xx,Ssd,dR,SsdSig0,d2SsdSig0,SsdSig1,d2SsdSig1)

        iblS =  i  
        iblD =  (j-1)*5  +1
        H(iblS,iblD  ) =  Vsd*sqrt3*cosl*cosm                      ! <s|xy> 
        H(iblS,iblD+1) =  Vsd*sqrt3*cosm*cosn                      ! <s|yz>
        H(iblS,iblD+2) =  Vsd*sqrt3*cosl*cosn                      ! <s|xz> 
        H(iblS,iblD+3) =  Vsd*sqrt3*(cosll-cosmm)*0.5d0            ! <s|x^2-y^2> 
        H(iblS,iblD+4) =  Vsd*(cosnn -0.5d0*(cosll+cosmm))         ! <s|3z^2-r^2>

        S(iblS,iblD  ) =  Ssd*sqrt3*cosl*cosm                      ! <s|xy> 
        S(iblS,iblD+1) =  Ssd*sqrt3*cosm*cosn                      ! <s|yz>
        S(iblS,iblD+2) =  Ssd*sqrt3*cosl*cosn                      ! <s|xz> 
        S(iblS,iblD+3) =  Ssd*sqrt3*(cosll-cosmm)*0.5d0            ! <s|x^2-y^2> 
        S(iblS,iblD+4) =  Ssd*(cosnn -0.5d0*(cosll+cosmm))         ! <s|3z^2-r^2>

      enddo
    enddo 

  end subroutine set_block_sd 


 !@
 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !@ p-orbitals:  x, y ,z
 !@ d-orbitals:  xy, yz, zx, x^2-y^2, 3*z^2-r^2
 !@
  subroutine set_block_PD(H,S,LDA,C1,iat1,Nat1,C2,iat2,Nat2) 
    implicit none
    integer       :: LDA
    real (kind=8) :: H(LDA,*),  S(LDA,*)
    real (kind=8) :: C1(3,*),   C2(3,*)
    real (kind=8) ::  x1,y1,z1,  x2, y2, z2
    integer       :: Iat1(*),  Iat2(*), Nat1, Nat2 
    integer       :: i,j,k, itype, jtype  ,j0
    real (kind=8) :: x,y,z, hh,dR, R,Ri, R2, cosl,cosm, cosn, cosll, cosmm, cosnn ,  coslmn, xx
    real (kind=8)       :: VpdSig0, VpdSig1, d2VpdSig0, d2VpdSig1, VpdSig
    real (kind=8)       :: SpdSig0, SpdSig1, d2SpdSig0, d2SpdSig1, SpdSig
    real (kind=8)       :: VpdPi0,  VpdPi1,  d2VpdPi0,  d2VpdPi1,  VpdPi
    real (kind=8)       :: SpdPi0,  SpdPi1,  d2SpdPi0,  d2SpdPi1,  SpdPi
    integer       :: iblP, iblD 
    integer       ::  iix,iiy,iiz,  iixy,iiyz,iixz, iix2y2, iiz2
    real (kind=8) ::  sqrt3, sqrt3h ,  w0, w1, aa, bb
    integer  ::  nGridPoints
    

    hh=1.0d-20   ! infinitesimal number for onsite  data (added to avoid dividing by R=zero, needed for cosines) 
     
    sqrt3  =  sqrt(3.0d0)
    sqrt3h =  0.5d0*sqrt3

    iix  = 0 ; iiy  =1 ; iiz  = 2 
    iixy = 0 ; iiyz =1 ; iixz = 2 ;  iix2y2= 3 ; iiz2 = 4
    
    do i=1,Nat1
      do j=1,Nat2 
        itype =iat1(i) 
        jtype =iat2(j) 
        dR = dRGrids(itype,jtype)
        NGridPoints = NGrids(itype,jtype)  
        x1 =C1(1,i)  ;   y1 =C1(2,i)  ;   z1 =C1(3,i)  ;   
        x2 =C2(1,j)  ;   y2 =C2(2,j)  ;   z2 =C2(3,j)  ;   
        x = x2 -x1 + hh 
        y = y2 -y1  
        z = z2 -z1  
        R2   = x*x +y*y +z*z 
        R    =  sqrt(R2)
        Ri   = 1.0d0/R
        cosl = x*Ri
        cosm = y*Ri 
        cosn = z*Ri 
        cosll  = cosl*cosl 
        cosmm  = cosm*cosm 
        cosnn  = cosn*cosn 
        coslmn = cosl*cosm*cosn

        j0= R/dR 
        j0= min(j0,NgridPoints+1) 
        xx=  R -j0*dR

        VpdSig0   =  sktab( 4, j0  ,itype,jtype)
        VpdSig1   =  sktab( 4, j0+1,itype,jtype) 
        d2VpdSig0 = sktab2( 4, j0  ,itype,jtype)
        d2VpdSig1 = sktab2( 4, j0+1,itype,jtype)
        call cubics(xx,VpdSig,dR,VpdSig0,d2VpdSig0,VpdSig1,d2VpdSig1)
 

        VpdPi0    =  sktab( 5, j0  ,itype,jtype)
        VpdPi1    =  sktab( 5, j0+1,itype,jtype) 
        d2VpdPi0  = sktab2( 5, j0  ,itype,jtype)
        d2VpdPi1  = sktab2( 5, j0+1,itype,jtype)
        call cubics(xx,VpdPi,dR,VpdPi0,d2VpdPi0,VpdPi1,d2VpdPi1)

        SpdSig0   =  sktab(14, j0  ,itype,jtype)
        SpdSig1   =  sktab(14, j0+1,itype,jtype) 
        d2SpdSig0 = sktab2(14, j0  ,itype,jtype)
        d2SpdSig1 = sktab2(14, j0+1,itype,jtype)
        call cubics(xx,SpdSig,dR,SpdSig0,d2SpdSig0,SpdSig1,d2SpdSig1)

        SpdPi0    =  sktab(15, j0  ,itype,jtype)
        SpdPi1    =  sktab(15, j0+1,itype,jtype) 
        d2SpdPi0  = sktab2(15, j0  ,itype,jtype)
        d2SpdPi1  = sktab2(15, j0+1,itype,jtype)
        call cubics(xx,SpdPi,dR,SpdPi0,d2SpdPi0,SpdPi1,d2SpdPi1)

        iblP = (i-1)*3 +1 
        iblD = (j-1)*5 +1 

! basis set order:
!    p: x,y,z
!    d: xy,yz, zx,  x2-y2,  3z2-r2

   ! <x|xy> :
        w0 = sqrt3*cosll*cosm  ;  w1 =  cosm*(1.d0-2.0d0*cosll)         !  <x|xy> 
             H(iblP+iix, iblD+iixy) = w0*VpdSig  + w1*VpdPi             
             S(iblP+iix, iblD+iixy) = w0*SpdSig  + w1*SpdPi             
   ! <x|xz>
        w0 = sqrt3*cosll*cosn  ;  w1 =  cosn*(1.d0-2.0d0*cosll)          !  <x|xz> 
             H(iblP+iix, iblD+iixz) = w0*VpdSig  + w1*VpdPi             
             S(iblP+iix, iblD+iixz) = w0*SpdSig  + w1*SpdPi             
   ! <y|yz>
        w0 = sqrt3*cosmm*cosn  ;  w1 =  cosn*(1.d0-2.0d0*cosmm)          !  <y|yz> 
             H(iblP+iiy, iblD+iiyz) = w0*VpdSig  + w1*VpdPi             
             S(iblP+iiy, iblD+iiyz) = w0*SpdSig  + w1*SpdPi             
   ! <y|xy>
        w0 = sqrt3*cosmm*cosl  ;  w1 =  cosl*(1.d0-2.0d0*cosmm)          !  <y|xy> 
             H(iblP+iiy, iblD+iixy) = w0*VpdSig  + w1*VpdPi             
             S(iblP+iiy, iblD+iixy) = w0*SpdSig  + w1*SpdPi             
   ! <z|xz>
        w0 = sqrt3*cosnn*cosl  ;  w1 =  cosl*(1.d0-2.0d0*cosnn)          !  <z|xz> 
             H(iblP+iiz, iblD+iixz) = w0*VpdSig  + w1*VpdPi             
             S(iblP+iiz, iblD+iixz) = w0*SpdSig  + w1*SpdPi             
   ! <z|xz>
        w0 = sqrt3*cosnn*cosm  ;  w1 =  cosm*(1.d0-2.0d0*cosnn)          !  <z|yz> 
             H(iblP+iiz, iblD+iiyz) = w0*VpdSig  + w1*VpdPi             
             S(iblP+iiz, iblD+iiyz) = w0*SpdSig  + w1*SpdPi             


!--- these 3 are identical:
        w0 =  sqrt3*coslmn      ;   w1  =             -2.0d0*coslmn          !  <x|yz>= <y|zx>=<z|xy>
            
        aa =   w0*VpdSig  + w1*VpdPi
        bb =   w0*SpdSig  + w1*SpdPi

   !  <x|yz>:
             H(iblP+iix, iblD+iiyz) = aa 
             S(iblP+iix, iblD+iiyz) = bb
   !  <y|xz>:
             H(iblP+iiy, iblD+iixz) = aa 
             S(iblP+iiy, iblD+iixz) = bb
   !  <z|xy>:
             H(iblP+iiz, iblD+iixy) = aa 
             S(iblP+iiz, iblD+iixy) = bb

!----


        aa = cosll - cosmm ;      bb  =  sqrt3h*aa ; 
   !  <x|x^2-y^2>:
        w0 =  cosl*bb                ;   w1  =   cosl*(1.0d0 - aa)                  !  <x|x^2-y^2>
             H(iblP+iix, iblD+iix2y2) = w0*VpdSig  + w1*VpdPi             
             S(iblP+iix, iblD+iix2y2) = w0*SpdSig  + w1*SpdPi             
   !  <y|x^2-y^2>:
        w0 =  cosm*bb                ;   w1  =  -cosm*(1.0d0 + aa)                  !  <y|x^2-y^2>
             H(iblP+iiy, iblD+iix2y2) = w0*VpdSig  + w1*VpdPi             
             S(iblP+iiy, iblD+iix2y2) = w0*SpdSig  + w1*SpdPi             
   !  <z|x^2-y^2>:
        w0 =  cosn*bb                ;   w1  =  -cosn*aa                            !  <z|x^2-y^2>
             H(iblP+iiz, iblD+iix2y2) = w0*VpdSig  + w1*VpdPi             
             S(iblP+iiz, iblD+iix2y2) = w0*SpdSig  + w1*SpdPi             
 
        aa =  cosll + cosmm ;      bb  = cosnn - 0.5d0 * aa  
   !  <x|2z^2-r^2>:
        w0 =  cosl*bb                ; w1 = -sqrt3*cosl*cosnn                       !  <x|2z^2-r^2>      
             H(iblP+iix, iblD+iiz2) = w0*VpdSig  + w1*VpdPi             
             S(iblP+iix, iblD+iiz2) = w0*SpdSig  + w1*SpdPi             
   !  <y|2z^2-r^2>:
        w0 =  cosm*bb                ; w1 = -sqrt3*cosm*cosnn                       !  <y|2z^2-r^2>
             H(iblP+iiy, iblD+iiz2) = w0*VpdSig  + w1*VpdPi             
             S(iblP+iiy, iblD+iiz2) = w0*SpdSig  + w1*SpdPi             
   !  <z|2z^2-r^2>:
        w0 =  cosn*bb                ; w1 =  sqrt3*cosn*aa                          !  <z|2z^2-r^2>
             H(iblP+iiz, iblD+iiz2) = w0*VpdSig  + w1*VpdPi             
             S(iblP+iiz, iblD+iiz2) = w0*SpdSig  + w1*SpdPi             

      enddo
    enddo

        
    


  end subroutine set_block_pd

 !@
 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !@ d-orbitals:  xy, yz, zx, x^2-y^2, 3*z^2-r^2
 !@
  subroutine set_block_DD(H,S,LDA,C1,iat1,Nat1,C2,iat2,Nat2) 
    !use  utilmod
    implicit none
    integer       :: LDA
    real (kind=8) :: H(LDA,*), S(LDA,*)
    real (kind=8) :: C1(3,*), C2(3,*)
    real (kind=8) :: x1,y1,z1, x2,y2,z2
    integer       :: Iat1(*),  Iat2(*), NAt1, Nat2
    integer       :: i,j,k, itype,jtype, j0
    real (kind=8) :: x,y,z, hh,dR,R,Ri,R2, cosl,cosm,cosn, cosll,cosmm,cosnn, coslm,cosln,cosmn,xx

    real (kind=8) :: sqrt3 , one,  half,  quart , three4
    real (kind=8)       :: VddSig0, VddSig1, d2VddSig0, d2VddSig1, VddSig
    real (kind=8)       :: SddSig0, SddSig1, d2SddSig0, d2SddSig1, SddSig
    real (kind=8)       :: VddPi0,  VddPi1,  d2VddPi0,  d2VddPi1,  VddPi
    real (kind=8)       :: SddPi0,  SddPi1,  d2SddPi0,  d2SddPi1,  SddPi
    real (kind=8)       :: VddDel0, VddDel1, d2VddDel0, d2VddDel1, VddDel
    real (kind=8)       :: SddDel0, SddDel1, d2SddDel0, d2SddDel1, SddDel

    real (kind=8)       ::  w0,w1,w2 , aa0,aa1,aa2, aa

    integer       :: iblD1, IblD2
    integer       :: iixy,iiyz,iixz, iix2y2, iiz2
    integer  ::  nGridPoints


    iixy = 0 ; iiyz =1 ; iixz = 2 ;  iix2y2= 3 ; iiz2 = 4     !basis set mapping
    sqrt3 = sqrt(3.0d0) ;  one=1.0d0;  half = 0.5d0 ;  quart =0.25d0 ;  three4= 0.75d0 

    hh=1.0d-20   ! infinitesimal number for onsite  data (added to avoid dividing by R=zero, needed for cosines) 

    do i=1,Nat1
      do j=1,Nat2
        itype = iat1(i) 
        jtype = iat2(j) 
        dR = dRGrids(itype,jtype)
        NGridPoints = NGrids(itype,jtype)  
        x1 =C1(1,i)  ;   y1 =C1(2,i)  ;   z1 =C1(3,i)  ;   
        x2 =C2(1,j)  ;   y2 =C2(2,j)  ;   z2 =C2(3,j)  ;   
        x = x2 -x1 + hh 
        y = y2 -y1  
        z = z2 -z1  
        R2   = x*x +y*y +z*z 
        R    =  sqrt(R2)
        Ri   = 1.0d0/R
        cosl = x*Ri
        cosm = y*Ri 
        cosn = z*Ri 
        cosll  = cosl*cosl 
        coslm  = cosl*cosm 
        cosln  = cosl*cosn 
        cosmm  = cosm*cosm 
        cosmn  = cosm*cosn 
        cosnn  = cosn*cosn 

        j0= R/dR 
        j0= min(j0,NgridPoints+1) 
        xx=  R -j0*dR

        VddSig0   =  sktab( 1, j0  ,itype,jtype)
        VddSig1   =  sktab( 1, j0+1,itype,jtype) 
        d2VddSig0 = sktab2( 1, j0  ,itype,jtype)
        d2VddSig1 = sktab2( 1, j0+1,itype,jtype)
        call cubics(xx,VddSig,dR,VddSig0,d2VddSig0,VddSig1,d2VddSig1)

        VddPi0    =  sktab( 2, j0  ,itype,jtype)
        VddPi1    =  sktab( 2, j0+1,itype,jtype) 
        d2VddPi0  = sktab2( 2, j0  ,itype,jtype)
        d2VddPi1  = sktab2( 2, j0+1,itype,jtype)
        call cubics(xx,VddPi,dR,VddPi0,d2VddPi0,VddPi1,d2VddPi1)

        VddDel0   =  sktab( 3, j0  ,itype,jtype)
        VddDel1   =  sktab( 3, j0+1,itype,jtype) 
        d2VddDel0 = sktab2( 3, j0  ,itype,jtype)
        d2VddDel1 = sktab2( 3, j0+1,itype,jtype)
        call cubics(xx,VddDel,dR,VddDel0,d2VddDel0,VddDel1,d2VddDel1)

        SddSig0   =  sktab(11, j0  ,itype,jtype)
        SddSig1   =  sktab(11, j0+1,itype,jtype) 
        d2SddSig0 = sktab2(11, j0  ,itype,jtype)
        d2SddSig1 = sktab2(11, j0+1,itype,jtype)
        call cubics(xx,SddSig,dR,SddSig0,d2SddSig0,SddSig1,d2SddSig1)

        SddPi0    =  sktab(12, j0  ,itype,jtype)
        SddPi1    =  sktab(12, j0+1,itype,jtype) 
        d2SddPi0  = sktab2(12, j0  ,itype,jtype)
        d2SddPi1  = sktab2(12, j0+1,itype,jtype)
        call cubics(xx,SddPi,dR,SddPi0,d2SddPi0,SddPi1,d2SddPi1)

        SddDel0   =  sktab(13, j0  ,itype,jtype)
        SddDel1   =  sktab(13, j0+1,itype,jtype) 
        d2SddDel0 = sktab2(13, j0  ,itype,jtype)
        d2SddDel1 = sktab2(13, j0+1,itype,jtype)
        call cubics(xx,SddDel,dR,SddDel0,d2SddDel0,SddDel1,d2SddDel1)

        iblD1 = (i-1)*5 +1 
        iblD2 = (j-1)*5 +1 

! basis set order:
!    d: xy,yz, zx,  x2-y2,  3z2-r2

   !--- diagonal:
   ! <xy|xy> :
         w0 = 3.0d0*cosll*cosmm  ;    w1 = cosll +cosmm -4.0d0*cosll*cosmm     ; w2 =  cosnn +cosll*cosmm 
              H(iblD1 +iixy  ,  iblD2 +iixy) =  w0*VddSig   + w1*VddPi   +  w2*VddDel  
              S(iblD1 +iixy  ,  iblD2 +iixy) =  w0*SddSig   + w1*SddPi   +  w2*SddDel  
   ! <yz|yz> :
         w0 = 3.0d0*cosmm*cosnn  ;    w1 = cosmm +cosnn -4.0d0*cosmm*cosnn     ; w2 =  cosll +cosmm*cosnn 
              H(iblD1 +iiyz  ,  iblD2 +iiyz) =  w0*VddSig   + w1*VddPi   +  w2*VddDel  
              S(iblD1 +iiyz  ,  iblD2 +iiyz) =  w0*SddSig   + w1*SddPi   +  w2*SddDel  
   ! <xz|xz> :
         w0 = 3.0d0*cosll*cosnn  ;    w1 = cosnn +cosll -4.0d0*cosll*cosnn     ; w2 =  cosmm +cosll*cosnn 
              H(iblD1 +iixz  ,  iblD2 +iixz) =  w0*VddSig   + w1*VddPi   +  w2*VddDel  
              S(iblD1 +iixz  ,  iblD2 +iixz) =  w0*SddSig   + w1*SddPi   +  w2*SddDel  

   !--- nondiagonal:
   ! <xy|yz> = <yz|xy> : 
         w0 = 3.0d0*cosln*cosmm  ;    w1 = cosln*(one-4.0d0*cosmm)     ; w2 =  cosln*(cosmm -one)
              aa0  =  w0*VddSig   + w1*VddPi   +  w2*VddDel
              H(iblD1 +iixy  ,  iblD2 +iiyz) =  aa0 
              H(iblD1 +iiyz  ,  iblD2 +iixy) =  aa0 
              aa0  =  w0*SddSig   + w1*SddPi   +  w2*SddDel
              S(iblD1 +iixy  ,  iblD2 +iiyz) =  aa0 
              S(iblD1 +iiyz  ,  iblD2 +iixy) =  aa0 

   ! <xy|xz> = <xz|xy> :
         w0 = 3.0d0*cosll*cosmn  ;    w1 = cosmn*(one-4.0d0*cosll)     ; w2 =  cosmn*(cosll -one)
              aa0  =  w0*VddSig   + w1*VddPi   +  w2*VddDel
              H(iblD1 +iixy  ,  iblD2 +iixz) =  aa0 
              H(iblD1 +iixz  ,  iblD2 +iixy) =  aa0 
              aa0  =  w0*SddSig   + w1*SddPi   +  w2*SddDel
              S(iblD1 +iixy  ,  iblD2 +iixz) =  aa0 
              S(iblD1 +iixz  ,  iblD2 +iixy) =  aa0 
            
   ! <yz|xz> = <xz|yz> :
         w0 = 3.0d0*cosnn*coslm  ;    w1 = coslm*(one-4.0d0*cosnn)     ; w2 =  coslm*(cosnn -one) 
              aa0  =  w0*VddSig   + w1*VddPi   +  w2*VddDel
              H(iblD1 +iiyz  ,  iblD2 +iixz) =  aa0 
              H(iblD1 +iixz  ,  iblD2 +iiyz) =  aa0 
              aa0  =  w0*SddSig   + w1*SddPi   +  w2*SddDel
              S(iblD1 +iiyz  ,  iblD2 +iixz) =  aa0 
              S(iblD1 +iixz  ,  iblD2 +iiyz) =  aa0 

   !-- next batch:
         aa =  (cosll-cosmm)  ; aa0= 1.5d0*aa ; aa1 =2.0d0*aa ; aa2 = 0.5d0*aa  

   ! <xy|x^2-y^2> = <x^2-y2|xy> : 
         w0 = aa0 * coslm  ;  w1 = -coslm*aa1          ;  w2 =  coslm*aa2  
              aa  =  w0*VddSig   +  w1*VddPi  +  w2*VddDel 
              H(iblD1 +iixy   , iblD2 +iix2y2) = aa 
              H(iblD1 +iix2y2 , iblD2 +iixy  ) = aa 
              aa  =  w0*SddSig   +  w1*SddPi  +  w2*SddDel 
              S(iblD1 +iixy   , iblD2 +iix2y2) = aa 
              S(iblD1 +iix2y2 , iblD2 +iixy  ) = aa 

   ! <yz|x^2-y^2> = <x^2-y^2|xz> :
         w0 = aa0 * cosmn  ;  w1 = -cosmn*(one +aa1)   ;  w2 =  cosmn*(one + aa2)
              aa  =  w0*VddSig   +  w1*VddPi  +  w2*VddDel 
              H(iblD1 +iiyz   , iblD2 +iix2y2) = aa 
              H(iblD1 +iix2y2 , iblD2 +iiyz  ) = aa 
              aa  =  w0*SddSig   +  w1*SddPi  +  w2*SddDel 
              S(iblD1 +iiyz   , iblD2 +iix2y2) = aa 
              S(iblD1 +iix2y2 , iblD2 +iiyz  ) = aa 

   ! <xz|x^2-y^2> = <x^2-y^2|xz> : 
         w0 = aa0 * cosln  ;  w1 =  cosln*(one -aa1)   ;  w2 = -cosln*(one - aa2)
              aa  =  w0*VddSig   +  w1*VddPi  +  w2*VddDel 
              H(iblD1 +iixz   , iblD2 +iix2y2) = aa 
              H(iblD1 +iix2y2 , iblD2 +iixz  ) = aa 
              aa  =  w0*SddSig   +  w1*SddPi  +  w2*SddDel 
              S(iblD1 +iixz   , iblD2 +iix2y2) = aa 
              S(iblD1 +iix2y2 , iblD2 +iixz  ) = aa 

   !-- next batch:
         aa = cosll +cosmm ; aa2= 0.5d0*aa ;  aa1 = aa -cosnn ;  aa0 = (cosnn -aa2)  

   ! <xy|3z^2-r2> = <3z^2-r2|xy> : 
         w0 = coslm*aa0    ;  w1 = -2.0d0*coslm*cosnn  ;  w2 =  0.5d0 * coslm*(one +cosnn)
              aa  =  (w0*VddSig   +  w1*VddPi  +  w2*VddDel )*sqrt3
              H(iblD1 +iixy   , iblD2 +iiz2) = aa 
              H(iblD1 +iiz2   , iblD2 +iixy) = aa 
              aa  =  (w0*SddSig   +  w1*SddPi  +  w2*SddDel )*sqrt3
              S(iblD1 +iixy   , iblD2 +iiz2) = aa 
              S(iblD1 +iiz2   , iblD2 +iixy) = aa 

   ! <yz|3z^2-r2> = <3z^2-r2|yz> : 
         w0 = cosmn*aa0    ;  w1 =  cosmn*aa1          ;  w2 = -cosmn*aa2
              aa  =  (w0*VddSig   +  w1*VddPi  +  w2*VddDel )*sqrt3
              H(iblD1 +iiyz   , iblD2 +iiz2) = aa 
              H(iblD1 +iiz2   , iblD2 +iiyz) = aa 
              aa  =  (w0*SddSig   +  w1*SddPi  +  w2*SddDel )*sqrt3
              S(iblD1 +iiyz   , iblD2 +iiz2) = aa 
              S(iblD1 +iiz2   , iblD2 +iiyz) = aa 

   ! <xz|3z^2-r2> = <3z^2-r2|xz> : 
         w0 = cosln*aa0    ;  w1 =  cosln*aa1          ;  w2 = -cosln*aa2
              aa  =  (w0*VddSig   +  w1*VddPi  +  w2*VddDel )*sqrt3
              H(iblD1 +iixz   , iblD2 +iiz2) = aa 
              H(iblD1 +iiz2   , iblD2 +iixz) = aa 
              aa  =  (w0*SddSig   +  w1*SddPi  +  w2*SddDel )*sqrt3
              S(iblD1 +iixz   , iblD2 +iiz2) = aa 
              S(iblD1 +iiz2   , iblD2 +iixz) = aa 

   !-- final batch:
        aa0= cosll+cosmm ; aa1 =  cosll-cosmm

   ! <x^2-y^2| x^2-y^2>:  
        w0  = 0.75d0 * aa1*aa1  ;   w1 = aa0 - aa1*aa1   ;  w2 = cosnn + 0.25d0*aa1*aa1   
              H(iblD1 +iix2y2 , iblD2 +iix2y2) =  w0*VddSig   +  w1*VddPi  +  w2*VddDel  
              S(iblD1 +iix2y2 , iblD2 +iix2y2) =  w0*SddSig   +  w1*SddPi  +  w2*SddDel 
        
   ! <x^2-y^2|3z^2-r2> =  <3z^2-r2| x2-y2> : 
        w0  = 0.5d0*aa1*(cosnn -0.5d0*aa0)   ;  w1 = -cosnn*aa1     ;     w2 = 0.25d0*(one+cosnn)*aa1
              aa =  sqrt3*(w0*VddSig +  w1*VddPi  + w2*VddDel) 
              H(iblD1 +iix2y2 , iblD2 +iiz2  ) = aa 
              H(iblD1 +iiz2   , iblD2 +iix2y2) = aa 
              aa =  sqrt3*(w0*SddSig +  w1*SddPi  + w2*SddDel) 
              S(iblD1 +iix2y2 , iblD2 +iiz2  ) = aa 
              S(iblD1 +iiz2   , iblD2 +iix2y2) = aa 

   ! <3z^2-r2|3z^2-r2> : 
        aa2 = cosnn - 0.5d0*aa0 
        w0  = aa2*aa2    ;  w1  = 3.0d0*cosnn*aa0  ;  w2 = 0.75d0*aa0*aa0 
              H(iblD1 +iiz2 , iblD2 +iiz2) =  w0*VddSig   +  w1*VddPi  +  w2*VddDel  
              S(iblD1 +iiz2 , iblD2 +iiz2) =  w0*SddSig   +  w1*SddPi  +  w2*SddDel 

        !write(*,*)'::**** Htemp=' ; call  dumpf(H,10)
        !write(*,*)'::**** Stemp=' ; call  dumpfa(S,10,10)     

      enddo
    enddo 


  end subroutine set_block_dd
!@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@
  subroutine sortatoms(Natoms,ntype,C0,iat0,C1,iat)
  !subroutine sortatoms (Natoms,ntype,C0,C1)
  ! mapat0   =:>   |iat |  lmax| ioff  |new
    implicit none 
    integer :: Natoms,ntype, i,j,k,itype
    real (kind=8) :: C0(3,Natoms), C1(3,Natoms)
    integer :: NatS,NatP,NatD,NatF,NatO
    integer :: IatS,IatP,IatD,IatF,IatO
    integer ::ioff, lmaxsh, lmaxmin,lmaxmax
    integer :: iat0(:),iat(:)
    integer :: i0,i1
    
    allocate(mapat0(Natoms,4),mapat1(Natoms,4))
         
    lmaxmin= lmax(1)
    lmaxmax= lmax(1)

     NatS =0 
     NatP =0 
     NatD =0 
     NatF =0 
     NatO = 0  ! other

    do i=1,ntype 
      lmaxmin =min(lmaxmin,lmax(i))
      lmaxmax =min(lmaxmax,lmax(i))
      enddo

    do i=1,Natoms 
      itype =lmax(iat0(i))
      if (itype==1)   then 
        NatS =NatS+1
      elseif (itype==2)   then 
        NatP =NatP+1
      elseif (itype==3)   then 
        NatD =NatD+1
      elseif (itype==4)   then 
        NatF =NatF+1
      else  
        NatO = NatO +1
      endif 
    enddo
    nbasis = NatS +4*NatP + 9*NatD +16*NatF
    allocate(mapbas(nbasis,2))
        
!  set iat offsets  for input 
   IatS = 1 
   IatP = IatS +NatS 
   IatD = IatP +NatP
   IatF = IatD +NatD

   ioff =1 
  
   do i=1,Natoms 
     mapat0(i,1)  = iat0(i)                             ! atom type for SK file
     lmaxsh       =  lmax(iat0(i))
     mapat0(i,2)  =  lmaxsh                            ! max  Shell type   /S=1,P=2, etc./
     mapat0(i,3)  = ioff                               ! basis  set offset in matrices   starts here 
     ioff = ioff + lmaxsh *lmaxsh 
     if (lmaxsh ==1)  then 
         mapat0(i,4)  = IatS 
         IatS = IatS +1
     elseif (lmaxsh==2) then     
         mapat0(i,4)  = IatP 
         IatP = IatP +1
     elseif (lmaxsh==3) then     
         mapat0(i,4)  = IatD 
         IatD = IatD +1
     elseif (lmaxsh==4) then     
         mapat0(i,4)  = IatF 
         IatF = IatF +1
     endif 
   enddo
  ! Ioff on exitiof above loop should  equal to nbasis
   if (nbasis/= ioff-1) write(*,*)'WARNING: check basis  number in sortatoms!,ioff=',ioff
!
!  set inverse map for sorted   atoms
   do i=1,Natoms 
     j = mapat0(i,4)  
     C1(1,j)  = C0(1,i)
     C1(2,j)  = C0(2,i)
     C1(3,j)  = C0(3,i)
     mapat1(j,1) =   mapat0(i,1)
     mapat1(j,2) =   mapat0(i,2)
     mapat1(j,3) =   mapat0(i,3)
     mapat1(j,4) =    i 
   enddo
!
!  return fixed iat
  iat = mapAt1(:,1)

!
! S-shell
   Ioff = 1
   do i=1, Natoms 
     mapbas(i,1) = mapat1(i,3)
     ioff = ioff+1 
   enddo
   
!
! P-shell
   do i= NatS+1, Natoms 
     mapbas(Ioff+0,1) = mapat1(i,3) +  1
     mapbas(Ioff+1,1) = mapat1(i,3) +  2
     mapbas(Ioff+2,1) = mapat1(i,3) +  3
     Ioff =Ioff +3
   enddo

!
! D-shell
   do i=NatS+NatP+1,  Natoms 
     mapbas(Ioff+0,1) = mapat1(i,3) +  4
     mapbas(Ioff+1,1) = mapat1(i,3) +  5
     mapbas(Ioff+2,1) = mapat1(i,3) +  6
     mapbas(Ioff+3,1) = mapat1(i,3) +  7
     mapbas(Ioff+4,1) = mapat1(i,3) +  8
     Ioff =Ioff +5
   enddo 

! F-shell
   do i=NatS+NatP+NatD+1,  Natoms 
     mapbas(Ioff+0,1) = mapat1(i,3) +  9
     mapbas(Ioff+1,1) = mapat1(i,3) + 10
     mapbas(Ioff+2,1) = mapat1(i,3) + 11
     mapbas(Ioff+3,1) = mapat1(i,3) + 12
     mapbas(Ioff+4,1) = mapat1(i,3) + 13
     mapbas(Ioff+5,1) = mapat1(i,3) + 14
     mapbas(Ioff+6,1) = mapat1(i,3) + 15
     Ioff =Ioff +7
   enddo 
   
! inverse map: input -> working
   do i=1,nbasis
     j =  mapbas(i,1) 
     mapbas(j,2) = i 
   enddo 

  end subroutine sortatoms 
!@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@
  subroutine trnsl2i(A0,A1,mapb,nbasis)
  implicit none 
  real (kind=8) ::  A0(nbasis,nbasis),A1(nbasis,nbasis)
  integer       ::  i0,j0,i1,j1, nbasis,natoms
  integer       ::  mapb(*)
  
  do i0= 1,Nbasis
     i1= MapB(i0)
     do j0= 1,nbasis 
        j1= MapB(j0)
        !write(*,*) i0,j0, ' -> ',i1,j1
        A1(i1,j1) =  A0(i0,j0)
    enddo
  enddo   
!  stop


  end subroutine
!@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  subroutine fill(A,nbasis)
    real (kind=8) :: A(nbasis,nbasis)
    integer   ::  i,j

    do i=1,nbasis
      do j=1,nbasis
        A(i,j) = 0.d0 + i + 0.01*j
      enddo
    enddo 
  end subroutine fill
!@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine   dumpf(A,m)
    integer :: i,i0,i1, j,j0,j1,m
    double precision, dimension(:,:) :: A
    100 format (100(f12.5,' '))
    i0 = lbound(A,1)  ; i1 = min(ubound(A,1), i0+m-1)
    j0 = lbound(A,2)  ; j1 = min(ubound(A,2), j0+m-1)
    do i=i0,i1
      write(*,100) (A(i,j),j=j0,j1)
     enddo
   !write(*,*)  'size', size(A)
  end subroutine   dumpf


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine   dumpfa(A,M,N)
  implicit none  
  integer  ::  M, N, i, j,k
  real (kind=8) ::  A(M,N) 
  do  i=1,M
    write(*,100) (A(i,j),j=1,N) 
    enddo 
  100  format (100(F12.5,' '))
 
  end subroutine   dumpfa

end module


