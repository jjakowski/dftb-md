
FTN= gfortran -g  -fbacktrace     -fbounds-check

OBJ = utilmod2.o datamod.o splinemod.o timing.o dochol.o  slkomod.o formsh.o  formg.o setdens.o   prog.o 


all:
#	pgfortran -o eldyn  utilmod2.f90 timing.f90 initmd.f90 dochol.f90 diagev.f90  commutp.f90 prog.f90  -lblas -llapack
#	#pgfortran -o newdftb utilmod2.f90 datamod.f90 splinemod.f90 timing.f90  slkomod.f90 formsh.f90  prog.f90 
#-- mMAIN:
	$(FTN)  -o newdftb utilmod2.f90  datamod.f90 splinemod.f90 timing.f90 dochol.f90  slkomod.f90 formsh.f90  formg.f90 setdens.f90 scf.f90 prog.f90  -lblas -llapack
#	 pgfortran  -o  gam  timing.f90 datamod.f90 utilmod2.f90 splinemod.f90  slkomod.f90  formsh.f90 formg.f90  prog-gamm.f90   -lblas -llapack
#	pgfortran -o newdftb ${OBJ}  -lblas -llapack 



#------
#  gamma:  just generate gamma matrix
gamma: timing.f90 datamod.f90 utilmod2.f90 splinemod.f90  slkomod.f90  formsh.f90 formg.f90  prog-gamm.f90    
	 pgfortran  -o  gam  timing.f90 datamod.f90 utilmod2.f90 splinemod.f90  slkomod.f90  formsh.f90 formg.f90  prog-gamm.f90   -lblas -llapack
	 #pgfortran  -o  gam datamod.f90 utilmod2.f90 splinemod.f90  slkomod.f90  formsh.f90 prog-gamm.f90  

#------
#   this one solves DIIS problem in Pulay and   diagonalizatin (Jacek's)  version
#   we use random  error matrice s B
#  The theory of Pulay's  diss is described  here   http://vergil.chemistry.gatech.edu/notes/diis/node2.html
#  
diis:   prog-diis.f90
	pgfortran  -o  diis utilmod2.f90  prog-diis.f90   -lblas -llapack 

#------
#  SCF  scenario  that employs diis routines
#   iter=0  :  Hcore(C), init/guess  P--> G(P), F=H+G(P), check the error= [P,F] 
#   iter=1  :  diag(F) ,  -->  P,  G(P), F=H+G(P), check the error= [P,F]
#   iter=2  :  diag(F) ,  -->  P,  G(P), F=H+G(P), check the error= [P,F]    
#   ....
#    
scf:   prog-diis-scf.f90
	#pgfortran  -o  scf utilmod2.f90  prog-diis-scf.f90   -lblas -llapack 
	$(FTN)   -o  scf utilmod2.f90  prog-diis-scf.f90   -lblas -llapack 



clean:
	rm ./*.o  ./*.mod


