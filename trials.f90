PROGRAM TRIAL
  	USE TENSOR_UTILS

	IMPLICIT NONE

	REAL*4, DIMENSION(3,3) :: smat1,smat2
	REAL*8, DIMENSION(3,3) :: dmat1,dmat2
	REAL*4, DIMENSION(9,9) :: sreskron
	REAL*8, DIMENSION(9,9) :: dreskron
	REAL*4, DIMENSION(9,3) :: sreskhrao
	REAL*8, DIMENSION(9,3) :: dreskhrao
	REAL*4, DIMENSION(3,3) :: sreshad
	REAL*8, DIMENSION(3,3) :: dreshad	
	INTEGER*4 :: ii
	REAL*4, DIMENSION(28*28,60000) :: mnist
	REAL*4, DIMENSION(2,2,2,2,2,2,2) :: prova
	TYPE(STENSOR) :: mytens

	smat2 = reshape((/ 1E0, 2E0, 3E0, 4E0, 5E0, 6E0, 7E0, 8E0, 9E0 /), shape(smat2))
	smat1 = reshape((/ 1E0, 4E0, 7E0, 2E0, 5E0, 8E0, 3E0, 6E0, 9E0 /), shape(smat1))
	dmat2 = reshape((/ 1D0, 2D0, 3D0, 4D0, 5D0, 6D0, 7D0, 8D0, 9D0 /), shape(dmat2))
	dmat1 = reshape((/ 1D0, 4D0, 7D0, 2D0, 5D0, 8D0, 3D0, 6D0, 9D0 /), shape(dmat1))

!	DO ii=1,SIZE(smat1,1)
!		PRINT*, smat1(ii,:)
!	END DO
!	DO ii=1,SIZE(dmat1,1)
!		PRINT*, dmat1(ii,:)
!	END DO
!	DO ii=1,SIZE(smat2,1)
!		PRINT*, smat2(ii,:)
!	END DO
!	DO ii=1,SIZE(dmat2,1)
!		PRINT*, dmat2(ii,:)
!	END DO

	sreskron = (smat1.KRON.smat2)
	dreskron = (dmat1.KRON.dmat2)
	sreskhrao = (smat1.KHRAO.smat2)
	dreskhrao = (dmat1.KHRAO.dmat2)
	sreshad = smat1.HAD.smat2	
	dreshad = dmat1.HAD.dmat2

!	DO ii=1,SIZE(sreskron,1)
!		PRINT*, sreskron(ii,:)
!	END DO
!	DO ii=1,SIZE(dreskron,1)
!		PRINT*, dreskron(ii,:)
!	END DO
!	DO ii=1,SIZE(sreskhrao,1)
!		PRINT*, sreskhrao(ii,:)
!	END DO
!	DO ii=1,SIZE(dreskhrao,1)
!		PRINT*, dreskhrao(ii,:)
!	END DO
!	DO ii=1,SIZE(sreshad,1)
!		PRINT*, sreshad(ii,:)
!	END DO
!	DO ii=1,SIZE(dreshad,1)
!		PRINT*, dreshad(ii,:)
!	END DO

	
!	mnist = READ_MNIST('./mnist_train.csv',.TRUE.)
	mytens = STENSOR(order=3, modes=(/2,2,2/))
	mytens%elems = (/1,3,5,7,2,4,6,8/)
	print*, mytens%elems
!	mytens%elems = PACK(mnist,.TRUE.)
!	print*, PACK(smat1,.TRUE.)
!	print*, smat1
	prova=0e0
	print*, prova(1,1,1,1,1,1,:)



	

END PROGRAM TRIAL
