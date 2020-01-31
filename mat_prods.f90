MODULE MAT_PRODS
	USE TENSOR_TYPES
	IMPLICIT NONE
	REAL*8, PARAMETER :: pi=ACOS(-1d0)



! INTERFACES FOR COMMON BINARY OPERATORS
	INTERFACE OPERATOR (.KRON.)
		MODULE PROCEDURE S_KRON, D_KRON, SVEC_KRON, DVEC_KRON
	END INTERFACE OPERATOR (.KRON.)

	INTERFACE OPERATOR (.KHRAO.)
		MODULE PROCEDURE S_KHRAO, D_KHRAO
	END INTERFACE OPERATOR (.KHRAO.)

	INTERFACE OPERATOR (.HAD.)
		MODULE PROCEDURE S_HADAMARD, D_HADAMARD
	END INTERFACE OPERATOR (.HAD.)




CONTAINS



	FUNCTION S_KRON(mat1, mat2)
	!====================================================================
	!Computes the single precision Kronecker product of the given single precision matrices.
	!INPUT:
	!- mat1, mat2	: (REAL*4) the input matrices
	!OUTPUT:
	!- S_KRON		: (REAL*4) the output matrix
	!====================================================================
		! INPUT ARGUMENTS
		REAL*4, DIMENSION(:,:), INTENT(IN) :: mat1,mat2
		! OUTPUT
		REAL*4, DIMENSION(SIZE(mat1,1)*SIZE(mat2,1) , SIZE(mat1,2)*SIZE(mat2,2)) :: S_KRON
		! UTILITY VARIABLES
		INTEGER*4 :: ii,jj,rows1,cols1,rows2,cols2
		! FUNCTION DEFINITION
		rows1=SIZE(mat1,1)
		cols1=SIZE(mat1,2)
		rows2=SIZE(mat2,1)
		cols2=SIZE(mat2,2)
		DO ii=1,rows1
			DO jj=1,cols1
				S_KRON((ii-1)*rows2+1:ii*rows2,(jj-1)*cols2+1:jj*cols2) = mat1(ii,jj)*mat2
			END DO
		END DO
		RETURN
	END FUNCTION S_KRON



	FUNCTION D_KRON(mat1, mat2)
	!====================================================================
	!Computes the double precision Kronecker product of the given double precision matrices.
	!INPUT:
	!- mat1, mat2	: (REAL*8) the input matrices
	!OUTPUT:
	!- S_KRON		: (REAL*8) the output matrix
	!====================================================================
		! INPUT ARGUMENTS
		REAL*8, DIMENSION(:,:), INTENT(IN) :: mat1,mat2
		! OUTPUT
		REAL*8, DIMENSION(SIZE(mat1,1)*SIZE(mat2,1) , SIZE(mat1,2)*SIZE(mat2,2)) :: D_KRON
		! UTILITY VARIABLES
		INTEGER*4 :: ii,jj,rows1,cols1,rows2,cols2
		! FUNCTION DEFINITION
		rows1=SIZE(mat1,1)
		cols1=SIZE(mat1,2)
		rows2=SIZE(mat2,1)
		cols2=SIZE(mat2,2)
		DO ii=1,rows1
			DO jj=1,cols1
				D_KRON((ii-1)*rows2+1:ii*rows2,(jj-1)*cols2+1:jj*cols2) = mat1(ii,jj)*mat2
			END DO
		END DO
		RETURN
	END FUNCTION D_KRON



	FUNCTION SVEC_KRON(vec1, vec2)
	!====================================================================
	!Computes the single precision Kronecker product of the given single precision vectors.
	!INPUT:
	!- vec1, vec2	: (REAL*4) the input vectors
	!OUTPUT:
	!- SVEC_KRON	: (REAL*4) the output vector
	!====================================================================
		! INPUT ARGUMENTS
		REAL*4, DIMENSION(:), INTENT(IN) :: vec1,vec2
		! OUTPUT
		REAL*4, DIMENSION(SIZE(vec1)*SIZE(vec2)) :: SVEC_KRON
		! UTILITY VARIABLES
		INTEGER*4 :: ii,jj,nele1,nele2
		! FUNCTION DEFINITION
		nele1=SIZE(vec1)
		nele2=SIZE(vec2)
		DO ii=1,nele1
			SVEC_KRON((ii-1)*nele2+1:ii*nele2) = vec1(ii)*vec2
		END DO
		RETURN
	END FUNCTION SVEC_KRON



	FUNCTION DVEC_KRON(vec1, vec2)
	!====================================================================
	!Computes the double precision Kronecker product of the given double precision vectors.
	!INPUT:
	!- vec1, vec2	: (REAL*8) the input vectors
	!OUTPUT:
	!- DVEC_KRON	: (REAL*8) the output vector
	!====================================================================
		! INPUT ARGUMENTS
		REAL*8, DIMENSION(:), INTENT(IN) :: vec1,vec2
		! OUTPUT
		REAL*8, DIMENSION(SIZE(vec1)*SIZE(vec2)) :: DVEC_KRON
		! UTILITY VARIABLES
		INTEGER*4 :: ii,jj,nele1,nele2
		! FUNCTION DEFINITION
		nele1=SIZE(vec1)
		nele2=SIZE(vec2)
		DO ii=1,nele1
			DVEC_KRON((ii-1)*nele2+1:ii*nele2) = vec1(ii)*vec2
		END DO
		RETURN
	END FUNCTION DVEC_KRON



	FUNCTION S_KHRAO(mat1, mat2)
	!====================================================================
	!Computes the single precision Khatri-Rao product of the given single precision matrices.
	!INPUT:
	!- mat1, mat2	: (REAL*4) the input matrices
	!OUTPUT:
	!- S_KHRAO		: (REAL*4) the output matrix
	!====================================================================
		! INPUT ARGUMENTS
		REAL*4, DIMENSION(:,:), INTENT(IN) :: mat1,mat2
		! OUTPUT
		REAL*4, DIMENSION(SIZE(mat1,1)*SIZE(mat2,1) , SIZE(mat1,2)*SIZE(mat2,2)) :: S_KHRAO
		! UTILITY VARIABLES
		INTEGER*4 :: ii,jj,rows1,cols1,rows2,cols2
		! FUNCTION DEFINITION		
		rows1=SIZE(mat1,1)
		cols1=SIZE(mat1,2)
		rows2=SIZE(mat2,1)
		cols2=SIZE(mat2,2)
		IF (cols1.EQ.cols2) THEN
			DO ii=1,cols1
				S_KHRAO(:,ii)=SVEC_KRON(mat1(:,ii),mat2(:,ii))
			END DO
		ELSE
			WRITE(*,*) "ERROR (S_KHRAO): expected equal second dimension in input matrices. Found ", cols1, cols2
			S_KHRAO = 0E0
		END IF
		RETURN
	END FUNCTION S_KHRAO


	FUNCTION D_KHRAO(mat1, mat2)
	!====================================================================
	!Computes the double precision Khatri-Rao product of the given double precision matrices.
	!INPUT:
	!- mat1, mat2	: (REAL*8) the input matrices
	!OUTPUT:
	!- D_KHRAO		: (REAL*8) the output matrix
	!====================================================================
		! INPUT ARGUMENTS
		REAL*8, DIMENSION(:,:), INTENT(IN) :: mat1,mat2
		! OUTPUT
		REAL*8, DIMENSION(SIZE(mat1,1)*SIZE(mat2,1) , SIZE(mat1,2)*SIZE(mat2,2)) :: D_KHRAO
		! UTILITY VARIABLES
		INTEGER*4 :: ii,jj,rows1,cols1,rows2,cols2
		! FUNCTION DEFINITION		
		rows1=SIZE(mat1,1)
		cols1=SIZE(mat1,2)
		rows2=SIZE(mat2,1)
		cols2=SIZE(mat2,2)
		IF (cols1.EQ.cols2) THEN
			DO ii=1,cols1
				D_KHRAO(:,ii)=DVEC_KRON(mat1(:,ii),mat2(:,ii))
			END DO
		ELSE
			WRITE(*,*) "ERROR (D_KHRAO): expected equal second dimension in input matrices. Found ", cols1, cols2
			D_KHRAO = 0D0
		END IF
		RETURN
	END FUNCTION D_KHRAO



	FUNCTION S_HADAMARD(mat1, mat2)
	!====================================================================
	!Computes the single precision Hadamard product of the given single precision matrices.
	!INPUT:
	!- mat1, mat2	: (REAL*4) the input matrices
	!OUTPUT:
	!- S_HADAMARD	: (REAL*4) the output matrix
	!====================================================================
		! INPUT ARGUMENTS
		REAL*4, DIMENSION(:,:), INTENT(IN) :: mat1,mat2
		! OUTPUT
		REAL*4, DIMENSION(SIZE(mat1,1),SIZE(mat2,2)) :: S_HADAMARD
		! UTILITY VARIABLES
		INTEGER*4 :: ii,jj,rows1,cols1,rows2,cols2
		! FUNCTION DEFINITION		
		rows1=SIZE(mat1,1)
		cols1=SIZE(mat1,2)
		rows2=SIZE(mat2,1)
		cols2=SIZE(mat2,2)
		IF ((cols1.EQ.cols2).AND.(rows1.EQ.rows2)) THEN
			S_HADAMARD = mat1*mat2
		ELSE
			WRITE(*,*) "ERROR (S_HADAMARD): expected equal dimensions in input matrices. Found (", rows1,cols1, ") and (", rows2,cols2, ")"
			S_HADAMARD = 0E0
		END IF
		RETURN
	END FUNCTION S_HADAMARD

	

	FUNCTION D_HADAMARD(mat1, mat2)
	!====================================================================
	!Computes the double precision Hadamard product of the given double precision matrices.
	!INPUT:
	!- mat1, mat2	: (REAL*8) the input matrices
	!OUTPUT:
	!- S_HADAMARD	: (REAL*8) the output matrix
	!====================================================================
		! INPUT ARGUMENTS
		REAL*8, DIMENSION(:,:), INTENT(IN) :: mat1,mat2
		! OUTPUT
		REAL*8, DIMENSION(SIZE(mat1,1),SIZE(mat2,2)) :: D_HADAMARD
		! UTILITY VARIABLES
		INTEGER*4 :: ii,jj,rows1,cols1,rows2,cols2
		! FUNCTION DEFINITION		
		rows1=SIZE(mat1,1)
		cols1=SIZE(mat1,2)
		rows2=SIZE(mat2,1)
		cols2=SIZE(mat2,2)
		IF ((cols1.EQ.cols2).AND.(rows1.EQ.rows2)) THEN
			D_HADAMARD = mat1*mat2
		ELSE
			WRITE(*,*) "ERROR (D_HADAMARD): expected equal dimensions in input matrices. Found (", rows1,cols1, ") and (", rows2,cols2, ")"
			D_HADAMARD = 0E0
		END IF
		RETURN
	END FUNCTION D_HADAMARD



!	FUNCTION S_MODE_PRODUCT(tensor1,mm,tensor2)
!	!====================================================================
!	!Computes the single precision mm-mode product of the given single precision tensors.
!	!INPUT:
!	!- tensor1, tensor2	: (REAL*4) the input tensors
!	!OUTPUT:
!	!- S_MODE_PRODUCT	: (REAL*4) the output tensor
!	!====================================================================
!		! INPUT ARGUMENTS
!		REAL*8, DIMENSION(:,:), INTENT(IN) :: tensor1,tensor2
!		! OUTPUT
!		REAL*8, DIMENSION(SIZE(mat1,1),SIZE(mat2,2)) :: D_HADAMARD
!		! UTILITY VARIABLES
!		INTEGER*4 :: ii,jj,rows1,cols1,rows2,cols2
!		! FUNCTION DEFINITION		
!		rows1=SIZE(mat1,1)
!		cols1=SIZE(mat1,2)
!		rows2=SIZE(mat2,1)
!		cols2=SIZE(mat2,2)
!		IF ((cols1.EQ.cols2).AND.(rows1.EQ.rows2)) THEN
!			D_HADAMARD = mat1*mat2
!		ELSE
!			WRITE(*,*) "ERROR (D_HADAMARD): expected equal dimensions in input matrices. Found (", rows1,cols1, ") and (", rows2,cols2, ")"
!			D_HADAMARD = 0E0
!		END IF
!		RETURN
!	END FUNCTION S_MODE_PRODUCT


	FUNCTION SMODE(tensor,NN)
	!====================================================================	
	!This function takes a stensor and a number, and 
	!returns a version of the tensor in its NN-mode.
	!INPUT:
	!- tensor	: (STENSOR) the input tensors
	!- NN		: (INTEGER*4) the output tensor	
	!OUTPUT:
	!- SMODE	: (REAL*4) the NN-mode representation of the tensor
	!====================================================================
		! INPUT ARGUMENTS
		TYPE(STENSOR), INTENT(IN) :: tensor
		INTEGER*4, INTENT(IN) :: NN
		! OUTPUT
		REAL*4, DIMENSION(tensor%modes(NN),PRODUCT(tensor%modes)/tensor%modes(NN)) :: SMODE
		! UTILITY VARIABLES
		INTEGER*4 :: ii
		! FUNCTION DEFINITION		
		SMODE=0D0
	END FUNCTION SMODE

	
	
!	FUNCTION FROM_FILE(filename)
!	!====================================================================	
!	!This function takes the MNIST file and turns it into an array.
!	!INPUT:
!	!- filename		: (STRING) the input file name	
!	!OUTPUT:
!	!- FROM_FILE	: (REAL*4) the array with the MNIST numbers
!	!====================================================================
!		! INPUT ARGUMENTS
!		CHARACTER(LEN=*) :: filename
!		! OUTPUT
!		REAL*4, DIMENSION(:,:), ALLOCATABLE :: READ_MNIST
!		! UTILITY VARIABLES
!		INTEGER*4 :: ii, jj, numrows
!		CHARACTER(LEN=28*28) :: tmp	
!		! FUNCTION DEFINITION
!		! select if the file is train or test
!		IF (train) THEN
!			numrows=60000
!		ELSE
!			numrows=10000
!		END IF
!		! allocate the appropriate memory
!		ALLOCATE(READ_MNIST(28*28,numrows))
!		! start reading the file
!		OPEN(UNIT=9, FILE=filename)
!		DO ii=1,numrows
!			READ(9, *) READ_MNIST(:,ii)
!!			DO jj=1,28*28		
!!				READ_MNIST(jj,ii)=tmp(jj)
!!			END DO
!		END DO
!		RETURN
!	END FUNCTION READ_MNIST
!	

		



END MODULE MAT_PRODS


















