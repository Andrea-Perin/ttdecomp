MODULE PARTINV
	
	IMPLICIT NONE
		
CONTAINS


    SUBROUTINE SVD(AA,UU,SIG,VVT,INFO)
    !====================================================================
	!Returns the SVD of the matrix AA.
	!INPUT/OUTPUT:
	!- AA 		: (REAL*8) the input matrix (size M*N)
	!- UU 		: (REAL*8) the left singular vectors matrix (size M*M)
	!- SIG  	: (REAL*8) the singular values vector (size min(M,N))
	!- VVT 		: (REAL*8) the TRANSPOSED right singular vectors matrix (size N*N)
	!- INFO		: (INTEGER*4) the info about how the subroutine worked out
	!====================================================================
	    ! INOUT VARIABLES
        REAL*8 :: AA(:,:)
        REAL*8, ALLOCATABLE :: UU(:,:), VVT(:,:), SIG(:)
        INTEGER*4 :: INFO
        ! UTILITY VARIABLES
        INTEGER*4 :: ii, MM, NN, LWORK
        REAL*8, ALLOCATABLE :: WORK(:),COPY(:,:)
        ! FUNCTION DEFINITION
        ! ALLOCATE THE MATRICES
        MM = SIZE(AA,1)
        NN = SIZE(AA,2)
        ALLOCATE(UU(MM,MM))
        ALLOCATE(VVT(NN,NN))
        ALLOCATE(SIG(MIN(MM,NN)))
        ALLOCATE(COPY(MM,NN))
        COPY = AA
        ! QUERY WORKSPACE
        LWORK = -1
        ALLOCATE(WORK(1))
        CALL DGESVD('A', 'A', MM, NN, AA, MM, SIG, UU, MM, VVT, NN, WORK, LWORK, INFO)
        LWORK = WORK(1)
        DEALLOCATE(WORK)
        ALLOCATE(WORK(LWORK))
        ! PERFORM SVD
        CALL DGESVD('A', 'A', MM, NN, AA, MM, SIG, UU, MM, VVT, NN, WORK, LWORK, INFO)
        !CHECK RESULTS
		IF (INFO.LT.0) THEN
			WRITE (*,*) 'ERROR (SVD): Illegal argument at position ', -INFO
		ELSEIF (INFO.GT.0) THEN
			WRITE (*,*) 'WARNING (SVD): Failed convergence in ', INFO, &
						'superdiagonals of an intermediate bidiagonal form.'
		END IF
		AA = COPY
    END SUBROUTINE SVD


    

	FUNCTION MTML(AA,BB)
	!====================================================================
	!Performs matrix-matrix multiplication between AA and BB.
	!INPUT:
	!- AA 		: (REAL*8) the input matrix (size M*N)
	!- BB 		: (REAL*8) the input matrix (size N*L)
	!OUTPUT:
	!- MTML		: (REAL*8) the output matrix (size M*L)
	!====================================================================
		! INPUT ARGUMENTS
		REAL*8 :: AA(:,:), BB(:,:)
		!OUTPUT VARIABLES
		REAL*8 :: MTML(SIZE(AA,1),SIZE(BB,2))
		! UTILITY VARIABLES
		INTEGER*4 :: ii, cols1, rows2
		! FUNCTION DEFINITION
		!CHECK DIMENSIONS
		cols1=SIZE(AA,2)
		rows2=SIZE(BB,1)
		IF (cols1.NE.rows2) THEN
			WRITE (*,*) 'ERROR (MTML): expected common dimension. Found ',cols1, rows2
			MTML = 0D0			
		ELSE
			CALL DGEMM( 'N','N',SIZE(AA,1),SIZE(BB,2),cols1,1D0,AA,MAX(1,SIZE(AA,1)), &
						BB,MAX(1,cols1),0D0,MTML,MAX(1,SIZE(AA,1)))
		END IF
		RETURN
	END FUNCTION MTML



	FUNCTION PINV(AA)
	!====================================================================
	!Returns the Moore-Penrose inverse (aka, pseudoinverse) of the matrix AA.
	!INPUT:
	!- AA 		: (REAL*8) the input matrix to be inverted
	!OUTPUT:
	!- PINV 	: (REAL*8) the pseudoinverse of AA
	!====================================================================
		! INPUT ARGUMENTS
		REAL*8 :: AA(:,:)
		! OUTPUT ARGUMENTS		
		REAL*8 :: PINV(SIZE(AA,2),SIZE(AA,1))
		! UTILITY VARIABLES
		INTEGER*4 :: ii, jj, INFO
		REAL*8, ALLOCATABLE :: UU(:,:), VVT(:,:), SIGMA(:)
		REAL*8, ALLOCATABLE :: SIGMA_MAT(:,:)
		REAL*8 :: tol
		! FUNCTION DEFINITION
		!PERFORM SVD ON THE MATRIX AA
		CALL SVD(AA, UU, SIGMA, VVT, INFO)
		IF (INFO.EQ.0) THEN !NO ERRORS: GO ON WITH PINV
			!PSEUDOINVERSE OF RECTANGULAR DIAGONAL MATRIX
			tol = EPSILON(1D0)*MAX(SIZE(AA,1),SIZE(AA,2))*MAXVAL(SIGMA)
			WHERE (SIGMA.LE.tol)
				SIGMA = 0D0
			ELSEWHERE
				SIGMA = 1/SIGMA
			END WHERE
			!TURN THE VECTOR TO MATRIX
			ALLOCATE(SIGMA_MAT(SIZE(AA,2),SIZE(AA,1)))
            SIGMA_MAT = 0D0
			DO ii=1,MINVAL(SHAPE(SIGMA_MAT))
                SIGMA_MAT(ii,ii)=SIGMA(ii)
            END DO
			!COMPOSE THE WHOLE PSEUDO INVERSE
			PINV = MTML(MTML(TRANSPOSE(VVT),SIGMA_MAT), TRANSPOSE(UU))
		ELSE !SOME KIND OF ERROR: ABORT, RETURN TRANSPOSED MATRIX
			WRITE (*,*) 'ERROR (PINV): Encountered an error in SVD computation.'
			PINV=TRANSPOSE(AA)
		END IF
		RETURN
	END FUNCTION PINV


END MODULE PARTINV
