MODULE MAT_UTILS
  
  ! IN THIS MODULE:
  !- KRONECKER PRODUCT
  !- KHATRI-RAO PRODUCT
  !- HADAMARD PRODUCT
  !- MATRIX AND VECTOR VISUALIZATION
  !- MATRIX RANDOM INITIALIZATION
  !- CUSTOM MATMUL
  !- SVD
  !- MOORE-PENROSE INVERSE
  !- MATRIX-DIAGONAL VECTOR MUTIPLICATION
  !- COLUMN-NORMALIZATION

  
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

  INTERFACE SHOW
     MODULE PROCEDURE SHOW_M
     MODULE PROCEDURE SHOW_V
  END INTERFACE SHOW

  
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
    REAL*4, DIMENSION(SIZE(mat1,1)*SIZE(mat2,1) , SIZE(mat1,2)) :: S_KHRAO
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
    REAL*8, DIMENSION(SIZE(mat1,1)*SIZE(mat2,1) , SIZE(mat1,2)) :: D_KHRAO
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
       WRITE(*,*) "ERROR (S_HADAMARD): expected equal dimensions in input matrices. Found ", rows1,cols1, rows2,cols2
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
       WRITE(*,*) "ERROR (D_HADAMARD): expected equal dimensions in input matrices. Found ", rows1,cols1, rows2,cols2
       D_HADAMARD = 0E0
    END IF
    RETURN
  END FUNCTION D_HADAMARD


  !----------------------SHOW MATRIX SUBROUTINE
  SUBROUTINE SHOW_M(FF,AA)	!Write array A in row,column order.
    INTEGER::FF	!Output file unit number.
    REAL*8::AA(:,:)	!The 2-D array, lower bound one.
    INTEGER::RR	!The row stepper.
    DO RR = 1,SIZE(AA,1)	!Each row gets its own line.
       WRITE (FF,1) AA(RR,:)		!Write all the columns of that row.
1      FORMAT (*(G10.3:x))		!This suffices for the example.
    END DO			!On to the next row.
  END SUBROUTINE SHOW_M	!WRITE (F,*) A or similar would show A as if transposed.


  !----------------------SHOW VECTOR SUBROUTINE
  SUBROUTINE SHOW_V(FF,AA)
    INTEGER::FF
    REAL*8::AA(:)
    WRITE (FF,1) AA
1   FORMAT (*(G10.3:x))
  END SUBROUTINE SHOW_V


  !----------------------GAUSSIAN RANDOM INITIALIZATION SUBROUTINE
  SUBROUTINE RAND_INIT(mat,sigma)
    REAL*8, ALLOCATABLE::mat(:,:)
    REAL*8::uni_1,uni_2,norm_1,norm_2,sigma
    INTEGER::ii,jj,dim1,dim2,new_dim1
    dim1=SIZE(mat,1)
    dim2=SIZE(mat,2)
    new_dim1=dim1
    IF (MOD(new_dim1,2).EQ.1) THEN
       new_dim1=new_dim1-1
    END IF
    DO jj=1,dim2
       DO ii=1,new_dim1,2
          CALL RANDOM_NUMBER(uni_1)
          CALL RANDOM_NUMBER(uni_2)
          mat(ii,jj)  =sigma*SQRT(-2*LOG(uni_1))*COS(8.D0*ATAN(1.D0)*uni_2)
          mat(ii+1,jj)=sigma*SQRT(-2*LOG(uni_1))*SIN(8.D0*ATAN(1.D0)*uni_2)
       END DO
       IF (MOD(dim1,2).EQ.1) THEN
          CALL RANDOM_NUMBER(uni_1)
          CALL RANDOM_NUMBER(uni_2)
          mat(dim1,jj)=sigma*SQRT(-2*LOG(uni_1))*COS(8.D0*ATAN(1.D0)*uni_2)
       END IF
    END DO
  END SUBROUTINE RAND_INIT


  SUBROUTINE SVD(AA,UU,SIG,VVT,INFO)
    !====================================================================
    !Returns the SVD of the matrix AA.
    !INPUT/OUTPUT:
    !- AA 		: (REAL*8) the input matrix (size M*N)
    !- UU 		: (REAL*8) the left singular vectors matrix (size M*M)
    !- SIG              : (REAL*8) the singular values vector (size min(M,N))
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


  FUNCTION MTDG(matrix,diagonal)
!!$===========================================================================
!!$    This function computes the product of a matrix and a diagonal matrix,
!!$    stored as an array.
!!$    INPUT ARGUMENTS:
!!$    - matrix    : (REAL*8) the input matrix 
!!$    - norms     : (REAL*8) the input vector which represents the
!!$                  diagonal matrix
!!$    OUTPUT:
!!$    - mtdg      : (REAL*8) the result of the matrix-matrix multiplication
!!$===========================================================================
    ! INOUT VARIABLES
    REAL*8 :: matrix(:,:), diagonal(:)
    REAL*8 :: mtdg(SIZE(matrix,1),SIZE(matrix,2))
    ! UTILITY VARIABLES
    INTEGER*4 :: ii

    IF (SIZE(matrix,2).NE.SIZE(diagonal)) THEN
       WRITE(*,*) "WARNING (MTDG): sizes do not match!"
    END IF
    DO ii=1,SIZE(matrix,2)
       mtdg(:,ii)=matrix(:,ii)*diagonal(ii)
    END DO
    RETURN
  END FUNCTION MTDG


  
  SUBROUTINE COL_NORM(matrix,norms)
!!$===========================================================================
!!$    Normalizes the columns of a given matrix, and stores the columnwise
!!$    norms inside norms.
!!$    INOUT:
!!$    - matrix    : (REAL*8) the input matrix to be normalized
!!$    - norms     : (REAL*8) the vector which contains the norms of the columns
!!$===========================================================================
    ! INOUT VARIABLES
    REAL*8 :: matrix(:,:), norms(:)
    ! UTILITY VARIABLES
    INTEGER*4 :: ii
    
    IF (SIZE(matrix,2).NE.SIZE(norms)) THEN
       WRITE(*,*) "WARNING (COL_NORM): sizes do not match!"
    END IF
    DO ii=1,SIZE(matrix,2)
       norms(ii)=SQRT(SUM(matrix(:,ii)**2))
       matrix(:,ii)=matrix(:,ii)/norms(ii)
    END DO
  END SUBROUTINE COL_NORM




END MODULE MAT_UTILS
