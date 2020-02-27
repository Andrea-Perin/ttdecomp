MODULE CPD_UTILS

  ! MODULES
  USE DEBUGGER
  USE TENSOR_TYPES
  USE MODE_N
  USE MAT_UTILS


  IMPLICIT NONE
! DEBUG FLAG
  LOGICAL :: deb=.TRUE.
! DEBUG FLAG


  ! INTERFACE FOR CPD
  INTERFACE CPD
     MODULE PROCEDURE CPD3
     MODULE PROCEDURE CPD4
  END INTERFACE CPD

  
CONTAINS



  ! ====================================
  ! ====================================
  ! CANONIC POLYADIC DECOMPOSITION
  ! ====================================
  ! ====================================
  

  SUBROUTINE CPD3(tensor, rank, factors, lambdas, error, thresh, numiter, verbose)
!!$==========================================================================================================     
!!$    Computes the Alternating Least Squares algorithm for the CPD decomposition.
!!$    Each factor matrix is updated as
!!$    B^(n) <- X_(n) (KHRAO_{k/=n} B^(k)) (HAD_{k/=n} B^(k)^T B^(k))^{INV}
!!$    After this, the mean squared error is computed; the algorithm goes on as long as the
!!$    decrease in the relative error is larger than some threshold.
!!$    INOUT:
!!$    - tensor       : (REAL*8) the input matrix to be inverted
!!$    - rank         : (INTEGER*4) the rank of the decomposition (nunmber of factors in the linear combination)
!!$    - factors      : (matrix_list) the list of the factor matrices. Should be 3 elements long
!!$    - lambdas      : (REAL*8) a vector containing the lambdas   
!!$    - error        : (REAL*8) the reconstruction error
!!$    - (O) thresh   : (REAL*8, OPTIONAL) the threshold on the error decrease. DEFAULT: 5D-5
!!$    - (O) numiter  : (INTEGER*4, OPTIONAL) the maximum number of iterations. DEFAULT: 500
!!$    - (O) verbose  : (LOGICAL, OPTIONAL) whether to print the error at each iteration. DEFAULT: .FALSE.  
!!$==========================================================================================================
    ! INOUT VARIABLES
    TYPE(DTENSOR3) :: tensor
    INTEGER*4 :: rank
    REAL*8 :: error
    TYPE(matrix_list) :: factors(3)
    REAL*8, ALLOCATABLE :: lambdas(:)
    REAL*8, OPTIONAL :: thresh
    LOGICAL, OPTIONAL :: verbose
    INTEGER*4, OPTIONAL :: numiter
    !UTILITY VARIABLES
    INTEGER*4 :: maxiter=500
    REAL*8 :: threshold=5D-5
    INTEGER*4 :: ii, jj, NN=SIZE(factors), idx(2), krsize, cnt 
    REAL*8, ALLOCATABLE :: V(:,:), B(:,:)
    REAL*8, ALLOCATABLE :: tensor_rec(:,:), B_rec(:,:)
    REAL*8 :: newerror, relative_error

    !========================================================================================================
    ! SET THE MAXIMUM NUMBER OF ITERATIONS
    IF (PRESENT(numiter).AND.(numiter.GE.1)) THEN
       maxiter = numiter
    END IF
    IF (PRESENT(thresh).AND.(thresh.GE.1D-15)) THEN
       threshold = thresh
    END IF
    ! ALLOCATE THE VECTOR OF LAMBDAS
    ALLOCATE(lambdas(rank))
    ! STORE THE DIMENSIONS OF THE TENSOR, ALLOCATE THE FACTOR MATRICES AND NORMALIZE THEM
    DO ii=1,NN
       ALLOCATE(factors(ii)%matr(tensor%modes(ii),rank))
       CALL RANDOM_NUMBER(factors(ii)%matr)
       CALL COL_NORM(factors(ii)%matr,lambdas)
    END DO
    ! ALLOCATE FOR V
    ALLOCATE(V(rank,rank))
    ! ALLOCATE THE RECONSTRUCTION STUFF
    ALLOCATE(tensor_rec(tensor%modes(1),PRODUCT(tensor%modes)/tensor%modes(1)))
    ALLOCATE(B_rec(PRODUCT(tensor%modes(2:)),rank))    
    ! REPEAT UNTIL CONVERGENCE (IN RELATIVE ERROR) OR EXCEEDED ITERATIONS
    cnt=0
    relative_error=1
    DO WHILE ((relative_error.GT.threshold).AND.(cnt.LT.maxiter))
       cnt=cnt+1
       ! LOOP OVER THE FACTOR MATRICES
       DO ii=1,NN
          ! DETERMINE THE INDICES TO SKIP
          idx = PACK( (/(jj, jj=1,NN,1)/) , (/(jj, jj=1,NN,1)/)/=ii )
          ! COMPUTE V
          V=1D0
          DO jj=1,SIZE(idx)
             V = V.HAD.MTML(TRANSPOSE(factors(idx(jj))%matr),factors(idx(jj))%matr)
          END DO
          ! COMPUTE B
          ALLOCATE(B(PRODUCT( (/(tensor%modes(idx(jj)), jj=1,SIZE(idx),1)/) ),rank))
          krsize = tensor%modes(idx(1))*tensor%modes(idx(2))
          B(1:krsize,:) = factors(idx(2))%matr.KHRAO.factors(idx(1))%matr
          DO jj=3,SIZE(idx)
             B(1:krsize*tensor%modes(idx(jj)),:) = factors(idx(jj))%matr.KHRAO.B(1:krsize,:)
             krsize=krsize*tensor%modes(idx(jj))
          END DO
          print*, "RANK SHOULD BE", rank
          print*, "Shape of V", SHAPE(V)
          print*, "Shape of B", SHAPE(B)
          print*, "Chosen indices", idx
          print*, "last columns of B", SUM(B(:,SIZE(B,2)))
          ! UPDATE THE ii-th FACTOR MATRIX 
          factors(ii)%matr = MTML( MTML( (tensor.MODE.ii),B ), PINV(V) )
          CALL COL_NORM(factors(ii)%matr, lambdas)
          ! DEALLOCATE B
          DEALLOCATE(B)
       END DO
       ! RECONSTRUCT THE TENSOR
       krsize = tensor%modes(2)*tensor%modes(3)
       B_rec(1:krsize,:)=factors(3)%matr.KHRAO.factors(2)%matr
       DO ii=4,NN
          B_rec(1:krsize*tensor%modes(ii),:)=factors(ii)%matr.KHRAO.B_rec(1:krsize,:)
          krsize=krsize*tensor%modes(ii)
       END DO
       tensor_rec = MTML(MTDG(factors(1)%matr,lambdas), TRANSPOSE(B_rec)) 
       ! COMPUTE THE ERROR (FROBENIUS NORM)
       newerror = SQRT(SUM(((tensor.MODE.1)-tensor_rec)**2))/SIZE(tensor%elems)
       relative_error = ABS((error-newerror)/error)
       IF (PRESENT(verbose).AND.verbose) THEN
          PRINT*, cnt, newerror, relative_error
       END IF
       error=newerror
    END DO
  END SUBROUTINE CPD3


  
  SUBROUTINE CPD4(tensor, rank, factors, lambdas, error, thresh, numiter, verbose)
!!$==========================================================================================================     
!!$    Computes the Alternating Least Squares algorithm for the CPD decomposition.
!!$    Each factor matrix is updated as
!!$    B^(n) <- X_(n) (KHRAO_{k/=n} B^(k)) (HAD_{k/=n} B^(k)^T B^(k))^{INV}
!!$    After this, the mean squared error is computed; the algorithm goes on as long as the
!!$    decrease in the relative error is larger than some threshold.
!!$    INOUT:
!!$    - tensor       : (REAL*8) the input matrix to be inverted
!!$    - rank         : (INTEGER*4) the rank of the decomposition (nunmber of factors in the linear combination)
!!$    - factors      : (matrix_list) the list of the factor matrices. Should be 4 elements long
!!$    - lambdas      : (REAL*8) a vector containing the lambdas   
!!$    - error        : (REAL*8) the reconstruction error
!!$    - (O) thresh   : (REAL*8, OPTIONAL) the threshold on the error decrease. DEFAULT: 5D-5
!!$    - (O) numiter  : (INTEGER*4, OPTIONAL) the maximum number of iterations. DEFAULT: 500
!!$    - (O) verbose  : (LOGICAL, OPTIONAL) whether to print the error at each iteration. DEFAULT: .FALSE.  
!!$==========================================================================================================
    ! INOUT VARIABLES
    TYPE(DTENSOR4) :: tensor
    INTEGER*4 :: rank
    REAL*8 :: error
    TYPE(matrix_list) :: factors(4)
    REAL*8, ALLOCATABLE :: lambdas(:)
    REAL*8, OPTIONAL :: thresh
    LOGICAL, OPTIONAL :: verbose
    INTEGER*4, OPTIONAL :: numiter
    !UTILITY VARIABLES
    INTEGER*4 :: maxiter=500
    REAL*8 :: threshold=5D-5
    INTEGER*4 :: ii, jj, NN=SIZE(factors), idx(3), krsize, cnt 
    REAL*8, ALLOCATABLE :: V(:,:), B(:,:)
    REAL*8, ALLOCATABLE :: tensor_rec(:,:), B_rec(:,:)
    REAL*8 :: newerror, relative_error

    !========================================================================================================
    ! SET THE MAXIMUM NUMBER OF ITERATIONS
    IF (PRESENT(numiter).AND.(numiter.GE.1)) THEN
       maxiter = numiter
    END IF
    IF (PRESENT(thresh).AND.(thresh.GE.1D-15)) THEN
       threshold = thresh
    END IF
    ! ALLOCATE THE VECTOR OF LAMBDAS
    ALLOCATE(lambdas(rank))
    ! STORE THE DIMENSIONS OF THE TENSOR, ALLOCATE THE FACTOR MATRICES AND NORMALIZE THEM
    DO ii=1,NN
       ALLOCATE(factors(ii)%matr(tensor%modes(ii),rank))
       CALL RANDOM_NUMBER(factors(ii)%matr)
       CALL COL_NORM(factors(ii)%matr,lambdas)
    END DO
    ! ALLOCATE THE RECONSTRUCTION STUFF
    ALLOCATE(tensor_rec(tensor%modes(1),PRODUCT(tensor%modes)/tensor%modes(1)))
    ALLOCATE(B_rec(PRODUCT(tensor%modes(2:)),rank))    
    ! REPEAT UNTIL CONVERGENCE (IN RELATIVE ERROR) OR EXCEEDED ITERATIONS
    cnt=0
    relative_error=1
    DO WHILE ((relative_error.GT.threshold).AND.(cnt.LT.maxiter))
       cnt=cnt+1
       ! LOOP OVER THE FACTOR MATRICES
       DO ii=1,NN
          ! DETERMINE THE INDICES TO SKIP
          idx = PACK( (/(jj, jj=1,NN,1)/) , (/(jj, jj=1,NN,1)/)/=ii )
          ! COMPUTE V
          ALLOCATE(V(rank,rank))
          V=1D0
          DO jj=1,SIZE(idx)
             V = V.HAD.MTML(TRANSPOSE(factors(idx(jj))%matr),factors(idx(jj))%matr)
          END DO
          ! COMPUTE B
          ALLOCATE(B(PRODUCT( (/(tensor%modes(idx(jj)), jj=1,SIZE(idx),1)/) ),rank))
          krsize = tensor%modes(idx(1))*tensor%modes(idx(2))
          B(1:krsize,:) = factors(idx(2))%matr.KHRAO.factors(idx(1))%matr
          DO jj=3,SIZE(idx)
             B(1:krsize*tensor%modes(idx(jj)),:) = factors(idx(jj))%matr.KHRAO.B(1:krsize,:)
             krsize=krsize*tensor%modes(idx(jj))
          END DO
          ! UPDATE THE ii-th FACTOR MATRIX 
          factors(ii)%matr = MTML( MTML( (tensor.MODE.ii),B ), PINV(V) )
          CALL COL_NORM(factors(ii)%matr, lambdas)
          ! DEALLOCATE B AND V
          DEALLOCATE(B)
          DEALLOCATE(V)
       END DO
       ! RECONSTRUCT THE TENSOR
       krsize = tensor%modes(2)*tensor%modes(3)
       B_rec(1:krsize,:)=factors(3)%matr.KHRAO.factors(2)%matr
       DO ii=4,NN
          B_rec(1:krsize*tensor%modes(ii),:)=factors(ii)%matr.KHRAO.B_rec(1:krsize,:)
          krsize=krsize*tensor%modes(ii)
       END DO
       tensor_rec = MTML(MTDG(factors(1)%matr,lambdas), TRANSPOSE(B_rec)) 
       ! COMPUTE THE ERROR (FROBENIUS NORM)
       newerror = SQRT(SUM(((tensor.MODE.1)-tensor_rec)**2))/SIZE(tensor%elems)
       relative_error = ABS((error-newerror)/error)
       IF (PRESENT(verbose).AND.verbose) THEN
          PRINT*, cnt, newerror, relative_error
       END IF
       error=newerror
    END DO
  END SUBROUTINE CPD4


  
  ! ====================================
  ! ====================================
  ! TO IDENTITY
  ! ====================================
  ! ====================================

  FUNCTION TO_IDENTITY(vector, order)
    ! ======================================================
    ! This function takes a vector, containing the lambdas
    ! of the CPD decomposition, and turns it into the
    ! mode-1 representation of the corresponding core.
    ! From there, it can be easily
    ! turned into a tensor.
    ! INPUT ARGUMENTS
    ! - vector       : (REAL*8) the vector of the lambdas.
    ! - order        : (INTEGER*4) the order of the tensor
    ! Its size is equal to the rank.
    ! OUTPUT ARGUMENTS
    ! - TO_IDENTITY  : (REAL*8) the mode-1 representation
    !                  of the diagonal tensor of lambdas.
    ! ======================================================    
    ! INOUT VARIABLES
    REAL*8 :: vector(:)
    INTEGER*4 :: order
    REAL*8, ALLOCATABLE :: TO_IDENTITY(:,:)
    ! UTILITY VARIABLES
    INTEGER*4 :: ii, RR
    
    RR = SIZE(vector)
    ALLOCATE(TO_IDENTITY(RR,RR**(order-1)))
    TO_IDENTITY = 0D0
    DO ii=1,RR
       TO_IDENTITY(ii,RR*(ii-1)+ii) = vector(ii)
    END DO
    RETURN
  END FUNCTION TO_IDENTITY




  ! ====================================
  ! ====================================
  ! new cpd maybe
  ! ====================================
  ! ====================================

  SUBROUTINE NEWCPD3(tensor, rank, factors, lambdas, error, thresh, numiter, verbose)
    TYPE(DTENSOR3) :: tensor
    INTEGER*4 :: rank
    TYPE(MATRIX_LIST) :: factors(3)
    REAL*8, ALLOCATABLE :: lambdas(:)
    REAL*8 :: error
    REAL*8, OPTIONAL :: thresh
    INTEGER*4, OPTIONAL :: numiter
    LOGICAL, OPTIONAL :: verbose
    ! UTILITY VARIABLES
    INTEGER*4 :: ii,jj,rel_err=1,cnt=0,NN=SIZE(factors)
    INTEGER*4 :: idx(SIZE(factors)-1)
    REAL*8 :: tol=5D-5
    INTEGER*4 :: maxiter=100
    REAL*8, ALLOCATABLE :: V(:,:), B(:,:)
    REAL*8, ALLOCATABLE :: B_rec(:,:), tens_rec(:,:)
    REAL*8 :: newerror, relative_error
    
    ! FUNCTION BODY
    ! set optional parameters
    IF (PRESENT(numiter).AND.(numiter.GE.1)) THEN
       maxiter = numiter
    END IF
    IF (PRESENT(thresh).AND.(thresh.GE.1D-15)) THEN
       tol = thresh
    END IF
    ! allocate the vector of lambdas
    ALLOCATE(lambdas(rank))
    ! randomly initialize the factor matrices
    DO ii=1,NN
       ! factor matrices are (I_n)x(R)
       ALLOCATE( factors(ii)%matr(tensor%modes(ii),rank) )
       ! random init
       CALL RANDOM_NUMBER(factors(ii)%matr)
       ! and normalize just in case
       CALL COL_NORM(factors(ii)%matr,lambdas)
    END DO
    ! allocate V, B and B_rec, tens_rec.
    ! - V has always the same shape, rank X rank
    ! - B changes shape. using the properties of allocatable variables,
    !   it is allocated only once, and reshape/reset when needed.
    !   Initially, it is (I_2 x I_3 x ... x I_N) X (rank)
    ! - B_rec is used to compute the reconstruciton error.
    !   It is the cumulative khatri rao product from 2 onwards, so always the same shape
    ! - tens_rec stores the reconstruction, and has always the same shape
    ALLOCATE(V(rank,rank))
    ALLOCATE(B(PRODUCT(tensor%modes(2:)),rank))
    ALLOCATE(B_rec(PRODUCT(tensor%modes(2:)),rank))
    ALLOCATE(tens_rec(tensor%modes(1), PRODUCT(tensor%modes(2:))))
    ! entry condition: relative error and number of iterations
    DO WHILE ((rel_err.GT.tol).AND.(cnt.LT.maxiter))
       ! update the iterations counter
       cnt=cnt+1
       ! now, loop over the factor matrices
       DO ii=1,NN
          ! create an array where the needed indices are stored
          idx = PACK( (/(jj,jj=1,NN,1)/), (/(jj,jj=1,NN,1)/)/=1 )
          ! compute the matrix V
          V = 1D0
          DO jj=1,NN-1
             V=V.HAD.MTML(TRANSPOSE(factors(idx(jj))%matr),factors(idx(jj))%matr)
          END DO
          ! compute the matrix B, which is a cumulative khatri rao
          ! the first passage is just the KR product of the first two matrices
          B = RESHAPE(B, (/ tensor%modes(idx(1))*tensor%modes(idx(2)), rank /) )
          B = factors(idx(2))%matr.KHRAO.factors(idx(1))%matr
          DO jj=3,NN-1
             B = factors(idx(jj))%matr.KHRAO.B
          END DO
          ! now, store the results in the right factor matrix
          factors(ii)%matr = MTML( tensor.MODE.ii, MTML( B,PINV(V) ) )
          ! normalize the columns of the factor matrix, and store the norms in lambdas
          CALL COL_NORM(factors(ii)%matr,lambdas)         
       END DO
       ! now, compute the error on the reconstruction
       ! first, compute the transposed cumulative khatri rao
       B_rec = RESHAPE(B_rec, (/tensor%modes(3)*tensor%modes(2),rank/) )
       B_rec = factors(3)%matr.KHRAO.factors(2)%matr
       DO ii=4,NN
          B_rec = factors(ii)%matr.KHRAO.B_rec
       END DO
       ! then combine everything 
       tens_rec = MTML( factors(1)%matr.MDDOT.lambdas, TRANSPOSE(B))
       newerror = SQRT(SUM( ((tensor.MODE.1)-tens_rec)**2 ))/SIZE(tensor%elems)
       ! and compute the relative error
       relative_error = ABS((error-newerror)/error)
       IF (PRESENT(verbose).AND.verbose) THEN
          PRINT*, cnt, newerror, relative_error
       END IF
       error=newerror
    END DO
  END SUBROUTINE NEWCPD3


  
  
END MODULE CPD_UTILS
