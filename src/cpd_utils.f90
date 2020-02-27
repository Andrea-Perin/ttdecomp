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



  SUBROUTINE CPD3(tensor, rank, factors, lambdas, error, thresh, numiter, verbose, svdinit)
!!$==========================================================================================================     
!!$    Computes the Alternating Least Squares algorithm for the CPD decomposition.
!!$    Each factor matrix is updated as
!!$    B^(n) <- X_(n) (KHRAO_{k/=n} B^(k)) (HAD_{k/=n} B^(k)^T B^(k))^{INV}
!!$    After this, the mean squared error is computed; the algorithm goes on as long as the
!!$    decrease in the relative error is larger than some threshold.
!!$    INOUT:
!!$    - tensor       : (DTENSOR3) the input tensor
!!$    - rank         : (INTEGER*4) the rank of the decomposition (nunmber of factors in the linear combination)
!!$    - factors      : (matrix_list) the list of the factor matrices. Should be 3 elements long
!!$    - lambdas      : (REAL*8) a vector containing the lambdas   
!!$    - error        : (REAL*8) the reconstruction error
!!$    - (O) thresh   : (REAL*8, OPTIONAL) the threshold on the error decrease. DEFAULT: 1D-5
!!$    - (O) numiter  : (INTEGER*4, OPTIONAL) the maximum number of iterations. DEFAULT: 100
!!$    - (O) verbose  : (LOGICAL, OPTIONAL) whether to print the error at each iteration. DEFAULT: .FALSE.
!!$    - (O) svdinit  : (LOGICAL, OPTIONAL) whether to initialize the factors with svd. DEFAULT: .FALSE.    
!!$==========================================================================================================
    TYPE(DTENSOR3) :: tensor
    INTEGER*4 :: rank
    TYPE(MATRIX_LIST) :: factors(3)
    REAL*8, ALLOCATABLE :: lambdas(:)
    REAL*8 :: error
    REAL*8, OPTIONAL :: thresh
    INTEGER*4, OPTIONAL :: numiter
    LOGICAL, OPTIONAL :: verbose
    LOGICAL, OPTIONAL :: svdinit
    ! UTILITY VARIABLES
    INTEGER*4 :: ii,jj,rel_err=1,cnt=0,NN=SIZE(factors), INFO
    INTEGER*4 :: idx(SIZE(factors)-1)
    REAL*8 :: tol=1D-5
    INTEGER*4 :: maxiter=100
    REAL*8, ALLOCATABLE :: V(:,:), B(:,:)
    REAL*8, ALLOCATABLE :: B_rec(:,:), tens_rec(:,:)
    REAL*8 :: newerror, relative_error
    REAL*8, ALLOCATABLE :: UU(:,:), SIG(:), VVT(:,:)
    
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
    ! initialize the factor matrices
    DO ii=1,NN
       ! factor matrices are (I_n)x(R)
       ALLOCATE( factors(ii)%matr(tensor%modes(ii),rank) )
       ! choose if random or svd init
       IF (PRESENT(svdinit).AND.svdinit) THEN          
          ! TSVD init
          IF (ALLOCATED(UU)) DEALLOCATE(UU)
          IF (ALLOCATED(SIG)) DEALLOCATE(SIG)
          IF (ALLOCATED(VVT)) DEALLOCATE(VVT)
          CALL TSVD(tensor.MODE.ii, UU, SIG, VVT, INFO)
          IF (rank.GT.tensor%modes(ii)) THEN
             factors(ii)%matr(:,1:tensor%modes(ii)) = UU
             CALL RANDOM_NUMBER(factors(ii)%matr(:,tensor%modes(ii):)) 
          ELSE
             factors(ii)%matr = UU(:,1:rank)
          END IF
       ELSE          
          ! random init
          CALL RANDOM_NUMBER(factors(ii)%matr)
       END IF
       ! and normalize just in  case 
       CALL COL_NORM(factors(ii)%matr,lambdas)
    END DO
    ! allocate V, B and B_rec, tens_rec.
    ! - V has always the same shape, rank X rank
    ! - B changes shape. Initially, it is (I_2 x I_3 x ... x I_N) X (rank)
    ! - B_rec is used to compute the reconstruciton error.
    !   It is the cumulative khatri rao product from 2 onwards, so always the same shape
    ! - tens_rec stores the reconstruction, and has always the same shape
    ALLOCATE(V(rank,rank))
    ALLOCATE(B_rec(PRODUCT(tensor%modes(2:)),rank))
    ALLOCATE(tens_rec(tensor%modes(1), PRODUCT(tensor%modes(2:))))
    ! entry condition: relative error and number of iterations
    DO WHILE ((rel_err.GT.tol).AND.(cnt.LT.maxiter))
       ! update the iterations counter
       cnt=cnt+1
       ! now, loop over the factor matrices
       DO ii=1,NN
          ! create an array where the needed indices are stored
          idx = PACK( (/(jj,jj=1,NN,1)/), (/(jj,jj=1,NN,1)/)/=ii )
          ! compute the matrix V
          V = 1D0
          DO jj=1,SIZE(idx)
             V=V.HAD.MTML(TRANSPOSE(factors(idx(jj))%matr),factors(idx(jj))%matr)
          END DO
          ! compute the matrix B, which is a cumulative khatri rao
          ! its size is the following
          ALLOCATE(B( PRODUCT(tensor%modes(idx(:))), rank ))
          ! the first passage is just the KR product of the first two matrices
          B(1:PRODUCT(tensor%modes(idx(:2))),:) = factors(idx(2))%matr.KHRAO.factors(idx(1))%matr
          DO jj=3,SIZE(idx)
             B(1:PRODUCT(tensor%modes(idx(:jj))),:) = factors(idx(jj))%matr.KHRAO.B(1:PRODUCT(tensor%modes(idx(:jj-1))),:)
          END DO
          ! now, store the results in the right factor matrix
          factors(ii)%matr = MTML( tensor.MODE.ii, MTML( B,PINV(V) ) )
          DEALLOCATE(B)
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
       tens_rec = MTML( factors(1)%matr,  lambdas.MDDOT.TRANSPOSE(B_rec))
       newerror = SQRT(SUM( ((tensor.MODE.1)-tens_rec)**2 ))/SIZE(tensor%elems)
       ! and compute the relative error
       relative_error = ABS((error-newerror)/error)
       IF (PRESENT(verbose).AND.verbose) THEN
          PRINT*, cnt, newerror, relative_error
       END IF
       error=newerror
    END DO
  END SUBROUTINE CPD3



  SUBROUTINE CPD4(tensor, rank, factors, lambdas, error, thresh, numiter, verbose, svdinit)
!!$==========================================================================================================     
!!$    Computes the Alternating Least Squares algorithm for the CPD decomposition.
!!$    Each factor matrix is updated as
!!$    B^(n) <- X_(n) (KHRAO_{k/=n} B^(k)) (HAD_{k/=n} B^(k)^T B^(k))^{INV}
!!$    After this, the mean squared error is computed; the algorithm goes on as long as the
!!$    decrease in the relative error is larger than some threshold.
!!$    INOUT:
!!$    - tensor       : (DTENSOR4) the input tensor to be inverted
!!$    - rank         : (INTEGER*4) the rank of the decomposition (number of factors in the linear combination)
!!$    - factors      : (matrix_list) the list of the factor matrices. Should be 4 elements long
!!$    - lambdas      : (REAL*8) a vector containing the lambdas   
!!$    - error        : (REAL*8) the reconstruction error
!!$    - (O) thresh   : (REAL*8, OPTIONAL) the threshold on the error decrease. DEFAULT: 1D-5
!!$    - (O) numiter  : (INTEGER*4, OPTIONAL) the maximum number of iterations. DEFAULT: 100
!!$    - (O) verbose  : (LOGICAL, OPTIONAL) whether to print the error at each iteration. DEFAULT: .FALSE.
!!$    - (O) svdinit  : (LOGICAL, OPTIONAL) whether to initialize the factors with svd. DEFAULT: .FALSE.    
!!$==========================================================================================================
    TYPE(DTENSOR4) :: tensor
    INTEGER*4 :: rank
    TYPE(MATRIX_LIST) :: factors(4)
    REAL*8, ALLOCATABLE :: lambdas(:)
    REAL*8 :: error
    REAL*8, OPTIONAL :: thresh
    INTEGER*4, OPTIONAL :: numiter
    LOGICAL, OPTIONAL :: verbose
    LOGICAL, OPTIONAL :: svdinit
    ! UTILITY VARIABLES
    INTEGER*4 :: ii,jj,rel_err=1,cnt=0,NN=SIZE(factors), INFO
    INTEGER*4 :: idx(SIZE(factors)-1)
    REAL*8 :: tol=1D-5
    INTEGER*4 :: maxiter=100
    REAL*8, ALLOCATABLE :: V(:,:), B(:,:)
    REAL*8, ALLOCATABLE :: B_rec(:,:), tens_rec(:,:)
    REAL*8 :: newerror, relative_error
    REAL*8, ALLOCATABLE :: UU(:,:), SIG(:), VVT(:,:)
    
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
    ! initialize the factor matrices
    DO ii=1,NN
       ! factor matrices are (I_n)x(R)
       ALLOCATE( factors(ii)%matr(tensor%modes(ii),rank) )
       ! choose if random or svd init
       IF (PRESENT(svdinit).AND.svdinit) THEN          
          ! TSVD init
          IF (ALLOCATED(UU)) DEALLOCATE(UU)
          IF (ALLOCATED(SIG)) DEALLOCATE(SIG)
          IF (ALLOCATED(VVT)) DEALLOCATE(VVT)
          CALL TSVD(tensor.MODE.ii, UU, SIG, VVT, INFO)
          IF (rank.GT.tensor%modes(ii)) THEN
             factors(ii)%matr(:,1:tensor%modes(ii)) = UU
             CALL RANDOM_NUMBER(factors(ii)%matr(:,tensor%modes(ii):)) 
          ELSE
             factors(ii)%matr = UU(:,1:rank)
          END IF
       ELSE          
          ! random init
          CALL RANDOM_NUMBER(factors(ii)%matr)
       END IF
       ! and normalize just in  case 
       CALL COL_NORM(factors(ii)%matr,lambdas)
    END DO
    ! allocate V, B and B_rec, tens_rec.
    ! - V has always the same shape, rank X rank
    ! - B changes shape. Initially, it is (I_2 x I_3 x ... x I_N) X (rank)
    ! - B_rec is used to compute the reconstruciton error.
    !   It is the cumulative khatri rao product from 2 onwards, so always the same shape
    ! - tens_rec stores the reconstruction, and has always the same shape
    ALLOCATE(V(rank,rank))
    ALLOCATE(B_rec(PRODUCT(tensor%modes(2:)),rank))
    ALLOCATE(tens_rec(tensor%modes(1), PRODUCT(tensor%modes(2:))))
    ! entry condition: relative error and number of iterations
    DO WHILE ((rel_err.GT.tol).AND.(cnt.LT.maxiter))
       ! update the iterations counter
       cnt=cnt+1
       ! now, loop over the factor matrices
       DO ii=1,NN
          ! create an array where the needed indices are stored
          idx = PACK( (/(jj,jj=1,NN,1)/), (/(jj,jj=1,NN,1)/)/=ii )
          ! compute the matrix V
          V = 1D0
          DO jj=1,SIZE(idx)
             V=V.HAD.MTML(TRANSPOSE(factors(idx(jj))%matr),factors(idx(jj))%matr)
          END DO
          ! compute the matrix B, which is a cumulative khatri rao
          ! its size is the following
          ALLOCATE(B( PRODUCT(tensor%modes(idx(:))), rank ))
          ! the first passage is just the KR product of the first two matrices
          B(1:PRODUCT(tensor%modes(idx(:2))),:) = factors(idx(2))%matr.KHRAO.factors(idx(1))%matr
          DO jj=3,SIZE(idx)
             B(1:PRODUCT(tensor%modes(idx(:jj))),:) = factors(idx(jj))%matr.KHRAO.B(1:PRODUCT(tensor%modes(idx(:jj-1))),:)
          END DO
          ! now, store the results in the right factor matrix
          factors(ii)%matr = MTML( tensor.MODE.ii, MTML( B,PINV(V) ) )
          DEALLOCATE(B)
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
       tens_rec = MTML( factors(1)%matr,  lambdas.MDDOT.TRANSPOSE(B_rec))
       newerror = SQRT(SUM( ((tensor.MODE.1)-tens_rec)**2 ))/SIZE(tensor%elems)
       ! and compute the relative error
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




  
  
END MODULE CPD_UTILS
