MODULE TUCKER
  
  ! IN THIS MODULE:
  USE DEBUGGER
  USE MODE_N
  USE MAT_UTILS
  USE TENSOR_TYPES
  

  IMPLICIT NONE

  ! INTERFACE FOR HOOI
  INTERFACE HOOI
     MODULE PROCEDURE HOOI3
     MODULE PROCEDURE HOOI4
  END INTERFACE HOOI

  ! INTERFACE FOR HOSVD
  INTERFACE HOSVD
     MODULE PROCEDURE HOSVD3
     MODULE PROCEDURE HOSVD4
  END INTERFACE HOSVD

  ! INTERFACE FOR TUCKER-RELATED CORE COMPUTATION
  INTERFACE TCORE
     MODULE PROCEDURE TCORE3
     MODULE PROCEDURE TCORE4
  END INTERFACE TCORE

  ! INTERFACE FOR TUCKER-RELATED RECONSTRUCTION
  INTERFACE RECO
     MODULE PROCEDURE RECO3
     MODULE PROCEDURE RECO4
  END INTERFACE RECO


CONTAINS


  !======================================================= 
  !======================================================= 
  ! RECOSTRUCTIONS
  !======================================================= 
  !=======================================================
  
  FUNCTION RECO3(core, factors)
    !=================================================================
    ! This function computes the reconstruction of the full
    ! tensor that is decomposed with the Tucker decomposition.
    ! The full tensor is returned.
    ! INPUT
    ! - core     : (DTENSOR3) the core tensor
    ! - factors  : (MATRIX_LIST) the list of factor matrices
    !              (in this case, Tucker factors)
    ! OUTPUT
    ! - RECO3    : (DTENSOR3) the full tensor
    !=================================================================
    TYPE(DTENSOR3) :: core,RECO3
    TYPE(MATRIX_LIST) :: factors(3)
    ! UTILITY VARIABLES
    REAL*8, ALLOCATABLE :: Xhat(:,:), Xtilde(:,:)
    INTEGER*4 :: ranks(SIZE(factors)), newmodes(SIZE(factors))
    INTEGER*4 :: NN=SIZE(factors), ii


    ! FILL THE RANKS
    DO ii=1,NN
       ranks(ii) = SIZE(factors(ii)%matr,1)
    END DO
    ! From now on:
    ! - Xhat contains the mode_n representation of the enhanced core
    ! - Xtilde contains the matrix product of Xhat with a factor matrix
    ! In the beginning,
    newmodes = core%modes
    ALLOCATE(Xhat(core%modes(1), PRODUCT(core%modes(2:))))
    Xhat=core.MODE.1
    DO ii=1,NN-1
       ! Update nemwodes: the ii-th entry becomes rank(ii)
       newmodes(ii)=ranks(ii)
       ! In Xhat2, store the matrix product of the mode_n tensor and
       ! the factor matrix
       ALLOCATE(Xtilde(newmodes(ii), PRODUCT(newmodes)/newmodes(ii)))
       Xtilde = MTML(factors(ii)%matr, Xhat)
       DEALLOCATE(Xhat)
       ALLOCATE(Xhat(newmodes(ii+1),PRODUCT(newmodes)/newmodes(ii+1) ))
       Xhat = TENSOR3(newmodes,Xtilde,ii).MODE.(ii+1)
       ! Now, Xhat has become the n-mode representation of the running
       ! mode-n product.
       ! Xtilde can be deallocated.
       DEALLOCATE(Xtilde)
    END DO
    ! The last Xtilde will contain the mode_NN representation
    ! of the tensor
    newmodes(NN)=ranks(NN)
    ALLOCATE(Xtilde(newmodes(NN),PRODUCT(newmodes(1:NN-1))) )
    Xtilde = MTML(factors(NN)%matr, Xhat)
    RECO3 = TENSOR3(newmodes,Xtilde,NN)
    RETURN
  END FUNCTION RECO3


  FUNCTION RECO4(core, factors)
    !=================================================================
    ! This function computes the reconstruction of the full tensor
    ! that is decomposed with the Tucker decomposition.
    ! The full tensor is returned.
    ! INPUT
    ! - core     : (DTENSOR4) the core tensor
    ! - factors  : (MATRIX_LIST) the list of factor matrices
    !              (in this case, Tucker factors)
    ! OUTPUT
    ! - RECO4    : (DTENSOR4) the full tensor
    !=================================================================
    TYPE(DTENSOR4) :: core,RECO4
    TYPE(MATRIX_LIST) :: factors(4)
    ! UTILITY VARIABLES
    REAL*8, ALLOCATABLE :: Xhat(:,:), Xtilde(:,:)
    INTEGER*4 :: ranks(SIZE(factors)), newmodes(SIZE(factors))
    INTEGER*4 :: NN=SIZE(factors), ii


    ! FILL THE RANKS
    DO ii=1,NN
       ranks(ii) = SIZE(factors(ii)%matr,1)
    END DO
    ! From now on:
    ! - Xhat contains the mode_n representation of the enhanced core
    ! - Xtilde contains the matrix product of Xhat with a factor matrix
    ! In the beginning,
    newmodes = core%modes
    ALLOCATE(Xhat(core%modes(1), PRODUCT(core%modes(2:))))
    Xhat=core.MODE.1
    DO ii=1,NN-1
       ! Update nemwodes: the ii-th entry becomes rank(ii)
       newmodes(ii)=ranks(ii)
       ! In Xhat2, store the matrix product of the mode_n tensor and
       ! the factor matrix
       ALLOCATE(Xtilde(newmodes(ii), PRODUCT(newmodes)/newmodes(ii)))
       Xtilde = MTML(factors(ii)%matr, Xhat)
       DEALLOCATE(Xhat)
       ALLOCATE(Xhat(newmodes(ii+1),PRODUCT(newmodes)/newmodes(ii+1) ))
       Xhat = TENSOR4(newmodes,Xtilde,ii).MODE.(ii+1)
       ! Now, Xhat has become the n-mode representation of the running
       ! mode-n product.
       ! Xtilde can be deallocated.
       DEALLOCATE(Xtilde)
    END DO
    ! The last Xtilde will contain the mode_NN representation
    ! of the tensor
    newmodes(NN)=ranks(NN)
    ALLOCATE(Xtilde(newmodes(NN),PRODUCT(newmodes(1:NN-1))) )
    Xtilde = MTML(factors(NN)%matr, Xhat)
    RECO4 = TENSOR4(newmodes,Xtilde,NN)
    RETURN
  END FUNCTION RECO4

  
  !======================================================= 
  !======================================================= 
  ! TUCKER CORES
  !======================================================= 
  !======================================================= 


  FUNCTION TCORE3(tensor, factors) RESULT(Xtilde)
    !=================================================================
    ! This function computes the reconstruction of a core tensor that
    ! is employed in the Tucker decomposition. The NN-mode (in this
    ! case, 3) of the core tensor is returned.
    ! If the core tensor is needed, just employ the proper TENSOR
    ! function.
    ! INPUT
    ! - tensor     : (DTENSOR3) the input tensor
    ! - factors    : (MATRIX_LIST) the list of factor matrices
    !                (in this case, Tucker factors)
    ! OUTPUT
    ! - Xtilde     : (REAL*8) the matrix of size
    !                       (R_NN,prod_{n=1}^NN-1 R_i)
    !=================================================================
    TYPE(DTENSOR3) :: tensor
    TYPE(MATRIX_LIST) :: factors(3)
    ! UTILITY VARIABLES
    REAL*8, ALLOCATABLE :: Xhat(:,:), Xtilde(:,:)
    INTEGER*4 :: ranks(SIZE(factors)), newmodes(SIZE(factors))
    INTEGER*4 :: NN=SIZE(factors), ii

    ! FILL THE RANKS
    DO ii=1,NN
       ranks(ii) = SIZE(factors(ii)%matr,2)
    END DO
    ! From now on:
    ! - Xhat contains the mode_n representation of the tensor
    ! - Xtilde contains the matrix product of Xhat with a factor matrix
    ! In the beginning,
    newmodes = tensor%modes
    ALLOCATE(Xhat(tensor%modes(1), PRODUCT(tensor%modes(2:))))
    Xhat=tensor.MODE.1
    DO ii=1,NN-1
       ! Update nemwodes: the ii-th entry becomes rank(ii)
       newmodes(ii)=ranks(ii)
       ! In Xhat2, store the matrix product of the mode_n tensor and
       ! the factor matrix
       ALLOCATE(Xtilde(newmodes(ii), PRODUCT(newmodes)/newmodes(ii)))
       Xtilde = MTML(TRANSPOSE(factors(ii)%matr), Xhat)
       DEALLOCATE(Xhat)
       ALLOCATE(Xhat(newmodes(ii+1),PRODUCT(newmodes)/newmodes(ii+1) ))
       Xhat = TENSOR3(newmodes,Xtilde,ii).MODE.(ii+1)
       ! Now, Xhat has become the n-mode representation of the running
       ! mode-n product.
       ! Xtilde can be deallocated.
       DEALLOCATE(Xtilde)
    END DO
    ! The last Xtilde will contain the mode_NN representation
    ! of the tensor
    newmodes(NN)=ranks(NN)
    ALLOCATE(Xtilde(newmodes(NN),PRODUCT(newmodes(1:NN-1))) )
    Xtilde = MTML(TRANSPOSE(factors(NN)%matr), Xhat)
    RETURN
  END FUNCTION TCORE3

  
  FUNCTION TCORE4(tensor, factors) RESULT(Xtilde)
    !=================================================================
    ! This function computes the reconstruction of a core tensor that
    ! is employed in the Tucker decomposition. The NN-mode (in this
    ! case, 4) of the core tensor is returned. If the core tensor is
    ! needed, just employ the proper TENSOR(N) function.
    ! INPUT
    ! - tensor    : (DTENSOR4) the input tensor
    ! - factors   : (MATRIX_LIST) the list of factor matrices
    !               (in this case, Tucker factors)
    ! OUTPUT
    ! - Xtilde    : (REAL*8) the matrix of size
    !                      (R_NN, prod_{n=1}^NN-1 R_i)
    !=================================================================
    TYPE(DTENSOR4) :: tensor
    TYPE(MATRIX_LIST) :: factors(4)
    ! UTILITY VARIABLES
    REAL*8, ALLOCATABLE :: Xhat(:,:), Xtilde(:,:)
    INTEGER*4 :: ranks(SIZE(factors)), newmodes(SIZE(factors))
    INTEGER*4 :: NN=SIZE(factors), ii

    ! FILL THE RANKS
    DO ii=1,NN
       ranks(ii) = SIZE(factors(ii)%matr,2)
    END DO
    ! From now on:
    ! - Xhat contains the mode_n representation of the tensor
    ! - Xtilde contains the matrix product of Xhat with a factor matrix
    ! In the beginning,
    newmodes = tensor%modes
    ALLOCATE(Xhat(tensor%modes(1), PRODUCT(tensor%modes(2:))))
    Xhat=tensor.MODE.1
    DO ii=1,NN-1
       ! Update nemwodes: the ii-th entry becomes rank(ii)
       newmodes(ii)=ranks(ii)
       ! In Xhat2, store the matrix product of the mode_n tensor and
       ! the factor matrix
       ALLOCATE(Xtilde(newmodes(ii), PRODUCT(newmodes)/newmodes(ii)))
       Xtilde = MTML(TRANSPOSE(factors(ii)%matr), Xhat)
       DEALLOCATE(Xhat)
       ALLOCATE(Xhat(newmodes(ii+1),PRODUCT(newmodes)/newmodes(ii+1) ))
       Xhat = TENSOR4(newmodes,Xtilde,ii).MODE.(ii+1)
       ! Now, Xhat has become the n-mode representation of the
       ! running mode-n product.
       ! Xtilde can be deallocated.
       DEALLOCATE(Xtilde)
    END DO
    ! The last Xtilde will contain the mode_NN representation
    ! of the tensor
    newmodes(NN)=ranks(NN)
    ALLOCATE(Xtilde(newmodes(NN),PRODUCT(newmodes(1:NN-1))) )
    Xtilde = MTML(TRANSPOSE(factors(NN)%matr), Xhat)
    RETURN
  END FUNCTION TCORE4


  !======================================================= 
  !======================================================= 
  ! HOSVD
  !======================================================= 
  !======================================================= 
    
  
  SUBROUTINE HOSVD3(tens,ranks,core,factors)
    !=================================================================
    !Returns core and factors of the Tucker Decomposition using HOSVD.
    !INPUT/OUTPUT:
    !- tens 		: (DTENSOR3) the input tensor (size S1,S2,S3)
    !- ranks 		: (INTEGER*4) the vector of ranks (size 3) 
    !- core             : (DTENSOR3) the core tensor (size S1,S2,S3)
    !- factors 		: (MATRIX_LIST) vector (size 3) of factor
    !                     matrices (size Si,Ri)
    !=================================================================
    ! INOUT VARIABLES
    TYPE(DTENSOR3) :: tens,core
    INTEGER*4 :: ranks(3)
    TYPE(MATRIX_LIST) :: factors(3)
    ! UTILITY VARIABLES
    INTEGER*4 :: ii, info, NN=SIZE(factors)
    INTEGER*4 :: new(SIZE(factors)) ! new modes 
    REAL*8, ALLOCATABLE :: SIG(:), UU(:,:), VVT(:,:), res(:,:)
    ! ALLOCATE THE FACTORS
    DO ii=1,NN
       ALLOCATE(factors(ii)%matr(tens%modes(ii),ranks(ii)))
    END DO
    ! PERFORM SVD ON ALL POSSIBLE MODES
    DO ii=1,NN
       IF (ALLOCATED(UU)) DEALLOCATE(UU)
       IF (ALLOCATED(SIG)) DEALLOCATE(SIG)
       IF (ALLOCATED(VVT)) DEALLOCATE(VVT)
       CALL SVD(tens.MODE.ii,UU,SIG,VVT,info)
       factors(ii)%matr=UU(:,1:ranks(ii))
    END DO
    ! ALLOCATE CORE TENSOR
    core%modes=ranks
    ALLOCATE(core%elems(ranks(1),ranks(2),ranks(3)))
    ! COMPUTE THE CORE TENSOR
    res = TCORE(tens,factors) ! TCORE returns the nn mode
    core = TENSOR3(ranks,res,NN) ! go from nn mode to tensor
  END SUBROUTINE HOSVD3

  
  SUBROUTINE HOSVD4(tens,ranks,core,factors)
    !=================================================================
    !Returns core and factors of the Tucker Decomposition using HOSVD.
    !INPUT/OUTPUT:
    !- tens 		: (DTENSOR4) the input tensor (size S1,S2,S3,S4)
    !- ranks 		: (INTEGER*4) the vector of ranks (size 4) 
    !- core             : (DTENSOR4) the core tensor (size S1,S2,S3,S4)
    !- factors 		: (MATRIX_LIST) vector (size 4) of factor
    !                     matrices (size Si,Ri)
    !=================================================================
    ! INOUT VARIABLES
    TYPE(DTENSOR4) :: tens,core
    INTEGER*4 :: ranks(4)
    TYPE(MATRIX_LIST) :: factors(4)
    ! UTILITY VARIABLES
    INTEGER*4 :: ii, info, NN=SIZE(factors)
    INTEGER*4 :: new(SIZE(factors)) ! new modes 
    REAL*8, ALLOCATABLE :: SIG(:), UU(:,:), VVT(:,:), res(:,:)
    ! ALLOCATE THE FACTORS
    DO ii=1,NN
       ALLOCATE(factors(ii)%matr(tens%modes(ii),ranks(ii)))
    END DO
    ! PERFORM SVD ON ALL POSSIBLE MODES
    DO ii=1,NN
       IF (ALLOCATED(UU)) DEALLOCATE(UU)
       IF (ALLOCATED(SIG)) DEALLOCATE(SIG)
       IF (ALLOCATED(VVT)) DEALLOCATE(VVT)
       CALL SVD(tens.MODE.ii,UU,SIG,VVT,info)
       factors(ii)%matr=UU(:,1:ranks(ii))
    END DO
    ! ALLOCATE CORE TENSOR
    core%modes=ranks
    ALLOCATE(core%elems(ranks(1),ranks(2),ranks(3),ranks(4)))
    ! COMPUTE THE CORE TENSOR
    res = TCORE(tens,factors) ! TCORE returns the nn mode
    core = TENSOR4(ranks,res,NN) ! go from nn mode to tensor
  END SUBROUTINE HOSVD4
  

  
  !======================================================= 
  !======================================================= 
  ! HOOI
  !======================================================= 
  !======================================================= 
  

  SUBROUTINE HOOI3(tensor, ranks, core, factors, error, verbose, numiter, thresh)
    !=================================================================
    !Returns core and factors of the Tucker Decomposition using HOOI.
    !INPUT/OUTPUT:
    !- tensor 		: (DTENSOR3) the input tensor
    !- ranks 		: (INTEGER*4) the vector of ranks  
    !- core             : (DTENSOR3) the core tensor that is
    !                     retrieved by HOOI
    !- factors 		: (MATRIX_LIST) the list of factor matrices
    !                     (size Si,Ri)
    !- error            : (REAL*8, OPTIONAL) the reconstruction error
    !                     (normalized on the size of the tensor) 
    !- verbose          : (LOGICAL, OPTIONAL) whether to print the
    !                     error on screen at each iteration
    !- numiter          : (INTEGER*4, OPTIONAL) the maximum number of
    !                     iterations
    !- thresh           : (REAL*8, OPTIONAL) the threshold on the
    !                     relative error decrease
    !=================================================================
    ! INOUT VARIABLES
    TYPE(DTENSOR3) :: tensor, core
    TYPE(MATRIX_LIST) :: factors(3)
    INTEGER*4 :: ranks(3)
    REAL*8 :: error
    REAL*8, OPTIONAL :: thresh
    INTEGER*4, OPTIONAL :: numiter
    LOGICAL, OPTIONAL :: verbose
    ! UTILITY VARIABLES
    !matrices for mode-n product
    REAL*8, ALLOCATABLE :: Xhat(:,:), Xtilde(:,:)
    !matrices for SVD
    REAL*8, ALLOCATABLE :: UU(:,:), VV(:,:), SIGMA(:)   
    INTEGER*4 :: ii, jj, NN=SIZE(factors), cnt, idx(SIZE(factors)), INFO
    INTEGER*4 :: newmodes(SIZE(factors))
    REAL*8 :: relative_error, newerror
    INTEGER*4 :: maxiter=500
    REAL*8 :: threshold=5D-6
    
    ! ACTUAL FUNCTION
    ! SET OPTIONAL PARAMETERS
    IF (PRESENT(thresh).AND.(thresh.GT.1D-15)) THEN
       threshold=thresh
    END IF
    IF (PRESENT(numiter).AND.(numiter.GT.1)) THEN
       maxiter=numiter
    END IF   
    ! INITIALIZE THE FACTOR MATRICES WITH HOSVD
    CALL HOSVD(tensor,ranks,core,factors)
    ! REPEAT UNTIL CONVERGENCE
    !do ii=1,NN
       !print*, SHAPE(factors(ii)%matr)
    !end do
    cnt=0
    relative_error=1
    DO WHILE ((relative_error.GT.threshold).AND.(cnt.LT.maxiter))
       cnt=cnt+1
       DO ii=1,NN
          ! DETERMINE THE INDICES TO SKIP
          idx = PACK( (/(jj, jj=1,NN,1)/) , (/(jj, jj=1,NN,1)/)/=ii )
          idx(NN) = ii
          ! KEEP TRACK OF THE CHANGES IN THE MODES
          newmodes = tensor%modes
          ! PERFORM THE MODE-N PRODUCTS
          ! WITHOUT SAVING ANY HIGH ORDER TENSOR
          ! 
          ! From now on:
          ! - Xhat1 will store the tensor in one of its
          !   mode_n representations.
          ! - Xhat2 will store the matrix product of Xhat1
          !   with one of the factor matrices.
          !
          ALLOCATE(Xhat(newmodes(idx(1)),PRODUCT(newmodes)/newmodes(idx(1))))
          Xhat=tensor.MODE.idx(1)
          DO jj=1,NN-1
             !
             ! IMPORTANT: recall that we are reindexing, using idx,
             !            in order to skip the ii-th index.
             !            As such, all quantities are indexed by
             !            idx(jj), and not by jj!
             !
             ! On Xtilde, the matrix product of the transposed
             ! factor matrix and the mode_n tensor is stored.
             ! As such, it has size (R_n, W), where W denotes
             ! the current product of the other dimensions.
             ! Recall that, after each mode_n product, the n-th
             ! dimension is decreased, from I_n->R_n.
             ! The newmodes variable keeps track of these changes
             ! in dimension.
             newmodes(idx(jj)) = ranks(idx(jj))
             ALLOCATE(Xtilde(newmodes(idx(jj)), PRODUCT(newmodes)/newmodes(idx(jj))) )
             Xtilde=MTML(TRANSPOSE(factors(idx(jj))%matr), Xhat)
             ! Now that the matrix product has been computed, Xhat
             ! can be deallocated and reallocated. It now serves
             ! a different purpose: it must store the tensor in
             ! the next mode, so its shape has to change accordingly.
             ! This basically means that it will become a matrix
             ! of shape (I_{n}, W), where W is, as before, the
             ! product of the "running sizes".
             DEALLOCATE(Xhat)
             ALLOCATE(Xhat(newmodes(idx(jj+1)), PRODUCT(newmodes)/newmodes(idx(jj+1))))
             Xhat=TENSOR3(newmodes,Xtilde,idx(jj)).MODE.idx(jj+1)
             ! Now that Xhat has become, once again, the matrix
             ! representation of the "current iteration" of the tensor,
             ! Xtilde can be safely deallocated, in order to prepare
             ! it to the next iteration, where it will contain
             ! the matrix product of the next step. 
             DEALLOCATE(Xtilde)
          END DO
          ! GET THE FIRST R_N LEADING LEFT SINGULAR VECTORS OF Y_ii
          ! (which has not been actually created!)
          IF (ALLOCATED(UU)) DEALLOCATE(UU)
          IF (ALLOCATED(SIGMA)) DEALLOCATE(SIGMA)
          IF (ALLOCATED(VV)) DEALLOCATE(VV)
          !print*, "dc", SIZE(XhAt,1)
          CALL SVD(Xhat,UU,SIGMA,VV,INFO)
          factors(ii)%matr=UU(:,1:ranks(ii))
          DEALLOCATE(Xhat)
       END DO
       !RECONSTRUCT THE TENSOR AND COMPUTE THE ERROR
       ! Compute the error 
       newerror = SQRT(SUM(tensor%elems**2)-SUM(TCORE(tensor,factors)**2))/SIZE(tensor%elems)
       relative_error = ABS((error-newerror)/error)
       error = newerror
       ! If verbose is on, output the current error
       IF (PRESENT(verbose).AND.verbose) THEN
          WRITE(*,*) cnt, error, relative_error
       END IF
    END DO
    ! COMPUTE THE CORE TENSOR
    core = TENSOR3(ranks,TCORE(tensor,factors),NN)
  END SUBROUTINE HOOI3


  SUBROUTINE HOOI4(tensor, ranks, core, factors, error, verbose, numiter, thresh)
    !=================================================================
    !Returns core and factors of the Tucker Decomposition using HOOI.
    !INPUT/OUTPUT:
    !- tensor 		: (DTENSOR4) the input tensor
    !- ranks 		: (INTEGER*4) the vector of ranks  
    !- core             : (DTENSOR4) the core tensor that is retrieved
    !                     by HOOI
    !- factors 		: (MATRIX_LIST) the list of factor matrices
    !                     (size Si,Ri)
    !- error            : (REAL*8, OPTIONAL) the reconstruction error
    !                     (normalized on the size of the tensor) 
    !- verbose          : (LOGICAL, OPTIONAL) whether to print the
    !                     error on screen at each iteration
    !- numiter          : (INTEGER*4, OPTIONAL) the maximum number
    !                     of iterations
    !- thresh           : (REAL*8, OPTIONAL) the threshold on the
    !                     relative error decrease
    !=================================================================
    ! INOUT VARIABLES
    TYPE(DTENSOR4) :: tensor, core
    TYPE(MATRIX_LIST) :: factors(4)
    INTEGER*4 :: ranks(4)
    REAL*8 :: error
    REAL*8, OPTIONAL :: thresh
    INTEGER*4, OPTIONAL :: numiter
    LOGICAL, OPTIONAL :: verbose
    ! UTILITY VARIABLES
    !matrices for mode-n product
    REAL*8, ALLOCATABLE :: Xhat(:,:), Xtilde(:,:)
    !matrices for SVD
    REAL*8, ALLOCATABLE :: UU(:,:), VV(:,:), SIGMA(:)   
    INTEGER*4 :: ii, jj, NN=SIZE(factors), cnt, idx(SIZE(factors)), INFO
    INTEGER*4 :: newmodes(SIZE(factors))
    REAL*8 :: relative_error, newerror
    INTEGER*4 :: maxiter=500
    REAL*8 :: threshold=5D-6
    
    ! ACTUAL FUNCTION
    ! SET OPTIONAL PARAMETERS
    IF (PRESENT(thresh).AND.(thresh.GT.1D-15)) THEN
       threshold=thresh
    END IF
    IF (PRESENT(numiter).AND.(numiter.GT.1)) THEN
       maxiter=numiter
    END IF   
    ! INITIALIZE THE FACTOR MATRICES WITH HOSVD
    CALL HOSVD(tensor,ranks,core,factors)
    ! REPEAT UNTIL CONVERGENCE
    cnt=0
    relative_error=1
    DO WHILE ((relative_error.GT.threshold).AND.(cnt.LT.maxiter))
       cnt=cnt+1
       DO ii=1,NN
          ! DETERMINE THE INDICES TO SKIP
          idx = PACK( (/(jj, jj=1,NN,1)/) , (/(jj, jj=1,NN,1)/)/=ii )
          idx(NN) = ii
          ! KEEP TRACK OF THE CHANGES IN THE MODES
          newmodes = tensor%modes
          ! PERFORM THE MODE-N PRODUCTS
          ! WITHOUT SAVING ANY HIGH ORDER TENSOR
          ! 
          ! From now on:
          ! - Xhat1 will store the tensor in one of its
          !   mode_n representations.
          ! - Xhat2 will store the matrix product of Xhat1
          !   with one of the factor matrices.
          !
          ALLOCATE(Xhat(newmodes(idx(1)),PRODUCT(newmodes)/newmodes(idx(1))))
          Xhat=tensor.MODE.idx(1)
          DO jj=1,NN-1
             !
             ! IMPORTANT: recall that we are reindexing, using idx,
             ! in order to skip the ii-th index. As such, all
             ! quantities are indexed by idx(jj), and not by jj!
             !
             ! On Xtilde, the matrix product of the transposed
             ! factor matrix and the mode_n tensor is stored.
             ! As such, it has size (R_n, W), where W denotes
             ! the current product of the other dimensions.
             ! Recall that, after each mode_n product, the n-th
             ! dimension is decreased, from I_n->R_n.
             ! The newmodes variable keeps track of these
             ! changes in dimension.
             newmodes(idx(jj)) = ranks(idx(jj))
             ALLOCATE(Xtilde(newmodes(idx(jj)), PRODUCT(newmodes)/newmodes(idx(jj))) )
             Xtilde=MTML(TRANSPOSE(factors(idx(jj))%matr), Xhat)
             ! Now that the matrix product has been computed, Xhat
             ! can be deallocated and reallocated. It now serves
             ! a different purpose: it must store the tensor in
             ! the next mode, so its shape has to change accordingly.
             ! This basically means that it will become a matrix
             ! of shape (I_{n}, W), where W is, as before, the
             ! product of the "running sizes".
             DEALLOCATE(Xhat)
             ALLOCATE(Xhat(newmodes(idx(jj+1)), PRODUCT(newmodes)/newmodes(idx(jj+1))))
             Xhat=TENSOR4(newmodes,Xtilde,idx(jj)).MODE.idx(jj+1)
             ! Now that Xhat has become, once again, the matrix
             ! representation of the "current iteration" of the tensor,
             ! Xtilde can be safely deallocated, in order to prepare
             ! it to the next iteration, where it will contain
             ! the matrix product of the next step. 
             DEALLOCATE(Xtilde)
          END DO
          ! GET THE FIRST R_N LEADING LEFT SINGULAR VECTORS OF Y_ii
          ! (which has not been actually created!)
          IF (ALLOCATED(UU)) DEALLOCATE(UU)
          IF (ALLOCATED(SIGMA)) DEALLOCATE(SIGMA)
          IF (ALLOCATED(VV)) DEALLOCATE(VV)
          CALL SVD(Xhat,UU,SIGMA,VV,INFO)
          factors(ii)%matr=UU(:,1:ranks(ii))
          DEALLOCATE(Xhat)
       END DO
       !RECONSTRUCT THE TENSOR AND COMPUTE THE ERROR
       ! Compute the error 
       newerror = SQRT(SUM(tensor%elems**2)-SUM(TCORE(tensor,factors)**2))/SIZE(tensor%elems)
       relative_error = ABS((error-newerror)/error)
       error = newerror
       ! If verbose is on, output the current error
       IF (PRESENT(verbose).AND.verbose) THEN
          WRITE(*,*) cnt, error, relative_error
       END IF
    END DO
    ! COMPUTE THE CORE TENSOR
    core = TENSOR4(ranks,TCORE(tensor,factors),NN)
  END SUBROUTINE HOOI4







  ! NON NEGATIVE TUCKER DECOMPOSITION(?)

  ! SUBROUTINE MU_NTD3(tensor, ranks, core, factors)
  !   ! INOUT VARIABLES
  !   TYPE(DTENSOR3) :: tensor
  !   INTEGER*4 :: ranks(3)
  !   TYPE(DTENSOR3) :: core
  !   TYPE(MATRIX_LIST) :: factors(3)
  !   ! UTILITY VARIABLES
  !   TYPE(MATRIX_LIST) :: fac_tilde(3)
  !   TYPE(DTENSOR3) :: core_tilde, XX, XX_tilde
  !   INTEGER*4 :: ii, cnt, cnv_cnt, maxiter=100, NN=SIZE(factors)
  !   REAL*8 :: error, relative_error, threshold=1D-6

  !   ! INITIALIZATION OF TILDED VARIABLES USING HOSVD
  !   CALL HOSVD(tensor,ranks,core_tilde,fac_tilde)
  !   ! INITIALIZATION OF UNTILDED VARIABLES
  !   DO ii=1,NN
  !      ALLOCATE(factors(ii)%matr(tensor%modes(ii), ranks(ii)))
  !      CALL RANDOM_NUMBER(factors(ii)%matr)
  !   END DO
  !   core = TENSOR3()


    
  !   ! CONTINUE UP UNTIL CONVERGENCE
  !   relative_error = 1D0
  !   cnt = 0
  !   DO WHILE ((cnt.LT.maxiter).AND.(relative_error.GT.threshold))
  !      cnt=cnt+1
  !      ! LOOP OVER THE FACTOR MATRICES
  !      DO ii=1,NN
  !         ! UNTIL CONVERGENCE, UPDATE X, XTILDE, A^(ii)
  !         cnv_cnt=0
  !         DO WHILE (().AND.())
  !            cnv_cnt=cnv_cnt+1
             


             
  !         END DO

  !         ! UNTIL CONVERGENCE, UPDATE CORE
  !         DO WHILE (().AND.())

  !         END DO


          
  !      END DO
    
    






  ! END SUBROUTINE MU_NTD3



  
END MODULE TUCKER
