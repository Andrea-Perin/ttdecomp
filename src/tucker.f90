MODULE TUCKER
  
  ! IN THIS MODULE:
  USE DEBUGGER
  USE MODE_N
  USE MAT_UTILS
  USE TENSOR_TYPES
  

  IMPLICIT NONE

  INTERFACE HOSVD
     MODULE PROCEDURE HOSVD3
  END INTERFACE HOSVD


  ! INTERFACE FOR TUCKER-RELATED RECONSTRUCTION
  INTERFACE RECO
     MODULE PROCEDURE RECO3
     MODULE PROCEDURE RECO4
     MODULE PROCEDURE RECO5
  END INTERFACE RECO
     


CONTAINS


  !======================================================= 
  !======================================================= 
  ! RECONSTRUCTIONS
  !======================================================= 
  !======================================================= 
  

  FUNCTION RECO3(tensor, factors) RESULT(Xtilde)
    !==================================================================================
    ! This function computes the reconstruction of a core tensor that is employed in the Tucker decomposition.
    ! The NN-mode (in this case, 3) of the core tensor is returned. If the core tensor is needed, just employ
    ! the proper TENSOR(N) function.
    ! INPUT
    ! - tensor    : (DTENSOR3) the input tensor
    ! - factors   : (MATRIX_LIST) the list of factor matrices (in this case, Tucker factors)
    ! OUTPUT
    ! - Xtilde    : (REAL*8) the matrix of size (R_NN, prod_{n=1}^NN-1 R_i)
    !==================================================================================   
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
       ! In Xhat2, store the matrix product of the mode_n tensor and the factor matrix
       ALLOCATE(Xtilde(newmodes(ii), PRODUCT(newmodes)/newmodes(ii)))
       Xtilde = MTML(TRANSPOSE(factors(ii)%matr), Xhat)
       DEALLOCATE(Xhat)
       ALLOCATE(Xhat(newmodes(ii+1),PRODUCT(newmodes)/newmodes(ii+1) ))
       Xhat = (TENSOR3(newmodes,Xtilde,ii)).MODE.(ii+1)
       ! Now, Xhat has become the n-mode representation of the running mode-n product.
       ! Xtilde can be deallocated.
       DEALLOCATE(Xtilde)
    END DO
    ! The last Xtilde will contain the mode_NN representation of the tensor
    newmodes(NN)=ranks(NN)
    ALLOCATE(Xtilde(newmodes(NN),PRODUCT(newmodes(1:NN-1))) )
    Xtilde = MTML(TRANSPOSE(factors(NN)%matr), Xhat)
    RETURN
  END FUNCTION RECO3

  
  FUNCTION RECO4(tensor, factors) RESULT(Xtilde)
    !==================================================================================
    ! This function computes the reconstruction of a core tensor that is employed in the Tucker decomposition.
    ! The NN-mode (in this case, 4) of the core tensor is returned. If the core tensor is needed, just employ
    ! the proper TENSOR(N) function.
    ! INPUT
    ! - tensor    : (DTENSOR4) the input tensor
    ! - factors   : (MATRIX_LIST) the list of factor matrices (in this case, Tucker factors)
    ! OUTPUT
    ! - Xtilde    : (REAL*8) the matrix of size (R_NN, prod_{n=1}^NN-1 R_i)
    !==================================================================================   
    TYPE(DTENSOR3) :: tensor
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
       ! In Xhat2, store the matrix product of the mode_n tensor and the factor matrix
       ALLOCATE(Xtilde(newmodes(ii), PRODUCT(newmodes)/newmodes(ii)))
       Xtilde = MTML(TRANSPOSE(factors(ii)%matr), Xhat)
       DEALLOCATE(Xhat)
       ALLOCATE(Xhat(newmodes(ii+1),PRODUCT(newmodes)/newmodes(ii+1) ))
       Xhat = (TENSOR4(newmodes,Xtilde,ii)).MODE.(ii+1)
       ! Now, Xhat has become the n-mode representation of the running mode-n product.
       ! Xtilde can be deallocated.
       DEALLOCATE(Xtilde)
    END DO
    ! The last Xtilde will contain the mode_NN representation of the tensor
    newmodes(NN)=ranks(NN)
    ALLOCATE(Xtilde(newmodes(NN),PRODUCT(newmodes(1:NN-1))) )
    Xtilde = MTML(TRANSPOSE(factors(NN)%matr), Xhat)
    RETURN
  END FUNCTION RECO4


  FUNCTION RECO5(tensor, factors) RESULT(Xtilde)
    !==================================================================================
    ! This function computes the reconstruction of a core tensor that is employed in the Tucker decomposition.
    ! The NN-mode (in this case, 5) of the core tensor is returned. If the core tensor is needed, just employ
    ! the proper TENSOR(N) function.
    ! INPUT
    ! - tensor    : (DTENSOR5) the input tensor
    ! - factors   : (MATRIX_LIST) the list of factor matrices (in this case, Tucker factors)
    ! OUTPUT
    ! - Xtilde    : (REAL*8) the matrix of size (R_NN, prod_{n=1}^NN-1 R_i)
    !==================================================================================   
    TYPE(DTENSOR3) :: tensor
    TYPE(MATRIX_LIST) :: factors(5)
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
       ! In Xhat2, store the matrix product of the mode_n tensor and the factor matrix
       ALLOCATE(Xtilde(newmodes(ii), PRODUCT(newmodes)/newmodes(ii)))
       Xtilde = MTML(TRANSPOSE(factors(ii)%matr), Xhat)
       DEALLOCATE(Xhat)
       ALLOCATE(Xhat(newmodes(ii+1),PRODUCT(newmodes)/newmodes(ii+1) ))
       Xhat = (TENSOR5(newmodes,Xtilde,ii)).MODE.(ii+1)
       ! Now, Xhat has become the n-mode representation of the running mode-n product.
       ! Xtilde can be deallocated.
       DEALLOCATE(Xtilde)
    END DO
    ! The last Xtilde will contain the mode_NN representation of the tensor
    newmodes(NN)=ranks(NN)
    ALLOCATE(Xtilde(newmodes(NN),PRODUCT(newmodes(1:NN-1))) )
    Xtilde = MTML(TRANSPOSE(factors(NN)%matr), Xhat)
    RETURN
  END FUNCTION RECO5



  


  !======================================================= 
  !======================================================= 
  ! HOSVD
  !======================================================= 
  !======================================================= 
    
  
  SUBROUTINE HOSVD3(tens,ranks,core,factors)
    !==================================================================================
    !Returns core and factors of the Tucker Decomposition using HOSVD.
    !INPUT/OUTPUT:
    !- tens 		: (DTENSOR3) the input tensor (size S1,S2,S3)
    !- ranks 		: (INTEGER*4) the vector of ranks (size 3) 
    !- core             : (DTENSOR3) the core tensor (size S1,S2,S3)
    !- factors 		: (MATRIX_LIST) vector (size 3) of factor matrices (size Si,Ri)
    !==================================================================================
    ! INOUT VARIABLES
    TYPE(DTENSOR3) :: tens,core
    INTEGER*4 :: ranks(3)
    TYPE(MATRIX_LIST) :: factors(3)
    ! UTILITY VARIABLES
    INTEGER*4 :: ii, info
    INTEGER*4 :: new(3) ! new modes 
    REAL*8, ALLOCATABLE :: SIG(:,:), UU(:,:), VVT(:,:), res(:,:)
    ! ALLOCATE THE FACTORS
    DO ii=1,3
       ALLOCATE(factors(ii)%matr(tens%modes(ii),ranks(ii)))
    END DO
    ! PERFORM SVD ON ALL POSSIBLE MODES
    DO ii=1,3
       ALLOCATE(UU(tens%modes(ii),tens%modes(ii)))
       ALLOCATE(SIG(tens%modes(ii),PRODUCT(tens%modes)/tens%modes(ii)))
       ALLOCATE(VVT(PRODUCT(tens%modes)/tens%modes(ii),PRODUCT(tens%modes)/tens%modes(ii)))
       CALL SVD(tens%elems.MODE.ii,UU,SIG,VVT,info)
       factors(ii)%matr=UU(:,1:ranks(ii))
       DEALLOCATE(UU)
       DEALLOCATE(SIG)
       DEALLOCATE(VVT)
    END DO
    ! ALLOCATE CORE TENSOR
    core%modes=ranks
    ALLOCATE(core%elems(ranks(1),ranks(2),ranks(3)))
    ! COMPUTE THE CORE TENSOR
    res = RECO(tens,factors) ! RECO returns the nn mode
    core = TENSOR3(ranks,res,3) ! go from nn mode to tensor
  END SUBROUTINE HOSVD3
  

  !======================================================= 
  !======================================================= 
  ! HOOI
  !======================================================= 
  !======================================================= 
  

  SUBROUTINE HOOI3(tensor, ranks, core, factors, error, verbose, numiter, thresh)
    ! INOUT VARIABLES
    TYPE(DTENSOR3) :: tensor, core
    TYPE(MATRIX_LIST) :: factors(3)
    INTEGER*4 :: ranks(3)
    REAL*8 :: error
    REAL*8, OPTIONAL :: thresh
    INTEGER*4, OPTIONAL :: numiter
    LOGICAL, OPTIONAL :: verbose
    ! UTILITY VARIABLES
    REAL*8, ALLOCATABLE :: Xhat(:,:), Xtilde(:,:)       !matrices for mode-n product
    REAL*8, ALLOCATABLE :: UU(:,:), VV(:,:), SIGMA(:)   !matrices for SVD
    INTEGER*4 :: ii, jj, NN=SIZE(factors), cnt, idx(2), INFO
    INTEGER*4 :: newmodes(3)
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
          ! KEEP TRACK OF THE CHANGES IN THE MODES
          newmodes = tensor%modes
          ! PERFORM THE MODE-N PRODUCTS WITHOUT SAVING ANY HIGH ORDER TENSOR
          ! 
          ! From now on:
          ! - Xhat1 will store the tensor in one of its mode_n representations.
          ! - Xhat2 will store the matrix product of Xhat1 with one of the factor matrices.
          !
          ALLOCATE(Xhat(newmodes(idx(1))),PRODUCT(newmodes)/newmodes(idx(1)))
          Xhat=tensor.MODE.idx(1)
          DO jj=1,SIZE(idx)
             !
             ! IMPORTANT: recall that we are reindexing, using idx, in order to skip the ii-th index.
             ! As such, all quantities are indexed by idx(jj), and not by jj!
             !
             ! On Xtilde, the matrix product of the transposed factor matrix and the mode_n tensor is stored.
             ! As such, it has size (R_n, W), where W denotes the current product of the other dimensions.
             ! Recall that, after each mode_n product, the n-th dimension is decreased, from I_n->R_n.
             ! The newmodes variable keeps track of these changes in dimension.
             newmodes(idx(jj)) = ranks(idx(jj))
             ALLOCATE(Xtilde(newmodes(idx(jj)), PRODUCT(newmodes)/newmodes(idx(jj))) )
             Xtilde=MTML(TRANSPOSE(factors(idx(jj))%matr), Xhat)
             ! Now that the matrix product has been computed, Xhat can be deallocated and reallocated. It now serves
             ! a different purpose: it must store the tensor in the next mode, so its shape has to change accordingly.
             ! This basically means that it will become a matrix of shape (I_{n}, W), where W is, as before, the
             ! product of the "running sizes".
             DEALLOCATE(Xhat)
             ALLOCATE(Xhat(newmodes(idx(jj+1)), PRODUCT(newmodes)/newmodes(idx(jj+1))))
             Xhat=(TENSOR3(newmodes,Xtilde,idx(jj))).MODE.idx(jj+1)
             ! Now that Xhat has become, once again, the matrix representation of the "current iteration" of the tensor,
             ! Xtilde can be safely deallocated, in order to prepare it to the next iteration, where it will contain
             ! the matrix product of the next step. 
             DEALLOCATE(Xtilde)
          END DO
          ! GET THE FIRST R_N LEADING LEFT SINGULAR VECTORS OF Y_ii (which has not been actually created!)
          CALL SVD(Xhat,UU,SIGMA,VV,INFO)
          factors(ii)%matr=UU(:,ranks(ii))
          DEALLOCATE(Xhat)
       END DO
       !RECONSTRUCT THE TENSOR AND COMPUTE THE ERROR
       ! Compute the error 
       newerror = SQRT(SUM(tensor**2)-SUM(RECO(tensor,factors)**2))/SIZE(tensor)
       relative_error = ABS((error-newerror)/error)
       error = newerror
       ! If verbose is on, output the current error
       IF (PRESENT(verbose).AND.verbose) THEN
          WRITE(*,*) cnt, error, relative_error
       END IF
    END DO
    ! ALLOCATE CORE TENSOR
    core%modes=ranks
    ALLOCATE(core%elems(ranks(1),ranks(2),ranks(3)))
    ! COMPUTE THE CORE TENSOR
    core%modes = ranks
    core%elems = TENSOR(ranks,RECO(tensor,factors),NN)
  END SUBROUTINE HOOI3



END MODULE TUCKER
