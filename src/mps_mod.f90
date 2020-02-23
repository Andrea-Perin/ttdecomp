MODULE MPS_MOD

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
  INTERFACE MPS
     MODULE PROCEDURE MPS3
     MODULE PROCEDURE MPS4
  END INTERFACE MPS

  
CONTAINS

  ! ====================================
  ! ====================================
  ! MPS REPRESENTATION TO TENSOR
  ! ====================================
  ! ====================================

  FUNCTION MPS_TO_TENSOR3(GG)
    ! =========================================================================
    ! Reconstruct a 3-tensor from its MPS form, provided as the list of cores.
    ! INPUT ARGUMENTS
    ! - GG             : (TENSOR_LIST) the list of core tensors of the MPS form
    ! OUTPUT ARGUMENTS
    ! - MPS_TO_TENSOR3 : (DTENSOR3) reconstruction of the original tensor
    ! =========================================================================  
    ! INOUT VARIABLES
    TYPE(TENSOR_LIST) :: GG(3)
    TYPE(DTENSOR3) :: MPS_TO_TENSOR3
    ! UTILITY VARIABLES
    REAL*8, ALLOCATABLE :: MPS_TO_TENS(:,:)
    INTEGER*4 :: ii, order=SIZE(GG), prev_rank, next_rank
    INTEGER*4 :: og_ranks(SIZE(GG))
    REAL*8, ALLOCATABLE :: current_factor(:,:)
    
    ! CHECK IF THE NUMBER OF FACTORS IS RIGHT
    IF (order.NE.3) THEN
       WRITE(*,*) "WARNING (MPS_TO_TENSOR3): only 3 factors are allowed."
    END IF
    ! STORE THE SHAPE OF THE RESULTING TENSOR
    og_ranks = (/ (GG(ii)%cores%modes(2),ii=1,order,1) /)
    ! ALLOCATE MPS_TO_TENS
    ALLOCATE( MPS_TO_TENS(og_ranks(1), SIZE(GG(1)%cores%elems)/og_ranks(1) ) )
    MPS_TO_TENS = GG(1)%cores.MODE.2
    ! UPDATE, LOOPING OVER FACTORS
    DO ii=2,order
       ! STORE THE RANKS OF THE CURRENT FACTOR
       prev_rank = GG(ii)%cores%modes(1) 
       next_rank = GG(ii)%cores%modes(3)
       ! ALLOCATE THE CURRENT FACTOR
       ALLOCATE( current_factor(prev_rank, SIZE(GG(ii)%cores%elems)/prev_rank) )
       current_factor = GG(ii)%cores.MODE.1
       ! UPDATE MPS_TO_TENS
       MPS_TO_TENS = MTML(MPS_TO_TENS,current_factor)
       DEALLOCATE(current_factor)
       MPS_TO_TENS = RESHAPE(MPS_TO_TENS, (/ SIZE(MPS_TO_TENS)/next_rank, next_rank /) )
    END DO
    ! STORE EVERYTHING IN A TENSOR
    MPS_TO_TENSOR3%modes = og_ranks
    MPS_TO_TENSOR3%elems = RESHAPE( SPREAD(MPS_TO_TENS, 3, 1), og_ranks )
  END FUNCTION MPS_TO_TENSOR3
    

  
  FUNCTION MPS_TO_TENSOR4(GG)
    ! =========================================================================
    ! Reconstruct a 4-tensor from its MPS form, provided as the list of cores.
    ! INPUT ARGUMENTS
    ! - GG             : (TENSOR_LIST) the list of core tensors of the MPS form
    ! OUTPUT ARGUMENTS
    ! - MPS_TO_TENSOR4 : (DTENSOR4) reconstruction of the original tensor
    ! =========================================================================
    ! INOUT VARIABLES
    TYPE(TENSOR_LIST) :: GG(4)
    TYPE(DTENSOR4) :: MPS_TO_TENSOR4
    ! UTILITY VARIABLES
    REAL*8, ALLOCATABLE :: MPS_TO_TENS(:,:)
    INTEGER*4 :: ii, order=SIZE(GG), prev_rank, next_rank
    INTEGER*4 :: og_ranks(SIZE(GG))
    REAL*8, ALLOCATABLE :: current_factor(:,:)
    
    ! CHECK IF THE NUMBER OF FACTORS IS RIGHT
    IF (order.NE.4) THEN
       WRITE(*,*) "WARNING (MPS_TO_TENSOR4): only 4 factors are allowed."
    END IF
    ! STORE THE SHAPE OF THE RESULTING TENSOR
    og_ranks = (/ (GG(ii)%cores%modes(2),ii=1,order,1) /)
    ! ALLOCATE MPS_TO_TENS
    ALLOCATE( MPS_TO_TENS(og_ranks(1), SIZE(GG(1)%cores%elems)/og_ranks(1) ) )
    MPS_TO_TENS = GG(1)%cores.MODE.2
    ! UPDATE, LOOPING OVER FACTORS
    DO ii=2,order
       ! STORE THE RANKS OF THE CURRENT FACTOR
       prev_rank = GG(ii)%cores%modes(1) 
       next_rank = GG(ii)%cores%modes(3)
       ! ALLOCATE THE CURRENT FACTOR
       ALLOCATE( current_factor(prev_rank, SIZE(GG(ii)%cores%elems)/prev_rank) )
       current_factor = GG(ii)%cores.MODE.1
       ! UPDATE MPS_TO_TENS
       MPS_TO_TENS = MTML(MPS_TO_TENS,current_factor)
       DEALLOCATE(current_factor)
       MPS_TO_TENS = RESHAPE(MPS_TO_TENS, (/ SIZE(MPS_TO_TENS)/next_rank, next_rank /) )
    END DO
    ! STORE EVERYTHING IN A TENSOR
    MPS_TO_TENSOR4%modes = og_ranks
    MPS_TO_TENSOR4%elems = RESHAPE( SPREAD(SPREAD(MPS_TO_TENS, 3, 1), 4, 1) , og_ranks )
  END FUNCTION MPS_TO_TENSOR4

  
  ! ====================================
  ! ====================================
  ! 2->3 RESHAPE
  ! ====================================
  ! ====================================

  FUNCTION HYPERSHAPE(matrix, modes)
    ! =========================================================================
    ! Reshape a matrix into a 3-dimensional tensor of shape equal to modes.
    ! INPUT ARGUMENTS
    ! - matrix       : (REAL*8) the input matrix to be reshaped
    ! - modes        : (INTEGER*4) the shape (must be 3-dimensional)
    ! OUTPUT ARGUMENTS
    ! - HYPERSHAPE   : (DTENSOR3) the reshaped matrix into its tensor3 form
    ! =========================================================================
    ! INOUT VARIABLES
    REAL*8, ALLOCATABLE :: matrix(:,:)
    INTEGER*4 :: modes(3)
    TYPE(DTENSOR3) :: HYPERSHAPE

    ! SET THE MODES
    HYPERSHAPE%modes = modes
    ! PROPERLY RESHAPE
    HYPERSHAPE%elems = SPREAD(matrix,2,1)
    HYPERSHAPE%elems = RESHAPE(HYPERSHAPE%elems, modes)
    RETURN
  END FUNCTION HYPERSHAPE
    
    

  ! ====================================
  ! ====================================
  ! MATRIX PRODUCT STATES
  ! ====================================
  ! ====================================

  SUBROUTINE MPS3(tensor, GG, eps)
    ! =========================================================================
    ! Compute the MPS form of the input 3-tensor.
    ! INOUT ARGUMENTS
    ! - tensor      : (DTENSOR3) the input tensor
    ! - GG          : (TENSOR_LIST) the factors of the MPS
    ! - eps         : (REAL*8, OPTIONAL) the precision to be attained.
    !                 Numeber of parameters scale with decreasing eps.
    ! =========================================================================
    ! INOUT VARIABLES
    TYPE(DTENSOR3) :: tensor
    TYPE(TENSOR_LIST) :: GG(3)
    REAL*8, OPTIONAL :: eps
    ! UTILITY VARIABLES
    REAL*8 :: delta
    REAL*8, ALLOCATABLE :: Ctmp(:,:), UU(:,:), SIG(:), VVT(:,:), TMP(:,:)
    INTEGER*4 :: ranks(0:SIZE(tensor%modes)) !this array starts from 0 for clarity reasons
    INTEGER*4 :: kk, dd=SIZE(tensor%modes), newshape(2), coreshape(3), INFO
    
    ! COMPUTE DELTA
    IF (PRESENT(eps).AND.(eps.GT.EPSILON(1D0)) ) THEN
       delta = eps
    ELSE
       delta = 1D-3
    END IF
    delta = (delta/(dd-1)) * SQRT(SUM(tensor%elems**2))
    ! SET BOUNDARY CONDITIONS ON THE RANKS
    ranks(0) = 1
    ranks(dd) = 1
    ! STORE THE INITIAL Ctmp DIRECTLY AS A MATRIX
    ALLOCATE(Ctmp( ranks(0)*tensor%modes(1), SIZE(tensor%elems)/(ranks(0)*tensor%modes(1)) ))
    Ctmp = tensor.MODE.1
    ! LOOP OVER THE DIMENSIONS
    DO kk=1,dd-1
       ! RESHAPING STEP ON Ctmp
       newshape = (/ ranks(kk-1)*tensor%modes(kk), SIZE(Ctmp)/(ranks(kk-1)*tensor%modes(kk)) /)
       Ctmp = RESHAPE(Ctmp, newshape)
       ! COMPUTE TRUNCATED SVD
       IF (ALLOCATED(UU)) DEALLOCATE(UU)
       IF (ALLOCATED(SIG)) DEALLOCATE(SIG)
       IF (ALLOCATED(VVT)) DEALLOCATE(VVT)
       CALL TSVD(Ctmp,UU,SIG,VVT,INFO)
       ! DELTA-TRUNCATE SIGMA, ONLY KEEPING SINGVALS>DELTA
       WHERE (SIG.LT.delta)
          SIG = 0D0
       END WHERE
       ! DEFINE THE RANK, TRUNCATE EVERYTHING
       ranks(kk) = COUNT(SIG/=0D0)
       SIG = SIG(1:ranks(kk))
       UU = UU(:,1:ranks(kk))
       VVT = VVT(1:ranks(kk),:)
       ! COMMUNICATE EVENTUAL ERRORS
       IF (INFO.NE.0) THEN
          WRITE(*,*) "WARNING (MPS3): TSVD failed."
       END IF
       ! STORE THE kk-TH CORE
       coreshape = (/ ranks(kk-1), tensor%modes(kk), ranks(kk) /)
       UU = RESHAPE(UU, (/SIZE(UU,1), ranks(kk)/))
       GG(kk)%cores = HYPERSHAPE(UU, coreshape)
       ! COMPUTE THE NEW Ctmp
       Ctmp = TRANSPOSE(MTDG(TRANSPOSE(VVT),SIG))
    END DO
    ! SAVE THE LAST CORE
    coreshape = (/ ranks(dd-1), tensor%modes(dd), ranks(dd) /)
    GG(kk)%cores = HYPERSHAPE(Ctmp, coreshape)
  END SUBROUTINE MPS3


  
  SUBROUTINE MPS4(tensor, GG, eps)
    ! =========================================================================
    ! Compute the MPS form of the input 4-tensor.
    ! INOUT ARGUMENTS
    ! - tensor      : (DTENSOR4) the input tensor
    ! - GG          : (TENSOR_LIST) the factors of the MPS
    ! - eps         : (REAL*8, OPTIONAL) the precision to be attained.
    !                 Numeber of parameters scale with decreasing eps.
    ! =========================================================================
    ! INOUT VARIABLES
    TYPE(DTENSOR4) :: tensor
    TYPE(TENSOR_LIST) :: GG(4)
    REAL*8, OPTIONAL :: eps
    ! UTILITY VARIABLES
    REAL*8 :: delta
    REAL*8, ALLOCATABLE :: Ctmp(:,:), UU(:,:), SIG(:), VVT(:,:), TMP(:,:)
    INTEGER*4 :: ranks(0:SIZE(tensor%modes)) !this array starts from 0 for clarity reasons
    INTEGER*4 :: kk, dd=SIZE(tensor%modes), newshape(2), coreshape(3), INFO
    
    ! COMPUTE DELTA
    IF (PRESENT(eps).AND.(eps.GT.EPSILON(1D0)) ) THEN
       delta = eps
    ELSE
       delta = 1D-3
    END IF
    delta = (delta/(dd-1)) * SQRT(SUM(tensor%elems**2))
    ! SET BOUNDARY CONDITIONS ON THE RANKS
    ranks(0) = 1
    ranks(dd) = 1
    ! STORE THE INITIAL Ctmp DIRECTLY AS A MATRIX
    ALLOCATE(Ctmp( ranks(0)*tensor%modes(1), SIZE(tensor%elems)/(ranks(0)*tensor%modes(1)) ))
    Ctmp = tensor.MODE.1
    ! LOOP OVER THE DIMENSIONS
    DO kk=1,dd-1
       ! RESHAPING STEP ON Ctmp
       newshape = (/ ranks(kk-1)*tensor%modes(kk), SIZE(Ctmp)/(ranks(kk-1)*tensor%modes(kk)) /)
       Ctmp = RESHAPE(Ctmp, newshape)
       ! COMPUTE TRUNCATED SVD
       IF (ALLOCATED(UU)) DEALLOCATE(UU)
       IF (ALLOCATED(SIG)) DEALLOCATE(SIG)
       IF (ALLOCATED(VVT)) DEALLOCATE(VVT)
       CALL TSVD(Ctmp,UU,SIG,VVT,INFO)
       ! DELTA-TRUNCATE SIGMA, ONLY KEEPING SINGVALS>DELTA
       WHERE (SIG.LT.delta)
          SIG = 0D0
       END WHERE
       ! DEFINE THE RANK, TRUNCATE EVERYTHING
       ranks(kk) = COUNT(SIG/=0D0)
       SIG = SIG(1:ranks(kk))
       UU = UU(:,1:ranks(kk))
       VVT = VVT(1:ranks(kk),:)
       ! COMMUNICATE EVENTUAL ERRORS
       IF (INFO.NE.0) THEN
          WRITE(*,*) "WARNING (MPS4): TSVD failed."
       END IF
       ! STORE THE kk-TH CORE
       coreshape = (/ ranks(kk-1), tensor%modes(kk), ranks(kk) /)
       UU = RESHAPE(UU, (/SIZE(UU,1), ranks(kk)/))
       GG(kk)%cores = HYPERSHAPE(UU, coreshape)
       ! COMPUTE THE NEW Ctmp
       Ctmp = TRANSPOSE(MTDG(TRANSPOSE(VVT),SIG))
    END DO
    ! SAVE THE LAST CORE
    coreshape = (/ ranks(dd-1), tensor%modes(dd), ranks(dd) /)
    GG(kk)%cores = HYPERSHAPE(Ctmp, coreshape)
  END SUBROUTINE MPS4

  

  
  
END MODULE MPS_MOD
