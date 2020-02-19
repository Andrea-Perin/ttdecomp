MODULE MOD_TUCKER
  
  ! IN THIS MODULE:
  !- HOSVD

  IMPLICIT NONE

  INTERFACE HOSVD
     MODULE PROCEDURE HOSVD3
  END INTERFACE HOSVD


CONTAINS

  
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
    CALL TENSOR3(ranks,res,core,3) ! go from nn mode to tensor
  END SUBROUTINE HOSVD
