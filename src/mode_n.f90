MODULE MODE_N
  USE TENSOR_TYPES
  USE MAT_UTILS
  IMPLICIT NONE

  ! INTERFACE FOR ALL TENSOR TYPES
  INTERFACE OPERATOR (.MODE.)
     MODULE PROCEDURE DTENSOR1_MODE_N
     MODULE PROCEDURE DTENSOR2_MODE_N
     MODULE PROCEDURE DTENSOR3_MODE_N 
     MODULE PROCEDURE DTENSOR4_MODE_N 
  END INTERFACE OPERATOR (.MODE.)

  !INTERFACE NPROD
  !   MODULE PROCEDURE NPROD3
  !END INTERFACE NPROD

CONTAINS

  
  FUNCTION NPROD3(tens,matt,nn)
    !====================================================================
    !Computes the n mode product between a tensor and a matrix
    !INPUT:
    !- tens             : (DTENSOR3) the input tensor
    !- mat	        : (REAL*8) the input matrix
    !- nn               : (INTEGER*4) the mode of the product
    !OUTPUT:
    !- NPROD3   	: (DTENSOR3) the output tensor
    !====================================================================
    ! INPUT ARGUMENTS
    TYPE(DTENSOR3) :: tens
    REAL*8, ALLOCATABLE :: matt(:,:)
    REAL*8 ,ALLOCATABLE :: tmp_mat(:,:)
    INTEGER*4 :: nn
    ! OUTPUT
    TYPE(DTENSOR3) :: NPROD3
    ! UTILITY VARIABLES
    INTEGER*4 :: ii
    ! FUNCTION DEFINITION
    tmp_mat=MTML(matt,tens.MODE.nn)
    CALL TENSOR3(tens%modes,tmp_mat,NPROD3,nn)
    RETURN
  END FUNCTION NPROD3


  
  FUNCTION REIDX(modes,nn,idxs)
    !====================================================================
    !Returns the multi index.
    !INPUT:
    !- modes			: (INTEGER*4) the modes of the input tensor
    !- nn			: (INTEGER*4) the mode for which to compute the multi index
    !- idxs			: (INTEGER*4) the values of the indexes for each dimension
    !OUTPUT:
    !- REIDX			: (INTEGER*4) the multi index
    !====================================================================
    ! INPUT ARGUMENTS
    INTEGER*4, INTENT(IN) :: modes(:), idxs(:), nn
    ! OUTPUT
    INTEGER*4 :: REIDX
    ! UTILITY VARIABLES
    INTEGER*4 :: ii
    ! FUNCTION DEFINITION		
    REIDX = 1
    DO ii=1,nn-1
       REIDX = REIDX + (idxs(ii)-1)*PRODUCT(modes(:ii-1))
    END DO
    DO ii=nn+1,SIZE(modes)
       REIDX = REIDX + (idxs(ii)-1)*PRODUCT(modes(:nn-1))*PRODUCT(modes(nn+1:ii-1))
    END DO
    RETURN
  END FUNCTION REIDX



  FUNCTION DTENSOR1_MODE_N(tensor,mm)
    !====================================================================
    !Returns the mm-mode representation of tensor.
    !INPUT:
    !- tensor			: (DTENSOR1) the input tensor
    !- mm			: (INTEGER*4) the mode
    !OUTPUT:
    !- DTENSOR1_MODE_N	: (REAL*8) the output matrix
    !====================================================================
    ! INPUT ARGUMENTS
    TYPE(DTENSOR1), INTENT(IN) :: tensor
    INTEGER*4, INTENT(IN) :: mm
    ! OUTPUT
    REAL*8 :: DTENSOR1_MODE_N(SIZE(tensor%elems))
    ! FUNCTION DEFINITION		
    IF (mm.NE.1) THEN
       WRITE(*,*) "ERROR (DTENSOR1_MODE_N): mode must be 1. Found ",mm
    END IF
    DTENSOR1_MODE_N=tensor%elems
    RETURN
  END FUNCTION DTENSOR1_MODE_N



  FUNCTION DTENSOR2_MODE_N(tensor,mm)
    !====================================================================
    !Returns the mm-mode representation of tensor.
    !INPUT:
    !- tensor			: (DTENSOR2) the input tensor
    !- mm			: (INTEGER*4) the mode
    !OUTPUT:
    !- DTENSOR2_MODE_N	: (REAL*8) the output matrix
    !====================================================================
    ! INPUT ARGUMENTS
    TYPE(DTENSOR2), INTENT(IN) :: tensor
    INTEGER*4, INTENT(IN) :: mm
    ! OUTPUT
    REAL*8, ALLOCATABLE :: DTENSOR2_MODE_N(:,:)
    ! UTILITY VARIABLES
    INTEGER*4 :: ii
    ! FUNCTION DEFINITION		
    IF ((mm.LT.1).OR.(mm.GT.2)) THEN
       WRITE(*,*) "ERROR (DTENSOR2_MODE_N): mode must be in [1,2]. Found ",mm
    ELSE
       ALLOCATE(DTENSOR2_MODE_N(SIZE(tensor%elems,1),SIZE(tensor%elems,2)))
       IF (mm.EQ.1) THEN
          DTENSOR2_MODE_N=tensor%elems
       ELSE
          DTENSOR2_MODE_N=TRANSPOSE(tensor%elems)
       END IF
    END IF
    RETURN
  END FUNCTION DTENSOR2_MODE_N



  FUNCTION DTENSOR3_MODE_N(tensor,mm)
    !====================================================================
    !Returns the mm-mode representation of tensor.
    !INPUT:
    !- tensor			: (DTENSOR3) the input tensor
    !- mm			: (INTEGER*4) the mode
    !OUTPUT:
    !- DTENSOR2_MODE_N	: (REAL*8) the output matrix
    !====================================================================
    ! INPUT ARGUMENTS
    TYPE(DTENSOR3), INTENT(IN) :: tensor
    INTEGER*4, INTENT(IN) :: mm
    ! OUTPUT
    REAL*8, ALLOCATABLE :: DTENSOR3_MODE_N(:,:)
    ! UTILITY VARIABLES
    INTEGER*4 :: nrows, ncols, i1, i2, i3
    ! FUNCTION DEFINITION		
    IF ((mm.LT.1).OR.(mm.GT.3)) THEN
       WRITE(*,*) "ERROR (DTENSOR3_MODE_N): mode must be in [1,2,3]. Found ",mm
    ELSE
       nrows=tensor%modes(mm)
       ncols=PRODUCT(tensor%modes)/nrows
       ALLOCATE(DTENSOR3_MODE_N(nrows,ncols))
       SELECT CASE(mm)
       CASE(1)
          DO i3=1,tensor%modes(3)
             DO i2=1,tensor%modes(2)
                DTENSOR3_MODE_N(:,REIDX(tensor%modes,1,(/i1,i2,i3/))) = tensor%elems(:,i2,i3)
             END DO
          END DO
       CASE(2)
          DO i3=1,tensor%modes(3)
             DO i1=1,tensor%modes(1)						
                DTENSOR3_MODE_N(:,REIDX(tensor%modes,2,(/i1,i2,i3/))) = tensor%elems(i1,:,i3)
             END DO
          END DO
       CASE(3)
          DO i2=1,tensor%modes(2)
             DO i1=1,tensor%modes(1)
                DTENSOR3_MODE_N(:,REIDX(tensor%modes,3,(/i1,i2,i3/))) = tensor%elems(i1,i2,:)
             END DO
          END DO
       END SELECT
    END IF
    RETURN
  END FUNCTION DTENSOR3_MODE_N



  FUNCTION DTENSOR4_MODE_N(tensor,mm)
    !====================================================================
    !Returns the mm-mode representation of tensor.
    !INPUT:
    !- tensor			: (DTENSOR4) the input tensor
    !- mm			: (INTEGER*4) the mode
    !OUTPUT:
    !- DTENSOR4_MODE_N	: (REAL*8) the output matrix
    !====================================================================
    ! INPUT ARGUMENTS
    TYPE(DTENSOR4), INTENT(IN) :: tensor
    INTEGER*4, INTENT(IN) :: mm
    ! OUTPUT
    REAL*8, ALLOCATABLE :: DTENSOR4_MODE_N(:,:)
    ! UTILITY VARIABLES
    INTEGER*4 :: nrows, ncols, i1, i2, i3, i4
    ! FUNCTION DEFINITION		
    IF ((mm.LT.1).OR.(mm.GT.4)) THEN
       WRITE(*,*) "ERROR (DTENSOR4_MODE_N): mode must be in [1,2,3,4]. Found ",mm
    ELSE
       nrows=tensor%modes(mm)
       ncols=PRODUCT(tensor%modes)/nrows
       ALLOCATE(DTENSOR4_MODE_N(nrows,ncols))
       SELECT CASE(mm)
       CASE(1)
          DO i4=1,tensor%modes(4)
             DO i3=1,tensor%modes(3)
                DO i2=1,tensor%modes(2)
                   DTENSOR4_MODE_N(:,REIDX(tensor%modes,1,(/i1,i2,i3,i4/))) = tensor%elems(:,i2,i3,i4)
                END DO
             END DO
          END DO
       CASE(2)
          DO i4=1,tensor%modes(4)
             DO i3=1,tensor%modes(3)
                DO i1=1,tensor%modes(1)
                   DTENSOR4_MODE_N(:,REIDX(tensor%modes,2,(/i1,i2,i3,i4/))) = tensor%elems(i1,:,i3,i4)
                END DO
             END DO
          END DO
       CASE(3)
          DO i4=1,tensor%modes(4)
             DO i2=1,tensor%modes(2)
                DO i1=1,tensor%modes(1)
                   DTENSOR4_MODE_N(:,REIDX(tensor%modes,3,(/i1,i2,i3,i4/))) = tensor%elems(i1,i2,:,i4)
                END DO
             END DO
          END DO
       CASE(4)
          DO i3=1,tensor%modes(3)
             DO i2=1,tensor%modes(2)
                DO i1=1,tensor%modes(1)
                   DTENSOR4_MODE_N(:,REIDX(tensor%modes,4,(/i1,i2,i3,i4/))) = tensor%elems(i1,i2,i3,:)
                END DO
             END DO
          END DO
       END SELECT
    END IF
    RETURN
  END FUNCTION DTENSOR4_MODE_N



END MODULE MODE_N


