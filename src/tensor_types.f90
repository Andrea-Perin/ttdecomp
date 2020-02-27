MODULE TENSOR_TYPES
  IMPLICIT NONE


  ! LIST OF MATRICES (ARBITRARY DIMENSION)
  TYPE :: matrix_list
     REAL*8, ALLOCATABLE :: matr(:,:)
  END TYPE  matrix_list
  
  
  ! TENSOR TYPES (DOUBLE PRECISION)
  TYPE DTENSOR1
     INTEGER*4::modes(1)
     REAL*8,ALLOCATABLE::elems(:)
  END TYPE DTENSOR1

  TYPE DTENSOR2
     INTEGER*4::modes(2)
     REAL*8,ALLOCATABLE::elems(:,:)
  END TYPE DTENSOR2

  TYPE DTENSOR3
     INTEGER*4::modes(3)
     REAL*8,ALLOCATABLE::elems(:,:,:)
  END TYPE DTENSOR3

  TYPE DTENSOR4
     INTEGER*4::modes(4)
     REAL*8,ALLOCATABLE::elems(:,:,:,:)
  END TYPE DTENSOR4

  TYPE DTENSOR5
     INTEGER*4::modes(5)
     REAL*8,ALLOCATABLE::elems(:,:,:,:,:)
  END TYPE DTENSOR5

  TYPE DTENSOR6
     INTEGER*4::modes(6)
     REAL*8,ALLOCATABLE::elems(:,:,:,:,:,:)
  END TYPE DTENSOR6

  TYPE DTENSOR7
     INTEGER*4::modes(7)
     REAL*8,ALLOCATABLE::elems(:,:,:,:,:,:,:)
  END TYPE DTENSOR7


  ! LIST OF TENSORS
  TYPE :: tensor_list
     TYPE(DTENSOR3) :: cores
  END TYPE  tensor_list


  ! INTERFACE FOR ALL-ORTHOGONALITY CHECK
  INTERFACE IS_ALL_ORTHOGONAL
     MODULE PROCEDURE IS_ALL_ORTHOGONAL3
     MODULE PROCEDURE IS_ALL_ORTHOGONAL4
  END INTERFACE IS_ALL_ORTHOGONAL

  
  
CONTAINS

  !-----------------------TO TENSOR2 SUBROUTINE
  FUNCTION TENSOR2(modes,mat,cm)
    ! cm = current mode
    INTEGER::modes(2),cm
    REAL*8,ALLOCATABLE::mat(:,:)
    TYPE(DTENSOR2)::tensor2
    ALLOCATE(tensor2%elems(modes(1),modes(2)))
    tensor2%modes=modes
    IF (cm.EQ.1) THEN
       tensor2%elems=mat
    ELSE IF (cm.EQ.2) THEN 
       tensor2%elems=TRANSPOSE(mat)
    ELSE
       PRINT*, "ERROR"
    ENDIF
    RETURN
  END FUNCTION TENSOR2

  !-----------------------TO TENSOR3 SUBSROUTINE
  FUNCTION TENSOR3(modes,mat,cm)
    ! cm = current mode
    INTEGER::modes(3),cm,new_modes(2)
    REAL*8,ALLOCATABLE::mat(:,:)
    TYPE(DTENSOR3)::tensor3
    INTEGER::ii
    ALLOCATE(tensor3%elems(modes(1),modes(2),modes(3)))
    tensor3%modes=modes
    IF (cm.EQ.1) THEN
       new_modes=modes(2:3) 
       DO ii=1,modes(1)
          tensor3%elems(ii,:,:)=RESHAPE(mat(ii,:),new_modes)
       END DO
    ELSE IF (cm.EQ.2) THEN
       new_modes=(/modes(1),modes(3)/) 
       DO ii=1,modes(2)
          tensor3%elems(:,ii,:)=RESHAPE(mat(ii,:),new_modes)
       END DO
    ELSE IF (cm.EQ.3) THEN
       new_modes=modes(1:2) 
       DO ii=1,modes(3)
          tensor3%elems(:,:,ii)=RESHAPE(mat(ii,:),new_modes)
       END DO
    ELSE
       PRINT*, "ERROR"
    ENDIF
    RETURN
  END FUNCTION TENSOR3

  !-----------------------TO TENSOR4 SUBSROUTINE
  FUNCTION TENSOR4(modes,mat,cm)
    ! cm = current mode
    INTEGER::modes(4),cm,new_modes(3)
    REAL*8,ALLOCATABLE::mat(:,:)
    TYPE(DTENSOR4)::tensor4
    INTEGER::ii
    ALLOCATE(tensor4%elems(modes(1),modes(2),modes(3),modes(4)))
    tensor4%modes=modes

    IF (cm.EQ.1) THEN
       new_modes=modes(2:4) 
       DO ii=1,modes(1)
          tensor4%elems(ii,:,:,:)=RESHAPE(mat(ii,:),new_modes)
       END DO
    ELSE IF (cm.EQ.2) THEN
       new_modes=(/modes(1),modes(3),modes(4)/) 
       DO ii=1,modes(2)
          tensor4%elems(:,ii,:,:)=RESHAPE(mat(ii,:),new_modes)
       END DO
    ELSE IF (cm.EQ.3) THEN
       new_modes=(/modes(1),modes(2),modes(4)/) 
       DO ii=1,modes(3)
          tensor4%elems(:,:,ii,:)=RESHAPE(mat(ii,:),new_modes)
       END DO
    ELSE IF (cm.EQ.4) THEN
       new_modes=modes(1:3) 
       DO ii=1,modes(4)
          tensor4%elems(:,:,:,ii)=RESHAPE(mat(ii,:),new_modes)
       END DO
    ELSE
       PRINT*, "ERROR"
    ENDIF
    RETURN
  END FUNCTION TENSOR4








! ==================================================
! ==================================================
! CHECK PROPERTIES OF TENSORS
! ==================================================
! ==================================================

  FUNCTION IS_ALL_ORTHOGONAL3(tensor, eps)
    ! ===============================================================
    ! Checks if the DTENSOR3 given as input is all-orthogonal,
    ! following the definition given by Cichocki.
    ! INPUT ARGUMENTS
    ! - tensor             : (DTENSOR3) the tensor to check
    ! - eps                : (REAL*8, OPTIONAL) the tolerance for the "=0" checks
    ! OUTPUT ARGUMENTS
    ! - IS_ALL_ORTHOGONAL3 : (LOGICAL) the result of the check
    ! ===============================================================  
    ! INOUT VARIABLES
    TYPE(DTENSOR3) :: tensor
    REAL*8, OPTIONAL :: eps
    LOGICAL :: IS_ALL_ORTHOGONAL3
    ! UTILITY VARIABLES
    INTEGER*4 :: ii,jj,kk,blklen
    REAL*8 :: tol=1D-5

    ! FUNCTION DEFINITION
    !----------------------------------------
    ! All-orthogonality has 2 checks:
    ! 1. all slices are mutually ortohogonal
    ! 2. frobenius norms of slices in each mode are decreasing for increasing index
    !----------------------------------------
    IF (PRESENT(eps).AND.(eps.GT.EPSILON(1D0))) THEN
       tol=eps
    END IF
    IS_ALL_ORTHOGONAL3=.TRUE.
    ! slice 1: mutually orthogonal?
    DO ii=1,tensor%modes(1)
       DO jj=ii+1,tensor%modes(1)
          IF (ABS(SUM(tensor%elems(ii,:,:)*tensor%elems(jj,:,:))).GT.tol*SIZE(tensor%elems(ii,:,:))) THEN
             IS_ALL_ORTHOGONAL3=.FALSE.
             !print*, "slice 1 not orth"
          END IF
       END DO
    END DO
    ! slice 2: mutually orthogonal?
    DO ii=1,tensor%modes(2)
       DO jj=ii+1,tensor%modes(2)
          IF (ABS(SUM(tensor%elems(:,ii,:)*tensor%elems(:,jj,:))).GT.tol*SIZE(tensor%elems(:,ii,:))) THEN
             IS_ALL_ORTHOGONAL3=.FALSE.
             !print*, "slice 2 not orth"
          END IF
       END DO
    END DO
    ! slice 3: mutually orthogonal?
    DO ii=1,tensor%modes(3)
       DO jj=ii+1,tensor%modes(3)
          IF (ABS(SUM(tensor%elems(:,:,ii)*tensor%elems(:,:,jj))).GT.tol*SIZE(tensor%elems(:,:,ii))) THEN
             IS_ALL_ORTHOGONAL3=.FALSE.
             !print*, "slice 3 not orth"
          END IF
       END DO
    END DO
    ! CHECK THE ORDERING OF THE SLICES
    ! slice 1: ordered?
    DO ii=1,tensor%modes(1)
       IF (SUM(tensor%elems(ii,:,:)**2).LT.SUM(tensor%elems(ii+1,:,:)**2)) THEN
          IS_ALL_ORTHOGONAL3=.FALSE.
             !print*, "slice 1 not ordered"
       END IF
    END DO
    ! slice 2: ordered?
    DO ii=1,tensor%modes(2)
       IF (SUM(tensor%elems(:,ii,:)**2).LT.SUM(tensor%elems(:,ii+1,:)**2)) THEN
          IS_ALL_ORTHOGONAL3=.FALSE.
             !print*, "slice 2 not ordered"
       END IF
    END DO
    ! slice 3: ordered?
    DO ii=1,tensor%modes(3)
       IF (SUM(tensor%elems(:,:,ii)**2).LT.SUM(tensor%elems(:,:,ii+1)**2)) THEN
          IS_ALL_ORTHOGONAL3=.FALSE.
             !print*, "slice 3 not ordered"
       END IF
    END DO
    ! RETURN THE RESULT OF THE CHECK
    RETURN
  END FUNCTION IS_ALL_ORTHOGONAL3



  FUNCTION IS_ALL_ORTHOGONAL4(tensor, eps)
    ! ===============================================================
    ! Checks if the DTENSOR4 given as input is all-orthogonal,
    ! following the definition given by Cichocki.
    ! INPUT ARGUMENTS
    ! - tensor             : (DTENSOR4) the tensor to check
    ! - eps                : (REAL*8, OPTIONAL) the tolerance for the "=0" checks
    ! OUTPUT ARGUMENTS
    ! - IS_ALL_ORTHOGONAL3 : (LOGICAL) the result of the check
    ! ===============================================================  
    ! INOUT VARIABLES
    TYPE(DTENSOR4) :: tensor
    REAL*8, OPTIONAL :: eps
    LOGICAL :: IS_ALL_ORTHOGONAL4
    ! UTILITY VARIABLES
    INTEGER*4 :: ii,jj,kk,blklen
    REAL*8 :: tol=1D-5

    ! FUNCTION DEFINITION
    !----------------------------------------
    ! All-orthogonality has 2 checks:
    ! 1. all slices are mutually ortohogonal
    ! 2. frobenius norms of slices in each mode are decreasing for increasing index
    !----------------------------------------
    IF (PRESENT(eps).AND.(eps.GT.EPSILON(1D0))) THEN
       tol=eps
    END IF
    IS_ALL_ORTHOGONAL4=.TRUE.
    ! slice 1: mutually orthogonal?
    DO ii=1,tensor%modes(1)
       DO jj=ii+1,tensor%modes(1)
          IF (ABS(SUM(tensor%elems(ii,:,:,:)*tensor%elems(jj,:,:,:))).GT.tol*SIZE(tensor%elems(ii,:,:,:))) THEN
             IS_ALL_ORTHOGONAL4=.FALSE.
          END IF
       END DO
    END DO
    ! slice 2: mutually orthogonal?
    DO ii=1,tensor%modes(2)
       DO jj=ii+1,tensor%modes(2)
          IF (ABS(SUM(tensor%elems(:,ii,:,:)*tensor%elems(:,jj,:,:))).GT.tol*SIZE(tensor%elems(:,ii,:,:))) THEN
             IS_ALL_ORTHOGONAL4=.FALSE.
          END IF
       END DO
    END DO
    ! slice 3: mutually orthogonal?
    DO ii=1,tensor%modes(3)
       DO jj=ii+1,tensor%modes(3)
          IF (ABS(SUM(tensor%elems(:,:,ii,:)*tensor%elems(:,:,jj,:))).GT.tol*SIZE(tensor%elems(:,:,ii,:))) THEN
             IS_ALL_ORTHOGONAL4=.FALSE.
          END IF
       END DO
    END DO
    ! slice 4: mutually orthogonal?
    DO ii=1,tensor%modes(4)
       DO jj=ii+1,tensor%modes(4)
          IF (ABS(SUM(tensor%elems(:,:,:,ii)*tensor%elems(:,:,:,jj))).GT.tol*SIZE(tensor%elems(:,:,:,ii))) THEN
             IS_ALL_ORTHOGONAL4=.FALSE.
          END IF
       END DO
    END DO
    ! CHECK THE ORDERING OF THE SLICES
    ! slice 1: ordered?
    DO ii=1,tensor%modes(1)
       IF (SUM(tensor%elems(ii,:,:,:)**2).LT.SUM(tensor%elems(ii+1,:,:,:)**2)) THEN
          IS_ALL_ORTHOGONAL4=.FALSE.
       END IF
    END DO
    ! slice 2: ordered?
    DO ii=1,tensor%modes(2)
       IF (SUM(tensor%elems(:,ii,:,:)**2).LT.SUM(tensor%elems(:,ii+1,:,:)**2)) THEN
          IS_ALL_ORTHOGONAL4=.FALSE.
       END IF
    END DO
    ! slice 3: ordered?
    DO ii=1,tensor%modes(3)
       IF (SUM(tensor%elems(:,:,ii,:)**2).LT.SUM(tensor%elems(:,:,ii+1,:)**2)) THEN
          IS_ALL_ORTHOGONAL4=.FALSE.
       END IF
    END DO
    ! slice 4: ordered?
    DO ii=1,tensor%modes(4)
       IF (SUM(tensor%elems(:,:,:,ii)**2).LT.SUM(tensor%elems(:,:,:,ii+1)**2)) THEN
          IS_ALL_ORTHOGONAL4=.FALSE.
       END IF
    END DO
    ! RETURN THE RESULT OF THE CHECK
    RETURN
  END FUNCTION IS_ALL_ORTHOGONAL4


  
  
END MODULE TENSOR_TYPES
