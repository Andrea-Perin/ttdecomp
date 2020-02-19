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

  
END MODULE TENSOR_TYPES
