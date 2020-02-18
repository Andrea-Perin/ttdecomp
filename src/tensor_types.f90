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

  !-----------------------TO TENSOR2 SUBSROUTINE
  SUBROUTINE TENSOR2(modes,mat,tens)
    INTEGER::modes(2)
    REAL*8,ALLOCATABLE::mat(:,:)
    TYPE(DTENSOR2)::tens
    ALLOCATE(tens%elems(modes(1),modes(2)))
    tens%modes=modes
    tens%elems=mat
  END SUBROUTINE TENSOR2

  !-----------------------TO TENSOR3 SUBSROUTINE
  SUBROUTINE TENSOR3(modes,mat,tens)
    INTEGER::modes(3)
    REAL*8,ALLOCATABLE::mat(:,:)
    TYPE(DTENSOR3)::tens
    INTEGER::ii
    ALLOCATE(tens%elems(modes(1),modes(2),modes(3)))
    tens%modes=modes
    DO ii=1,modes(1)
       tens%elems(ii,:,:)=RESHAPE(mat(ii,:),modes(2:3))
    END DO
  END SUBROUTINE TENSOR3

  !-----------------------TO TENSOR4 SUBSROUTINE
  SUBROUTINE TENSOR4(modes,mat,tens)
    INTEGER::modes(4)
    REAL*8,ALLOCATABLE::mat(:,:)
    TYPE(DTENSOR4)::tens
    INTEGER::ii
    ALLOCATE(tens%elems(modes(1),modes(2),modes(3),modes(4)))
    tens%modes=modes
    DO ii=1,modes(1)
       tens%elems(ii,:,:,:)=RESHAPE(mat(ii,:),modes(2:4))
    END DO
  END SUBROUTINE TENSOR4

  
END MODULE TENSOR_TYPES
