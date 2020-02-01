!------------------------MAIN PROGRAM
PROGRAM MAIN
  USE PRODUCTS
  USE TENSOR_TYPES
  IMPLICIT NONE
  REAL*8,ALLOCATABLE::AA(:,:),BB(:,:),CC(:,:),DD(:,:),IN(:,:),DOM(:,:)
  REAL*8,ALLOCATABLE::V1(:),V2(:),V3(:)
  REAL*8::sigma
  INTEGER,ALLOCATABLE::dim(:)
  INTEGER,ALLOCATABLE::test(:,:)
  INTEGER::UU,nn,rows,cols,kk,ii
  TYPE(DTENSOR3)::my_tens
  OPEN(332,file='./mnist_1k.csv',status="old",action="read")
  READ (332, *) nn
  ALLOCATE(dim(nn))
  READ(332,*) dim
  rows=dim(1)
  cols=1
  DO kk=1,nn-1
     cols=cols*dim(1+kk)
  END DO
  ALLOCATE(IN(cols,rows))
  READ(332,*) IN
  IN=TRANSPOSE(IN)
  CLOSE(332)
  OPEN(unit=333,file='output.dat',status="unknown",action="write")
  ALLOCATE(AA(2,2))
  ALLOCATE(BB(2,2))
  ALLOCATE(V1(2))
  ALLOCATE(V2(2))
  AA = TRANSPOSE(RESHAPE((/ 1.d0, 2.d0, 3.d0, 4.d0 /),(/SIZE(AA,2),SIZE(AA,1)/)))
  BB = TRANSPOSE(RESHAPE((/ 5.d0, 6.d0, 7.d0, 8.d0 /),(/SIZE(BB,2),SIZE(BB,1)/)))
  V1=(/1.d0, 2.d0/)
  V2=(/3.d0, 4.d0/)
  CALL KRONECKER(BB,AA,CC)
  CALL KRONECKER(V1,V2,V3)
  CALL KHATRI_RAO(AA,BB,DD)
  OPEN(unit=333,file='output.dat',status="unknown",action="write")
  UU=333
  ALLOCATE(DOM(1000,1000))
  sigma=10.d0
  !CALL RAND_INIT(DOM,sigma)
  !CALL SHOW(UU,DOM)
  !WRITE (UU,*) "MATRIX:"
  !CALL SHOW(UU,CC)
  !WRITE (UU,*) "VECTOR:"
  !CALL SHOW(UU,V3)
  !WRITE (UU,*) "KHATRI-RAO:"
  !CALL SHOW(UU,DD)
  CLOSE(333)
  !WRITE (6,*) nn
  !WRITE (6,*) "info:",dim
  !WRITE (6,*) "MATRIX:"
  !CALL SHOW(6,IN)
!   CALL TENSOR3(dim,IN,my_tens)
!   DO ii=1,dim(2)
!      !WRITE(6,*) INT(my_tens%elems(1,ii,:))
!      WRITE(6,1) my_tens%elems(1,ii,:)
! 1      FORMAT (*(F4.0:x))
  !   END DO
END PROGRAM MAIN
