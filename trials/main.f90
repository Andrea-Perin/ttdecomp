PROGRAM MAIN

  USE DEBUGGER
  USE TENSOR_TYPES
  USE MODE_N
  USE SWAP_AXIS
  USE MAT_UTILS
  USE CPD_UTILS
  USE TUCKER
  
  IMPLICIT NONE

  INTEGER,ALLOCATABLE::dimm(:)
  REAL*8,ALLOCATABLE::IN(:,:),boh1(:,:),boh2(:,:)
  INTEGER,ALLOCATABLE::vec(:),ranks(:)
  INTEGER::nn,rows,cols,kk,ii
  TYPE(DTENSOR3)::my_tens,copy,approx,core
  !JUST TO TRY...
  TYPE(matrix_list) :: lista(3)
  INTEGER*4 :: rango
  REAL*8 :: threshold, error, array(3,3), cose(3)
  REAL*8, ALLOCATABLE :: lambdas(:)
  !OPEN(332,file='../data/mnist_1k.csv',status="old",action="read")
  !OPEN(332,file='../data/mnist_small.csv',status="old",action="read")
  !OPEN(332,file='../data/prova.csv',status="old",action="read")
  OPEN(332,file='../data/small_landscape.csv',status="old",action="read")
  READ (332, *) nn
  ALLOCATE(dimm(nn))
  READ(332,*) dimm
  rows=dimm(1)
  cols=1
  DO kk=2,nn
     cols=cols*dimm(kk)
  END DO
  ALLOCATE(IN(cols,rows))
  READ(332,*) IN
  IN=TRANSPOSE(IN)
  CLOSE(332)
  !WRITE (6,*) nn
  !WRITE (6,*) "info:",dimm
  !WRITE (6,*) "MATRIX:"
  my_tens=TENSOR3(dimm,IN,1) ! assuming mode 1
  !DO ii=1,dimm(2)
  !   WRITE(6,1) my_tens%elems(1,ii,:)
!1    FORMAT (*(F4.0:x)) 
  !END DO
  ! WRITE (6,*) "dimensioni del tensore"
  ! ALLOCATE(vec(3))
  ! vec = SHAPE(my_tens%elems)
  ! WRITE (6,*) vec
  ! ALLOCATE(copy%elems(vec(1),vec(2),vec(3)))
  ! copy%elems=my_tens%elems
  ! DO ii=1,dimm(2)
  !    WRITE(6,1) copy%elems(1,ii,:) 
  ! END DO
  ! ALLOCATE(boh1(2,2))
  ! ALLOCATE(boh2(2,2))
  ! boh1 = RESHAPE((/ 1, 2, 3, 4 /), SHAPE(boh1))
  ! PRINT*, "CIAO1"
  ! boh2 = NPROD3(my_tens,boh1,1).MODE.3
  ! PRINT*, "CIAO2"
  ! DO ii=1,SIZE(boh2,1)
  !    WRITE(6,*) boh2(ii,:) 
  ! END DO
  !rango=10
  !threshold=100d0
  !array = transpose(reshape((/ 1D0, 2D0, 3D0, 4D0, 5D0, 6D0, 7D0, 8D0, 9D0 /), shape(array)))
  !CALL CPD3(my_tens, rango, lista, lambdas, error, verbose=.TRUE.)
  !print*, error
  ALLOCATE(ranks(3))
  ranks = (/ 20,20,3 /)
  CALL HOSVD(my_tens,ranks,core,lista)
  !CALL HOOI3(my_tens,ranks,core,lista,error)
  approx = RECO(core,lista)
  !PRINT*, SHAPE(approx%elems)
  OPEN(333,file='../data/hosvd_small_20_20_3.csv',status="unknown",action="write")
  DO ii=1,100
     WRITE(333,*) approx%elems(ii,:,:) 
  END DO
  CLOSE(333)
  !DO ii=1,28
  !   WRITE(6,2) INT(approx%elems(1,ii,:))
!2    FORMAT (*(I3:x)) 
  !END DO
END PROGRAM MAIN
  
