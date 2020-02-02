PROGRAM MAIN

  USE DEBUGGER
  USE TENSOR_TYPES
  USE MODE_N
  USE SWAP_AXIS
  USE MAT_UTILS
  
  IMPLICIT NONE

  INTEGER,ALLOCATABLE::dim(:)
  REAL*8,ALLOCATABLE::IN(:,:)
  INTEGER::nn,rows,cols,kk,ii
  TYPE(DTENSOR3)::my_tens
  OPEN(332,file='../data/mnist_1k.csv',status="old",action="read")
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
  WRITE (6,*) nn
  WRITE (6,*) "info:",dim
  WRITE (6,*) "MATRIX:"
  CALL TENSOR3(dim,IN,my_tens)
  DO ii=1,dim(2)
     WRITE(6,1) my_tens%elems(1,ii,:)
1    FORMAT (*(F4.0:x)) 
  END DO
  
END PROGRAM MAIN
