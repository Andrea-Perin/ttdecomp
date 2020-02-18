PROGRAM MAIN

  USE DEBUGGER
  USE TENSOR_TYPES
  USE MODE_N
  USE SWAP_AXIS
  USE MAT_UTILS
  
  IMPLICIT NONE

  INTEGER,ALLOCATABLE::dimm(:)
  REAL*8,ALLOCATABLE::IN(:,:)
  INTEGER,ALLOCATABLE::vec(:)
  INTEGER::nn,rows,cols,kk,ii
  TYPE(DTENSOR3)::my_tens,copy
  !OPEN(332,file='../data/mnist_1k.csv',status="old",action="read")
  OPEN(332,file='../data/prova.csv',status="old",action="read")
  READ (332, *) nn
  ALLOCATE(dimm(nn))
  READ(332,*) dimm
  rows=dimm(1)
  cols=1
  DO kk=1,nn-1
     cols=cols*dimm(1+kk)
  END DO
  ALLOCATE(IN(cols,rows))
  READ(332,*) IN
  IN=TRANSPOSE(IN)
  CLOSE(332)
  WRITE (6,*) nn
  WRITE (6,*) "info:",dimm
  WRITE (6,*) "MATRIX:"
  CALL TENSOR3(dimm,IN,my_tens,3) ! assuming mode 1
  DO ii=1,dimm(2)
     WRITE(6,1) my_tens%elems(1,ii,:)
1    FORMAT (*(F4.0:x)) 
  END DO
  WRITE (6,*) "dimensioni del tensore"
  ALLOCATE(vec(3))
  vec = SHAPE(my_tens%elems)
  WRITE (6,*) vec
  ALLOCATE(copy%elems(vec(1),vec(2),vec(3)))
  copy%elems=my_tens%elems
  DO ii=1,dimm(2)
     WRITE(6,1) copy%elems(1,ii,:) 
  END DO
END PROGRAM MAIN
