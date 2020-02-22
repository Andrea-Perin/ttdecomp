PROGRAM MAIN

  USE DEBUGGER
  USE TENSOR_TYPES
  USE MODE_N
  USE SWAP_AXIS
  USE MAT_UTILS
  USE CPD_UTILS
  USE TUCKER
  USE MPS_MOD
  
  IMPLICIT NONE

  INTEGER,ALLOCATABLE::dimm(:)
  REAL*8,ALLOCATABLE::IN(:,:)
  INTEGER,ALLOCATABLE::ranks(:)
  INTEGER::nn,rows,cols,kk,ii
  ! FOR VIDEO
  TYPE(DTENSOR4):: video, approx_video, core_video
  ! FOR LANDSCAPE
  TYPE(DTENSOR3) :: landscape, approx_landscape, core_landscape
  ! FOR MNIST
  TYPE(DTENSOR3) :: MNIST, approx_MNIST, core_MNIST
  ! OTHER
  TYPE(matrix_list) :: lista(3)
  TYPE(tensor_list) :: tens_lista(3)
  
  INTEGER*4 :: rango
  REAL*8 :: threshold, error
  REAL*8, ALLOCATABLE :: lambdas(:)

  INTEGER*4 :: choose_file


  !======================================
  ! CHOOSE FILE
  ! - 1    : MNIST
  ! - 2    : landscape
  ! - 3    : video
  !======================================
  choose_file = 2

  
  !======================================
  ! SELECT FILE
  !======================================  
  IF (choose_file.EQ.1) THEN
     ! MNIST 1K (3D TENSOR)
     OPEN(332,file='../data/mnist_1k.csv',status="old",action="read")
  ELSEIF (choose_file.EQ.2) THEN
     ! LANDSCAPE IMAGE (3D TENSOR)
     OPEN(332,file='../data/land_112_240.csv',status="old",action="read")
  ELSEIF (choose_file.EQ.3) THEN
     ! SHORT VIDEO (4D TENSOR)
     OPEN(332,file='../data/original_144p.csv',status="old",action="read")
  END IF
     
  !======================================
  ! READ FILE
  !======================================
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

  !======================================
  ! STORE IN THE APPROPRIATE TENSOR 
  !======================================
  IF (choose_file.EQ.1) THEN
     MNIST = TENSOR3(dimm,IN,1) ! assuming mode 1
  ELSEIF (choose_file.EQ.2) THEN
     landscape = TENSOR3(dimm,IN,1) ! assuming mode 1
  ELSEIF (choose_file.EQ.3) THEN
     video = TENSOR4(dimm,IN,1) ! assuming mode 1
  END IF
  
  !======================================
  ! COMPRESS
  !======================================

  ! !TUCKER DECOMPOSITION
  ! ALLOCATE(ranks(4))
  ! ranks = (/ 64,64,64,3 /)
  ! CALL HOSVD(my_tens,ranks,core,lista)
  ! !CALL HOOI4(my_tens,ranks,core,lista,error)!,thresh=1D-10,randinit=.TRUE.)
  ! approx = RECO(core,lista)
  ! !OPEN(333,file='../data/land_112_240_hosvd_10_20_3.csv',status="unknown",action="write")
  ! OPEN(333,file='../data/original_144p_hosvd_64_64_64_3.csv',status="unknown",action="write")
  ! DO ii=1,dimm(1)
  !    WRITE(333,*) approx%elems(ii,:,:,:) 
  ! END DO
  ! CLOSE(333)

  ! !CP DECOMPOSITION
  ! rango=16
  ! CALL CPD(my_tens,rango,lista,lambdas,error,numiter=2)
  ! !RECONSTRUCT
  ! core = TENSOR4( (/rango,rango,rango,rango/) , TO_IDENTITY(lambdas,SIZE(my_tens%modes)) , 1 )
  ! approx = RECO4(core,lista)
  ! !SAVE ON FILE
  ! !OPEN(333,file='../data/land_112_240_cpd_10_iter_2.csv',status="unknown",action="write")
  ! OPEN(333,file='../data/original_144p_cpd_16_iter_2.csv',status="unknown",action="write")
  ! DO ii=1,dimm(1)
  !   WRITE(333,*) approx%elems(ii,:,:,:) 
  ! END DO
  ! CLOSE(333)

  !MPS DECOMPOSITION 
  CALL MPS(landscape,tens_lista,eps=1D-2)  ! roughly half the parameters for 1% error
  ! sizes match the theory (at least)
  print*, "Shape of core 1:", tens_lista(1)%cores%modes 
  print*, "Shape of core 2:", tens_lista(2)%cores%modes
  print*, "Shape of core 3:", tens_lista(3)%cores%modes
  ! try reconstruction
  approx_landscape = MPS_TO_TENSOR3(tens_lista)
  print*, "Approx size:", SIZE(approx_landscape%elems)
  print*, "True size:", PRODUCT(landscape%modes)
  ! show error
  print*, SQRT(SUM((landscape%elems-approx_landscape%elems)**2))/SQRT(SUM(approx_landscape%elems**2))
  !SAVE ON FILE
  OPEN(333,file='../data/land_112_240_MPS.csv',status="unknown",action="write")
  DO ii=1,SIZE(approx_landscape%elems, 1)
    WRITE(333,*) approx_landscape%elems(ii,:,:) 
  END DO
  CLOSE(333)

END PROGRAM MAIN
  
