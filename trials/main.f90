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
  INTEGER::nn,rows,cols,ii
  ! FOR VIDEO
  TYPE(DTENSOR4):: video, approx_video, core_video
  ! FOR LANDSCAPE
  TYPE(DTENSOR3) :: landscape, approx_landscape, core_landscape
  ! FOR MNIST
  TYPE(DTENSOR3) :: MNIST, approx_MNIST, core_MNIST
  ! OTHER
  TYPE(matrix_list),ALLOCATABLE :: lista(:)
  TYPE(tensor_list),ALLOCATABLE :: tens_lista(:)
  
  INTEGER*4 :: rango
  REAL*8 :: threshold, error
  REAL*8, ALLOCATABLE :: lambdas(:)

  CHARACTER(LEN=64)::in_file,out_file
  INTEGER*4 :: choose_file, choose_method
  INTEGER*4 :: num_images, x_pixels, y_pixels
  


  !======================================
  ! CHOOSE FILE
  ! - 1    : MNIST
  ! - 2    : landscape
  ! - 3    : video
  !======================================
  choose_file = 2


  !======================================
  ! CHOOSE FILE RESOLUTION
  ! - num_images  : images in MNIST
  ! - x_pixels    : x pixels of landscape
  ! - Y_pixels    : y pixels of landscape
  !======================================
  num_images = 1000
  x_pixels = 112
  y_pixels = 240
  
  
  !======================================
  ! SELECT FILE
  !======================================  
  IF (choose_file.EQ.1) THEN
     ! MNIST 1K (3D TENSOR)
     WRITE (in_file,"(A14,I0,A4)") "../data/mnist_",num_images,".csv"
  ELSEIF (choose_file.EQ.2) THEN
     ! LANDSCAPE IMAGE (3D TENSOR)
     WRITE (in_file,"(A13,I0,A1,IO,A4)") "../data/land_",x_pixels,"_",y_pixels,".csv"
  ELSEIF (choose_file.EQ.3) THEN
     ! SHORT VIDEO (4D TENSOR)
     WRITE (in_file,"(A25)") "../data/original_144p.csv"
  END IF

  OPEN(332,file=in_file,status="old",action="read")
     
  !======================================
  ! READ FILE
  !======================================
  READ (332, *) nn
  ALLOCATE(dimm(nn))
  READ(332,*) dimm
  rows=dimm(1)
  cols=1
  DO ii=2,nn
     cols=cols*dimm(ii)
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
     ALLOCATE(lista(3))
     ALLOCATE(tens_lista(3))
  ELSEIF (choose_file.EQ.2) THEN
     ALLOCATE(lista(3))
     ALLOCATE(tens_lista(3))
     landscape = TENSOR3(dimm,IN,1) ! assuming mode 1
  ELSEIF (choose_file.EQ.3) THEN
     ALLOCATE(lista(4))
     ALLOCATE(tens_lista(4))
     video = TENSOR4(dimm,IN,1) ! assuming mode 1
  END IF

  !======================================
  ! CHOOSE METHOD
  ! - 4    : MPS 
  ! - 5    : CP
  ! - 6    : TUCKER
  !======================================
  choose_method = 5
  
  
  !======================================
  ! COMPRESS RECONSTRUCT AND SAVE FILE
  !======================================


  ! IF (choose_method.EQ.4) THEN
  !    !MPS DECOMPOSITION
  !    IF (choose_file.EQ.1) THEN
  !       CALL MPS(MNIST,tens_lista,eps=5D-2)
  !    ELSEIF (choose_file.EQ.2) THEN
  !       CALL MPS(landscape,tens_lista,eps=5D-2)
  !    ELSEIF (choose_file.EQ.3) THEN
  !       CALL MPS(video,tens_lista,eps=5D-2)
  !    END IF
  !    DO ii=1,SIZE(tens_lista)
  !       print*, "Shape of core :",ii, tens_lista(ii)%cores%modes 
  !    END DO
  !    IF (choose_file.EQ.1) THEN
  !       approx_MNIST = MPS_TO_TENSOR3(tens_lista)
  !       WRITE (out_file,"(A4,I0,A4)") "../data/mnist_",N,".csv"
  !    ELSEIF (choose_file.EQ.2) THEN
  !       approx_landscape = MPS_TO_TENSOR3(tens_lista)
  !    ELSEIF (choose_file.EQ.3) THEN
  !       approx_video = MPS_TO_TENSOR4(tens_lista)
  !    END IF
     
  !    !SAVE ON FILE
  !    OPEN(333,file='../data/land_149_320_MPS_eps_5.csv',status="unknown",action="write")
  !    DO ii=1,SIZE(approx_landscape%elems, 1)
  !       WRITE(333,*) approx_landscape%elems(ii,:,:) 
  !    END DO
  !    CLOSE(333)
     

  ! ELSEIF (choose_method.EQ.5) THEN
  !    landscape = TENSOR3(dimm,IN,1) ! assuming mode 1
  ! ELSEIF (choose_method.EQ.6) THEN
  !    video = TENSOR4(dimm,IN,1) ! assuming mode 1
  ! END IF

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
  CALL MPS(landscape,tens_lista,eps=5D-2)  ! roughly half the parameters for 1% error
  ! sizes match the theory (at least)
  print*, "Shape of core 1:", tens_lista(1)%cores%modes 
  print*, "Shape of core 2:", tens_lista(2)%cores%modes
  print*, "Shape of core 3:", tens_lista(3)%cores%modes
  ! try reconstruction
  approx_landscape = MPS_TO_TENSOR4(tens_lista)
  !print*, "Approx size:", SIZE(approx_landscape%elems)
  !print*, "True size:", PRODUCT(landscape%modes)
  ! show error
  !print*, SQRT(SUM((landscape%elems-approx_landscape%elems)**2))/SQRT(SUM(approx_landscape%elems**2))
  !SAVE ON FILE
  OPEN(333,file='../data/land_149_320_MPS_eps_5.csv',status="unknown",action="write")
  DO ii=1,SIZE(approx_landscape%elems, 1)
    WRITE(333,*) approx_landscape%elems(ii,:,:) 
  END DO
  CLOSE(333)

END PROGRAM MAIN
