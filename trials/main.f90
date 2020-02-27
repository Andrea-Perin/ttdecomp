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

  ! INPUT VARIABLES
  INTEGER,ALLOCATABLE :: dimm(:)
  REAL*8,ALLOCATABLE :: IN(:,:)
  INTEGER :: nn,rows,cols,ii
  ! FOR MNIST
  TYPE(DTENSOR3) :: MNIST, approx_MNIST, core_MNIST
  ! FOR LANDSCAPE
  TYPE(DTENSOR3) :: land, approx_land, core_land
  ! FOR VIDEO
  TYPE(DTENSOR4) :: video, approx_video, core_video
  ! DECOMPOSITION VARIABLES
  TYPE(matrix_list),ALLOCATABLE :: lista(:)
  TYPE(tensor_list),ALLOCATABLE :: tens_lista(:)
  REAL*8 :: error
  REAL*8, ALLOCATABLE :: lambdas(:)
  ! VARIED WITH PYTHON SCRIPT
  CHARACTER(LEN=64) :: in_file,out_file
  INTEGER*4 :: choose_file, choose_method
  INTEGER*4 :: num_images, x_pixels, y_pixels
  REAL*8 :: epsilon
  INTEGER*4 :: rango
  INTEGER*4,ALLOCATABLE :: ranks(:)
  


  !============================================
  ! CHOOSE FILE
  ! - 1    : MNIST
  ! - 2    : landscape
  ! - 3    : video
  !============================================
  !choose_file = 2
  OPEN (126,file="choose_file.dat",status="old",action="read")
  READ (126, *) choose_file
  CLOSE(126)


  !============================================
  ! CHOOSE FILE RESOLUTION
  ! - num_images  : images in MNIST
  ! - x_pixels    : x pixels of landscape
  ! - Y_pixels    : y pixels of landscape
  !============================================
  !num_images = 100 ! 100, 250, 500, 1000
  !x_pixels = 112    ! 112, 149, 298, 478
  !y_pixels = 240    ! 240, 320, 640, 1024
  OPEN (127,file="file_resolution.dat",status="old",action="read")
  IF (choose_file.EQ.1) THEN
     ! MNIST (3D TENSOR)
     READ (127, *) num_images
  ELSEIF (choose_file.EQ.2) THEN
     ! LANDSCAPE IMAGE (3D TENSOR)
     READ (127, *) x_pixels, y_pixels
  END IF
  CLOSE(127)


  !============================================
  ! CHOOSE METHOD
  ! - 4    : MPS 
  ! - 5    : CP
  ! - 6    : HOSVD TUCKER
  ! - 7    : HOOI  TUCKER
  ! - 8    : HOOI  TUCKER RANDOM INITIALIZATION
  !============================================
  !choose_method = 4
  OPEN (128,file="choose_method.dat",status="old",action="read")
  READ (128, *) choose_method
  CLOSE(128)

  !============================================
  ! CHOOSE DECOMPOSITION PARAMETERS
  ! - epsilon     : tensor train epsilon
  ! - rango       : cpd rank
  ! - ranks(:)    : tucker ranks
  !============================================ 
  OPEN (129,file="decomposition_parameters.dat",status="old",action="read")
  IF (choose_method.EQ.4) THEN
     READ (129, *) epsilon
     !epsilon = 5D-2
  ELSEIF (choose_method.EQ.4) THEN
     READ (129, *) rango
     !rango = 16
  ELSEIF ((choose_method.GE.6).AND.(choose_method.LE.8)) THEN
     IF ((choose_file.EQ.1).OR.(choose_file.EQ.2)) THEN
        ALLOCATE(ranks(3))
        READ (129, *) ranks(:)
        !ranks = (/ 10,3,3 /)
     ELSEIF (choose_file.EQ.3) THEN
        ALLOCATE(ranks(4))
        READ (129, *) ranks(:)
        !ranks = (/ 64,64,64,3 /)
     END IF
  END IF
  CLOSE(129)
  
  
  !============================================
  ! SELECT FILE
  !============================================
  IF (choose_file.EQ.1) THEN
     ! MNIST (3D TENSOR)
     WRITE (in_file,"(A14,I0,A4)") "../data/mnist_",num_images,".csv"
  ELSEIF (choose_file.EQ.2) THEN
     ! LANDSCAPE IMAGE (3D TENSOR)
     WRITE (in_file,"(A13,I0,A1,I0,A4)") "../data/land_",x_pixels,"_",y_pixels,".csv"
  ELSEIF (choose_file.EQ.3) THEN
     ! SHORT VIDEO (4D TENSOR)
     WRITE (in_file,"(A22)") "../data/video_144p.csv"
  END IF
  
     
  !============================================
  ! READ FILE
  !============================================
  OPEN(332,file=in_file,status="old",action="read")
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

  !============================================
  ! STORE IN THE APPROPRIATE TENSOR 
  !============================================
  IF (choose_file.EQ.1) THEN
     ! MNIST (3D TENSOR)
     MNIST = TENSOR3(dimm,IN,1) ! assuming mode 1
     ALLOCATE(lista(3))
     ALLOCATE(tens_lista(3))
  ELSEIF (choose_file.EQ.2) THEN
     ! LANDSCAPE IMAGE (3D TENSOR)
     ALLOCATE(lista(3))
     ALLOCATE(tens_lista(3))
     land = TENSOR3(dimm,IN,1) ! assuming mode 1
  ELSEIF (choose_file.EQ.3) THEN
     ! SHORT VIDEO (4D TENSOR)
     ALLOCATE(lista(4))
     ALLOCATE(tens_lista(4))
     video = TENSOR4(dimm,IN,1) ! assuming mode 1
  END IF
  
  
  !============================================
  ! COMPRESS AND RECONSTRUCT
  !============================================


  IF (choose_method.EQ.4) THEN
     ! MPS DECOMPOSITION
     OPEN (130,file="decomposition_parameters.dat",status="old",action="write")
     IF (choose_file.EQ.1) THEN
        ! MNIST (3D TENSOR)
        CALL MPS(MNIST,tens_lista,eps=epsilon)
        approx_MNIST = MPS_TO_TENSOR3(tens_lista)
        WRITE (out_file,"(A14,I0,A5,I0,A1,I0,A1,I0,A1,I0,A4)") &
             "../data/mnist_",num_images,"_MPS_",tens_lista(1)%cores%modes(1),"_",tens_lista(2)%cores%modes(1),&
             "_",tens_lista(3)%cores%modes(1),"_",tens_lista(3)%cores%modes(3),".csv"
        WRITE(130,*) tens_lista(1)%cores%modes(1),tens_lista(2)%cores%modes(1),&
             tens_lista(3)%cores%modes(1),tens_lista(3)%cores%modes(3)
     ELSEIF (choose_file.EQ.2) THEN
        ! LANDSCAPE IMAGE (3D TENSOR)
        CALL MPS(land,tens_lista,eps=epsilon)
        approx_land = MPS_TO_TENSOR3(tens_lista)
        WRITE (out_file,"(A13,I0,A1,I0,A5,I0,A1,I0,A1,I0,A1,I0,A4)") "../data/land_",x_pixels,"_",y_pixels,"_MPS_",&
             tens_lista(1)%cores%modes(1),"_",tens_lista(2)%cores%modes(1),"_",tens_lista(3)%cores%modes(1),&
             "_",tens_lista(3)%cores%modes(3),".csv"
        WRITE(130,*) tens_lista(1)%cores%modes(1),tens_lista(2)%cores%modes(1),&
             tens_lista(3)%cores%modes(1),tens_lista(3)%cores%modes(3)
     ELSEIF (choose_file.EQ.3) THEN
        ! SHORT VIDEO (4D TENSOR)
        CALL MPS(video,tens_lista,eps=epsilon)
        approx_video = MPS_TO_TENSOR4(tens_lista)
        WRITE (out_file,"(A23,I0,A1,I0,A1,I0,A1,I0,A1,I0,A4)") "../data/video_144p_MPS_",tens_lista(1)%cores%modes(1),&
             "_",tens_lista(2)%cores%modes(1),"_",tens_lista(3)%cores%modes(1),"_",tens_lista(4)%cores%modes(1),&
             "_",tens_lista(4)%cores%modes(3),".csv"
        WRITE(130,*) tens_lista(1)%cores%modes(1),tens_lista(2)%cores%modes(1),&
             tens_lista(3)%cores%modes(1),tens_lista(4)%cores%modes(1),tens_lista(4)%cores%modes(3)
     END IF
     CLOSE(130)
  ELSEIF (choose_method.EQ.5) THEN
     ! CP DECOMPOSITION
     rango=16
     IF (choose_file.EQ.1) THEN
        ! MNIST (3D TENSOR)
        CALL CPD(MNIST,rango,lista,lambdas,error)
        core_MNIST = TENSOR3((/rango,rango,rango/),TO_IDENTITY(lambdas,SIZE(MNIST%modes)),1)
        approx_MNIST = RECO(core_MNIST,lista)
        WRITE (out_file,"(A14,I0,A5,I0,A4)") "../data/mnist_",num_images,"_cpd_",rango,".csv"
     ELSEIF (choose_file.EQ.2) THEN
        ! LANDSCAPE IMAGE (3D TENSOR)
        CALL CPD(land,rango,lista,lambdas,error)
        core_land = TENSOR3((/rango,rango,rango/),TO_IDENTITY(lambdas,SIZE(land%modes)),1)
        approx_land = RECO(core_land,lista)
        WRITE (out_file,"(A13,I0,A1,I0,A5,I0,A4)") "../data/land_",x_pixels,"_",y_pixels,"_cpd_",rango,".csv"
     ELSEIF (choose_file.EQ.3) THEN
        ! SHORT VIDEO (4D TENSOR)
        CALL CPD(video,rango,lista,lambdas,error)
        core_video = TENSOR4((/rango,rango,rango,rango/),TO_IDENTITY(lambdas,SIZE(video%modes)),1)
        approx_video = RECO(core_video,lista)
        WRITE (out_file,"(A23,I0,A4)") "../data/video_144p_cpd_",rango,".csv"
     END IF
  ELSEIF (choose_method.EQ.6) THEN
     ! HOSVD TUCKER DECOMPOSITION
     IF (choose_file.EQ.1) THEN
        ! MNIST (3D TENSOR)
        CALL HOSVD(MNIST,ranks,core_MNIST,lista)
        approx_MNIST = RECO(core_MNIST,lista)
        WRITE (out_file,"(A14,I0,A7,I0,A1,I0,A1,I0,A4)") &
             "../data/mnist_",num_images,"_hosvd_",ranks(1),"_",ranks(2),"_",ranks(3),".csv"
     ELSEIF (choose_file.EQ.2) THEN
        ! LANDSCAPE IMAGE (3D TENSOR)
        CALL HOSVD(land,ranks,core_land,lista)
        approx_land = RECO(core_land,lista)
        WRITE (out_file,"(A13,I0,A1,I0,A7,I0,A1,I0,A1,I0,A4)") &
             "../data/land_",x_pixels,"_",y_pixels,"_hosvd_",ranks(1),"_",ranks(2),"_",ranks(3),".csv"
     ELSEIF (choose_file.EQ.3) THEN
        ! SHORT VIDEO (4D TENSOR)
        CALL HOSVD(video,ranks,core_video,lista)
        approx_video = RECO(core_video,lista)
        WRITE (out_file,"(A25,I0,A1,I0,A1,I0,A4)") &
             "../data/video_144p_hosvd_",ranks(1),"_",ranks(2),"_",ranks(3),".csv"
     END IF
  ELSEIF (choose_method.EQ.7) THEN
     ! HOOI TUCKER DECOMPOSITION
     IF (choose_file.EQ.1) THEN
        ! MNIST (3D TENSOR)
        CALL HOOI(MNIST,ranks,core_MNIST,lista,error)
        approx_MNIST = RECO(core_MNIST,lista)
        WRITE (out_file,"(A14,I0,A6,I0,A1,I0,A1,I0,A4)") &
             "../data/mnist_",num_images,"_hooi_",ranks(1),"_",ranks(2),"_",ranks(3),".csv"
     ELSEIF (choose_file.EQ.2) THEN
        ! LANDSCAPE IMAGE (3D TENSOR)
        CALL HOOI(land,ranks,core_land,lista,error)
        approx_land = RECO(core_land,lista)
        WRITE (out_file,"(A13,I0,A1,I0,A6,I0,A1,I0,A1,I0,A4)") &
             "../data/land_",x_pixels,"_",y_pixels,"_hooi_",ranks(1),"_",ranks(2),"_",ranks(3),".csv"
     ELSEIF (choose_file.EQ.3) THEN
        ! SHORT VIDEO (4D TENSOR)
        CALL HOOI(video,ranks,core_video,lista,error)
        approx_video = RECO(core_video,lista)
        WRITE (out_file,"(A24,I0,A1,I0,A1,I0,A4)") &
             "../data/video_144p_hooi_",ranks(1),"_",ranks(2),"_",ranks(3),".csv"
     END IF
  ELSEIF (choose_method.EQ.8) THEN
     ! HOOI TUCKER DECOMPOSITION
     IF (choose_file.EQ.1) THEN
        ! MNIST (3D TENSOR)
        CALL HOOI(MNIST,ranks,core_MNIST,lista,error,randinit=.TRUE.)
        approx_MNIST = RECO(core_MNIST,lista)
        WRITE (out_file,"(A14,I0,A8,I0,A1,I0,A1,I0,A4)") &
             "../data/mnist_",num_images,"_random_",ranks(1),"_",ranks(2),"_",ranks(3),".csv"
     ELSEIF (choose_file.EQ.2) THEN
        ! LANDSCAPE IMAGE (3D TENSOR)
        CALL HOOI(land,ranks,core_land,lista,error,randinit=.TRUE.)
        approx_land = RECO(core_land,lista)
        WRITE (out_file,"(A13,I0,A1,I0,A8,I0,A1,I0,A1,I0,A4)") &
             "../data/land_",x_pixels,"_",y_pixels,"_random_",ranks(1),"_",ranks(2),"_",ranks(3),".csv"
     ELSEIF (choose_file.EQ.3) THEN
        ! SHORT VIDEO (4D TENSOR)
        CALL HOOI(video,ranks,core_video,lista,error,randinit=.TRUE.)
        approx_video = RECO(core_video,lista)
        WRITE (out_file,"(A26,I0,A1,I0,A1,I0,A4)") &
             "../data/video_144p_random_",ranks(1),"_",ranks(2),"_",ranks(3),".csv"
     END IF
  END IF

  !============================================
  ! SAVE DATA ON FILE
  !============================================
  ! OPEN(333,file=out_file,status="unknown",action="write")
  ! IF (choose_file.EQ.1) THEN
  !    ! MNIST (3D TENSOR)
  !    DO ii=1,SIZE(approx_MNIST%elems, 1)
  !       WRITE(333,*) approx_MNIST%elems(ii,:,:) 
  !    END DO
  ! ELSEIF (choose_file.EQ.2) THEN
  !    ! LANDSCAPE IMAGE (3D TENSOR)
  !    DO ii=1,SIZE(approx_land%elems, 1)
  !       WRITE(333,*) approx_land%elems(ii,:,:) 
  !    END DO
  ! ELSEIF (choose_file.EQ.3) THEN
  !    ! SHORT VIDEO (4D TENSOR)
  !    DO ii=1,SIZE(approx_video%elems, 1)
  !       WRITE(333,*) approx_video%elems(ii,:,:,:) 
  !    END DO
  ! END IF
  ! CLOSE(333)

  !============================================
  ! SAVE ERROR ON FILE
  !============================================
  OPEN(334,file="../data/error.dat",status="unknown",action="write")
  IF (choose_file.EQ.1) THEN
     ! MNIST (3D TENSOR)
     WRITE(334,*) SQRT(SUM((MNIST%elems-approx_MNIST%elems)**2))/SIZE(MNIST%elems)
  ELSEIF (choose_file.EQ.2) THEN
     ! LANDSCAPE IMAGE (3D TENSOR)
     WRITE(334,*) SQRT(SUM((land%elems-approx_land%elems)**2))/SIZE(land%elems)
  ELSEIF (choose_file.EQ.3) THEN
     ! SHORT VIDEO (4D TENSOR)
     WRITE(334,*) SQRT(SUM((video%elems-approx_video%elems)**2))/SIZE(video%elems)
  END IF
  CLOSE(334)

END PROGRAM MAIN
