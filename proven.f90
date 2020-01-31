PROGRAM MAIN
	USE TENSOR_TYPES
	USE MODE_N
	USE SWAP_AXIS
	USE PARTINV
	IMPLICIT NONE

	TYPE(DTENSOR3) :: mammamia
	TYPE(DTENSOR2) :: mammatua
	REAL*8 :: test1(4,3), mask(3,4)
	REAL*8, ALLOCATABLE :: UU(:,:), SIG(:), VV(:,:)
    INTEGER*4 :: ii,jj 

	test1 = reshape((/1,2,3,6,8,-8,1,1,-1,2,2,5/),shape(test1),order=(/2,1/))
	mask = PINV(test1)
	do ii=1,SIZE(mask,1)
		print*, mask(ii,:)
	end do



END PROGRAM
