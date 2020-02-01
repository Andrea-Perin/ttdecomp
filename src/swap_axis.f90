MODULE SWAP_AXIS
  USE TENSOR_TYPES

  IMPLICIT NONE

  ! INTERFACE FOR ALL TENSOR TYPES
  INTERFACE SWAP
     MODULE PROCEDURE SWAP_DTENSOR1
     MODULE PROCEDURE SWAP_DTENSOR2
     MODULE PROCEDURE SWAP_DTENSOR3 
     MODULE PROCEDURE SWAP_DTENSOR4 
  END INTERFACE SWAP

CONTAINS



  SUBROUTINE SWAP_DTENSOR1(tensor,nn,mm)
    !====================================================================
    !Does nothing.
    !INPUT:
    !- tensor (I/O)		: (DTENSOR1) the input tensor
    !- nn				: (INTEGER*4) nothing
    !- mm				: (INTEGER*4) nothing
    !====================================================================
    ! INPUT ARGUMENTS
    TYPE(DTENSOR1), INTENT(INOUT) :: tensor
    INTEGER*4, INTENT(IN) :: nn, mm
    RETURN
  END SUBROUTINE SWAP_DTENSOR1


  SUBROUTINE SWAP_DTENSOR2(tensor,nn,mm)
    !====================================================================
    !Modifies tensor, swapping the nn-th index with the mm-th.
    !INPUT:
    !- tensor (I/O)		: (DTENSOR2) the input tensor
    !- nn				: (INTEGER*4) the mode to swap with mm
    !- mm				: (INTEGER*4) the mode to swap with nn
    !====================================================================
    ! INPUT ARGUMENTS
    TYPE(DTENSOR2), INTENT(INOUT) :: tensor
    INTEGER*4, INTENT(IN) :: nn, mm
    ! UTILITY VARIABLES
    INTEGER*4 :: ii
    ! FUNCTION DEFINITION
    IF ((mm.LE.2).AND.(mm.GE.0).AND.(nn.LE.2).AND.(nn.GE.0)) THEN
       IF (mm.NE.nn) THEN
          tensor%modes = (/ tensor%modes(2),tensor%modes(1) /)		
          tensor%elems = TRANSPOSE(tensor%elems)
       END IF
    ELSE
       WRITE(*,*) "ERROR (SWAP_DTENSOR2): swapped modes must be in [1,2]. Found ",nn,mm
    END IF
    RETURN
  END SUBROUTINE SWAP_DTENSOR2



  SUBROUTINE SWAP_DTENSOR3(tensor,nn,mm)
    !====================================================================
    !Modifies tensor, swapping the nn-th index with the mm-th.
    !INPUT:
    !- tensor (I/O)		: (DTENSOR3) the input tensor
    !- nn				: (INTEGER*4) the mode to swap with mm
    !- mm				: (INTEGER*4) the mode to swap with nn
    !====================================================================
    ! INPUT ARGUMENTS
    TYPE(DTENSOR3), INTENT(INOUT) :: tensor
    INTEGER*4, INTENT(IN) :: nn, mm
    ! UTILITY VARIABLES
    INTEGER*4 :: ii, newmodes(3), neworder(3)
    ! FUNCTION DEFINITION
    IF ((mm.LE.3).AND.(mm.GE.0).AND.(nn.LE.3).AND.(nn.GE.0)) THEN
       newmodes = tensor%modes
       newmodes(nn) = tensor%modes(mm)		
       newmodes(mm) = tensor%modes(nn)
       neworder = (/(ii, ii=1,SIZE(tensor%modes),1)/)
       neworder(nn)=mm
       neworder(mm)=nn
       tensor%elems = RESHAPE(tensor%elems, shape=newmodes, order=neworder)
    ELSE
       WRITE(*,*) "ERROR (SWAP_DTENSOR3): swapped modes must be in [1,2,3]. Found ",nn,mm
    END IF
    RETURN
  END SUBROUTINE SWAP_DTENSOR3



  SUBROUTINE SWAP_DTENSOR4(tensor,nn,mm)
    !====================================================================
    !Modifies tensor, swapping the nn-th index with the mm-th.
    !INPUT:
    !- tensor (I/O)		: (DTENSOR4) the input tensor
    !- nn				: (INTEGER*4) the mode to swap with mm
    !- mm				: (INTEGER*4) the mode to swap with nn
    !====================================================================
    ! INPUT ARGUMENTS
    TYPE(DTENSOR4), INTENT(INOUT) :: tensor
    INTEGER*4, INTENT(IN) :: nn, mm
    ! UTILITY VARIABLES
    INTEGER*4 :: ii, newmodes(4), neworder(4)
    ! FUNCTION DEFINITION
    IF ((mm.LE.4).AND.(mm.GE.0).AND.(nn.LE.4).AND.(nn.GE.0)) THEN
       newmodes = tensor%modes
       newmodes(nn) = tensor%modes(mm)		
       newmodes(mm) = tensor%modes(nn)
       neworder = (/(ii, ii=1,SIZE(tensor%modes),1)/)
       neworder(nn)=mm
       neworder(mm)=nn
       tensor%elems = RESHAPE(tensor%elems, shape=newmodes, order=neworder)
    ELSE
       WRITE(*,*) "ERROR (SWAP_DTENSOR4): swapped modes must be in [1,2,3,4]. Found ",nn,mm
    END IF
    RETURN
  END SUBROUTINE SWAP_DTENSOR4


END MODULE SWAP_AXIS
