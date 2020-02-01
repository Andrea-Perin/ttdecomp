!------------------------PRODUCTS MODULE
MODULE PRODUCTS
  IMPLICIT NONE
  INTERFACE KRONECKER
     MODULE PROCEDURE KRONECKER_M
     MODULE PROCEDURE KRONECKER_V
  END INTERFACE KRONECKER
  INTERFACE SHOW
     MODULE PROCEDURE SHOW_M
     MODULE PROCEDURE SHOW_V
  END INTERFACE SHOW
CONTAINS
  
  !----------------------KRONECKER MATRIX SUBROUTINE
  SUBROUTINE KRONECKER_M(AA,BB,AB)
    REAL*8::AA(:,:),BB(:,:)
    REAL*8, ALLOCATABLE::AB(:,:)
    INTEGER::RR,RA,RB,CC,CA,CB,II,JJ
    RA = SIZE(AA,1)
    CA = SIZE(AA,2)
    RB = SIZE(BB,1)
    CB = SIZE(BB,2)
    IF (ALLOCATED(AB)) DEALLOCATE(AB)	!Discard any lingering storage.
    ALLOCATE (AB(RA*RB,CA*CB))		!Obtain the exact desired size.
    RR = 0		!Syncopation: start the row offset.
    DO II = 1,RA	!Step down the rows of A.
       CC = 0		!For each row, start the column offset.
       DO JJ = 1,CA		!Step along the columns of A.
          AB(RR + 1:RR + RB,CC + 1:CC + CB) = AA(II,JJ)*BB	!Place a block of B values.
          CC = CC + CB		!Advance a block of columns.
       END DO		!On to the next column of A.
       RR = RR + RB		!Advance a block of rows.
    END DO	!On to the next row of A.
  END SUBROUTINE KRONECKER_M
  
  !----------------------KRONECKER VECTOR SUBROUTINE
  SUBROUTINE KRONECKER_V(AA,BB,AB)
    REAL*8::AA(:),BB(:)
    REAL*8, ALLOCATABLE::AB(:)
    INTEGER::LL,LA,LB,II
    LA = SIZE(AA)
    LB = SIZE(BB)
    IF (ALLOCATED(AB)) DEALLOCATE(AB)
    ALLOCATE (AB(LA*LB))
    LL = 0
    DO II = 1,LA
       AB(LL + 1:LL + LB) = AA(II)*BB
       LL = LL + LB
    END DO
  END SUBROUTINE KRONECKER_V
  
  !----------------------SHOW MATRIX SUBROUTINE
  SUBROUTINE SHOW_M(FF,AA)	!Write array A in row,column order.
    INTEGER::FF	!Output file unit number.
    REAL*8::AA(:,:)	!The 2-D array, lower bound one.
    INTEGER::RR	!The row stepper.
    DO RR = 1,SIZE(AA,1)	!Each row gets its own line.
       WRITE (FF,1) AA(RR,:)		!Write all the columns of that row.
1      FORMAT (*(G10.3:x))		!This suffices for the example.
    END DO			!On to the next row.
  END SUBROUTINE SHOW_M	!WRITE (F,*) A or similar would show A as if transposed.
  
  !----------------------SHOW VECTOR SUBROUTINE
  SUBROUTINE SHOW_V(FF,AA)
    INTEGER::FF
    REAL*8::AA(:)
    WRITE (FF,1) AA
1   FORMAT (*(G10.3:x))
  END SUBROUTINE SHOW_V
  
  !----------------------KHATRI RAO SUBROUTINE
  SUBROUTINE KHATRI_RAO(AA,BB,AB)
    REAL*8::AA(:,:),BB(:,:)
    REAL*8, ALLOCATABLE::AB(:,:),COL(:)
    INTEGER::RA,RB,II,CA,CB
    RA = SIZE(AA,1)
    CA = SIZE(AA,2)
    RB = SIZE(BB,1)
    CB = SIZE(BB,2)
    IF (CA==CB) THEN
       IF (ALLOCATED(AB)) DEALLOCATE(AB)
       ALLOCATE (AB(RA*RB,CA))
       ALLOCATE(COL(CA))
       DO II = 1,CA
          CALL KRONECKER_V(AA(:,II),BB(:,II),COL)
          AB(:,II)=COL
       END DO
    ELSE
       PRINT*, "Dimension error in khatri rao product"
    END IF
  END SUBROUTINE KHATRI_RAO

  !----------------------GAUSSIAN RANDOM INITIALIZATION SUBROUTINE
  SUBROUTINE RAND_INIT(mat,sigma)
    REAL*8, ALLOCATABLE::mat(:,:)
    REAL*8::uni_1,uni_2,norm_1,norm_2,sigma
    INTEGER::ii,jj,dim1,dim2,new_dim1
    dim1=SIZE(mat,1)
    dim2=SIZE(mat,2)
    new_dim1=dim1
    IF (MOD(new_dim1,2).EQ.1) THEN
       new_dim1=new_dim1-1
    END IF
    DO jj=1,dim2
       DO ii=1,new_dim1,2
          CALL RANDOM_NUMBER(uni_1)
          CALL RANDOM_NUMBER(uni_2)
          mat(ii,jj)  =sigma*SQRT(-2*LOG(uni_1))*COS(8.D0*ATAN(1.D0)*uni_2)
          mat(ii+1,jj)=sigma*SQRT(-2*LOG(uni_1))*SIN(8.D0*ATAN(1.D0)*uni_2)
       END DO
       IF (MOD(dim1,2).EQ.1) THEN
          CALL RANDOM_NUMBER(uni_1)
          CALL RANDOM_NUMBER(uni_2)
          mat(dim1,jj)=sigma*SQRT(-2*LOG(uni_1))*COS(8.D0*ATAN(1.D0)*uni_2)
       END IF
    END DO
  END SUBROUTINE RAND_INIT
  
END MODULE PRODUCTS
