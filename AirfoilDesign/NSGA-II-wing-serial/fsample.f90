SUBROUTINE mysub(OBJS, NOBJ, VARS, NVAR)
	implicit none
	integer,intent(in):: NVAR, NOBJ
	real*8, intent(out):: OBJS(NOBJ), VARS(NVAR)
	integer::i
	real:: g, PI
	! Here we do create the ZDT3 example for testing the C++ <=> Fortran communication 	
	g = 0
	PI = 3.14159265
	OBJS(1) = VARS(1)
	do 111 i = 2, NVAR
		g = g + VARS(i)
	111 continue
	g = 1 + 9*g/(NVAR-1);
	OBJS(2) = g*(1 - sqrt(VARS(1)/g) - VARS(1)*sin(10*PI*VARS(1))/g);
!	print*,VARS(1),', ', VARS(2),', ', VARS(3),', ', VARS(4),', ', VARS(5) 
!	print*,'OBJ(1)=',OBJS(1), ', OBJ(2)=', OBJS(2)
!print*,"Hi from Fortran"
END SUBROUTINE mysub
