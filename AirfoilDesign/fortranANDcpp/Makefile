    .SUFFIXES: .cpp .f90 .o

TARGET = main

    .f90.o:
	gfortran -c $<
    .cpp.o:
	g++ -c -fopenmp $<


    OBJECTS= mainLoop.o xfoil.o iopol.o  naca.o  spline.o  sort.o  userio.o  xbl.o  xblsys.o xgdes.o xgeom.o xoper.o xpanel.o xpol.o xsolve.o xutils.o 

    programs: $(OBJECTS)
	g++ $(OBJECTS)  -o $(TARGET)  -lgfortran -fopenmp 


clean:
	rm *.o
	$(RM) $(TARGET)



