# Makefile for a trial file

# Defining variables
PROGNAME = main
F90 = gfortran
OFLAGS = -O0
COMPFLAGS = -llapack -lblas
DEBFLAGS = -g #-Wall
OBJ_DIR = ../objects
SRC_DIR = ../src
OBJS = debugger.o tensor_types.o mode_n.o swap_axis.o mat_utils.o cpd_utils.o



# Actual filemaking
$(PROGNAME).x: $(OBJS) $(PROGNAME).o
	$(F90) -o $@ $(OBJ_DIR)/*.o $(COMPFLAGS) 

$(PROGNAME).o: $(PROGNAME).f90 $(OBJS)
	$(F90) $(OFLAGS) -I$(OBJ_DIR) -J$(OBJ_DIR) -c $(PROGNAME).f90 $(COMPFLAGS) $(DEBFLAGS) -o $(OBJ_DIR)/$@    

# Needed modules 

debugger.o: $(SRC_DIR)/debugger.f90
	$(F90) $(OFLAGS) -I$(OBJ_DIR) -J$(OBJ_DIR) -c $(SRC_DIR)/debugger.f90 $(COMPFLAGS) $(DEBFLAGS) -o $(OBJ_DIR)/$@

tensor_types.o: $(SRC_DIR)/tensor_types.f90
	$(F90) $(OFLAGS) -I$(OBJ_DIR) -J$(OBJ_DIR) -c $(SRC_DIR)/tensor_types.f90 $(COMPFLAGS) $(DEBFLAGS) -o $(OBJ_DIR)/$@

mode_n.o: $(SRC_DIR)/mode_n.f90 tensor_types.o mat_utils.o
	$(F90) $(OFLAGS) -I$(OBJ_DIR) -J$(OBJ_DIR) -c $(SRC_DIR)/mode_n.f90 $(COMPFLAGS) $(DEBFLAGS) -o $(OBJ_DIR)/$@

swap_axis.o: $(SRC_DIR)/swap_axis.f90 tensor_types.o
	$(F90) $(OFLAGS) -I$(OBJ_DIR) -J$(OBJ_DIR) -c $(SRC_DIR)/swap_axis.f90 $(COMPFLAGS) $(DEBFLAGS) -o $(OBJ_DIR)/$@

mat_utils.o: $(SRC_DIR)/mat_utils.f90
	$(F90) $(OFLAGS) -I$(OBJ_DIR) -J$(OBJ_DIR) -c $(SRC_DIR)/mat_utils.f90 $(COMPFLAGS) $(DEBFLAGS) -o $(OBJ_DIR)/$@

cpd_utils.o: $(SRC_DIR)/cpd_utils.f90 tensor_types.o mat_utils.o mode_n.o debugger.o
	$(F90) $(OFLAGS) -I$(OBJ_DIR) -J$(OBJ_DIR) -c $(SRC_DIR)/cpd_utils.f90 $(COMPFLAGS) $(DEBFLAGS) -o $(OBJ_DIR)/$@



# For cleaning up the compilation files
clean:
	rm $(OBJ_DIR)/* 
	rm $(PROGNAME).x
