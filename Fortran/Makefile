SHELL = /bin/sh

OBJS_f90d1 = main_demo.o fluspect2021_dataSpec.o fluspect2017_dataSpec.o Fluspect-KM.o 

SRC_DIR_f90d1 = 

OBJS_DIR = obj/
EXE_DIR = bin/

EXE = Fluspect-KM
FC = gfortran
LD = gfortran
IDIR = 
CFLAGS = -Wall -g -Wtabs -J$(OBJS_DIR) $(IDIR)
LFLAGS = 
LIBS = 

VPATH = $(SRC_DIR_f90d1):$(OBJS_DIR)
OBJS = $(addprefix $(OBJS_DIR), $(OBJS_f90d1))

all : $(EXE)

$(EXE) : $(OBJS_f90d1)
	@mkdir -p $(EXE_DIR)
	$(LD) -o $(EXE_DIR)$(EXE) $(OBJS) $(LFLAGS) $(LIBS)

$(OBJS_f90d1):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d1)$(@:.o=.f90) -o $(OBJS_DIR)$@

clean :
	rm -f $(OBJS_DIR)*.*
	rm -f $(EXE_DIR)$(EXE)

# Dependencies of files
main_demo.o: main_demo.f90 Fluspect-KM.o
fluspect2021_dataSpec.o: fluspect2021_dataSpec.f90 fluspect2017_dataSpec.o
fluspect2017_dataSpec.o: fluspect2017_dataSpec.f90
Fluspect-KM.o: Fluspect-KM.f90 fluspect2021_dataSpec.o
