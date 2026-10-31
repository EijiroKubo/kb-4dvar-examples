FC      = gfortran
RM      = rm -f

ifeq ($(MAKECMDGOALS),debug)
  FFLAGS = -O0 -g -Wall -fbounds-check -fPIC -cpp -fno-range-check
else
  FFLAGS = -O3 -fPIC -cpp -fno-range-check
endif

LDFLAGS = -shared

TARGET  = calc_forward_ns.so

OBJS = \
	mod_variables.o \
	sub_initial.o \
	sub_utils.o \
	sub_bc_outer.o \
	sub_bc_wall.o \
	sub_rhs3rd_fwd.o \
	sub_hsmac_fwd.o \
	calc_forward_ns.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.F90
	$(FC) $(FFLAGS) -c $< -o $@

debug: $(TARGET)

clean:
	$(RM) *.o *.mod *.so

realclean: clean
