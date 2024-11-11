# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -g -I/share/apps/dfprograms/fftw3/include -I/home/lcgarman/gsl/include 
LDFLAGS = -Wl,-rpath,/home/lcgarman/gsl/lib -L/share/apps/dfprograms/fftw3/lib -L/home/lcgarman/gsl/lib
LIBS = -lgsl -lgslcblas -lm -lfftw3 -mcmodel=medium 

# Source files for mjhp2theta function
SRCS_2THETA = mjhp2theta_main.c globals.c structures.c allocate_memory.c cell_parameters.c rwa_functions.c read_binary_abinit.c fftw_functions.c twotheta_functions.c energy_functions.c mjhp2theta_local_potential.c psp_functions.c mjhp2theta_nonlocal_potential.c hkl_functions.c

# Source files for mjhpHKL function
SRCS_HKL = mjhpHKL_main.c globals.c structures.c allocate_memory.c cell_parameters.c rwa_functions.c read_binary_abinit.c fftw_functions.c hkl_functions.c energy_functions.c mjhpHKL_local_potential.c psp_functions.c mjhpHKL_nonlocal_potential.c 

# Source files for prepare2theta function
SRCS_PREP = prep_mjhp2theta_main.c globals.c structures.c allocate_memory.c cell_parameters.c rwa_functions.c read_binary_abinit.c fftw_functions.c recgrid_functions.c twotheta_functions.c hkl_functions.c 

# Object files
OBJS_2THETA = $(SRCS_2THETA:.c=.o)
OBJS_HKL = $(SRCS_HKL:.c=.o)
OBJS_PREP = $(SRCS_PREP:.c=.o)

# Executable name
TARGET_2THETA = mjhp2theta_analyze 
TARGET_HKL = mjhpHKL_analyze 
TARGET_PREP = prepare_mjhp2theta 

all: $(TARGET_2THETA) $(TARGET_HKL) $(TARGET_PREP)

# Linking rule to create the mjhp2theta_analyze 
$(TARGET_2THETA): $(OBJS_2THETA)
	$(CC) $(CFLAGS) $(OBJS_2THETA) $(LDFLAGS) $(LIBS) -o $(TARGET_2THETA)
# Linking rule to create the mjhpHKL_analyze 
$(TARGET_HKL): $(OBJS_HKL)
	$(CC) $(CFLAGS) $(OBJS_HKL) $(LDFLAGS) $(LIBS) -o $(TARGET_HKL)
# Linking rule to create the prepare_mjhp2theta 
$(TARGET_PREP): $(OBJS_PREP)
	$(CC) $(CFLAGS) $(OBJS_PREP) $(LDFLAGS) $(LIBS) -o $(TARGET_PREP)

# Compilation rule for object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS_2THETA) $(OBJS_HKL) $(OBJS_PREP)
