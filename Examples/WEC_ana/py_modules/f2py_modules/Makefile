


# compiler
     FC = f2py
    CXX = cpp

# compile flags
FCFLAGS = -DF2PY_REPORT_ON_ARRAY_COPY=1 -c --fcompiler=intelem --compiler=intelem --f90flags="-g -check all -CA -CB-CS"

# source files and objects
    SRCS = R_tools_fort.F
    RCS = $(SRCS:.F=.f)



# program name
PROGRAM = R_tools_fort


all: $(RCS)
	$(CXX) $(SRCS) $(RCS)
	$(FC) $(FCFLAGS) -m $(PROGRAM) $(RCS)


clean:
	rm -f *.f *.so  *.o *.mod
