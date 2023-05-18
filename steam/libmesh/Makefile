BASENAMES = mesh init equation_systems system mpi

# We want to make  _BASENAME.so files
default: $(addprefix _,$(addsuffix .so,$(BASENAMES)))

# Keep these intermediate files around
.SECONDARY: $(addsuffix .o,$(BASENAMES)) $(addsuffix _wrap.o,$(BASENAMES)) $(addsuffix _wrap.cxx,$(BASENAMES))

### Define variables to improve readability
CXX = `${LIBMESH_ROOT}/bin/libmesh-config --cxx`
CXXFLAGS := ${CXXFLAGS}
CXXFLAGS += `${LIBMESH_ROOT}/bin/libmesh-config --cppflags --cxxflags --include --libs`

# Make the .o from .C and .h
%.o : %.C %.h
	@echo "compiling $@..."
	@${CXX} ${CXXFLAGS} -fPIC -c $<

mpi4py: get_mpi_dir.py
	@echo "Linking mpi4py directory"
	@ln -sf `./get_mpi_dir.py` .

# Make _wrap.cxx from these things.  Also,  strip out the docstrings.
%_wrap.cxx: %.C %.i %.h mpi4py
	@echo "swigging $@..."
	@../../utils/strip_docs.pl $*
	@swig -c++ -python $*.i

### Ensure that we get the right python "include" directory
PYTHON_INCLUDE = `python -c "from sysconfig import get_paths;print(get_paths()['include'])"`
%_wrap.o: %_wrap.cxx %.o
	@echo "compiling $@..."
	@${CXX} ${CXXFLAGS} -fPIC -c $< -I${PYTHON_INCLUDE}

_%.so: %_wrap.o %.o
	@echo "linking $@..."
	@${CXX} ${CXXFLAGS} -shared $*.o $< -o $@ `${LIBMESH_ROOT}/bin/libmesh-config --ldflags`


PY_OUT = $(addsuffix .py,$(BASENAMES))
clean:
	rm -f *.cxx
	rm -f *.o
	rm -f *.so
	rm -f *.pyc $(PY_OUT)
	rm -fr __pycache__
	rm -f *~
	rm -f mpi4py
