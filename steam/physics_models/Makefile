BASENAMES = modified_newtonian locate_stagnation_point streamlines

# If they depend on other libraries
modified_newtonian_dep      = ../libmesh/mesh.o
locate_stagnation_point_dep = ../libmesh/mesh.o
streamlines_dep 				 = ../libmesh/mesh.o locate_stagnation_point.o

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

# Make _wrap.cxx from these things.  Also,  strip out the docstrings.
%_wrap.cxx: %.C %.i %.h
	@echo "swigging $@..."
	../../utils/strip_docs.pl $*
	swig -c++ -python $*.i


### Ensure that we get the right python "include" directory
PYTHON_INCLUDE = `python -c "from sysconfig import get_paths;print(get_paths()['include'])"`
%_wrap.o: %_wrap.cxx %.o
	@echo "compiling $@..."
	@${CXX} ${CXXFLAGS} -fPIC -c $< -I${PYTHON_INCLUDE}

# I need to use SECONDEXPANSION here so that I can pull in the dependancies from above
.SECONDEXPANSION:
_%.so: %_wrap.o %.o $($*_dep)
	@echo "linking $@..."
	@${CXX} ${CXXFLAGS} -shared $^ $($*_dep) -o $@ `${LIBMESH_ROOT}/bin/libmesh-config --ldflags`


PY_OUT = $(addsuffix .py,$(BASENAMES))
clean:
	rm -f *.cxx
	rm -f *.o
	rm -f *.so
	rm -f *.pyc $(PY_OUT)
	rm -fr __pycache__
	rm -f *~
