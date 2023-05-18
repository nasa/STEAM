# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = STEAM
SOURCEDIR     = source
BUILDDIR      = build

.PHONY: all help html doxybuild clean clean_source

TARGETS =  html

# Check for doxygen
HAVE_DOX  := $(shell command -v doxygen  2> /dev/null)

ifdef HAVE_DOX
TARGETS += doxygen
HTML_DEP = doxygen
else
HTML_DEP = 
endif

# Check for libmesh files (really, just check one)
LIBMESH := $(wildcard ../steam/libmesh/init.py)
ifneq ("$(LIBMESH)","")
SPHINXOPTS = -t libmesh
else
SPHINXOPTS = -t nolibmesh
endif

all: $(TARGETS)


html: $(HTML_DEP) clean_source
	#echo $(SPHINXOPTS)
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean: clean_source
	rm -rf build
	rm -rf doxygen/{html,xml,latex,doxygen.log,doxygen.out}

clean_source:
	rm -rf `ls -1d source/*/ | grep -v users_guide`

doxybuild: doxygen
	cd doxygen; doxygen > doxygen.log

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

#
## Catch-all target: route all unknown targets to Sphinx using the new
## "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
#%: Makefile
#	echo "Bar"
#	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
