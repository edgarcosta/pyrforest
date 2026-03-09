# This Makefile is for convenience as a reminder and shortcut for the most used commands

# Package folder
PACKAGE = pyrforest

# change to your sage command if needed
SAGE = sage

all: install test

build:
	$(SAGE) -python setup.py build_ext

install:
	$(SAGE) -pip install --no-build-isolation -e .

sdist:
	$(SAGE) -python -m build --sdist --no-isolation

test:
	$(SAGE) -t --force-lib $(PACKAGE)

uninstall:
	$(SAGE) -pip uninstall $(PACKAGE)

coverage:
	$(SAGE) -coverage $(PACKAGE)/*

doc:
	cd docs && $(SAGE) -sh -c "make html"

doc-pdf:
	cd docs && $(SAGE) -sh -c "make latexpdf"

clean:
	rm -rf build dist *.egg-info
	rm -rf $(PACKAGE)/*.c

.PHONY: all build install test coverage sdist uninstall clean doc doc-pdf
