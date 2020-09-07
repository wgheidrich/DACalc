.PHONY: all install 

PIP=pip-3.7
PYTHON=python3.7
SITE-PACKAGES = $(shell $(PIP) show notebook | grep Location | cut -d ' ' -f 2)
CODEMIRROR = $(SITE-PACKAGES)/notebook/static/components/codemirror
CODEMIRROR-DACALC = $(CODEMIRROR)/mode/dacalc
$(info CODEMIRROR-DACALC: $(CODEMIRROR-DACALC))

all: install


# install package, kernel and codemirror mode
install: pip-install kernel-install codemirror-install


pip-install:
	$(PIP) install .

# run after the module is installed
kernel-install: install
	$(PYTHON) -m dacalc.install

# may have to be run as root
codemirror-install: dacalc.js
	mkdir -p $(CODEMIRROR-DACALC)
	cp dacalc.js $(CODEMIRROR-DACALC)

