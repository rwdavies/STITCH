PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" STITCH/DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" STITCH/DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: rd check2 clean

rd:
	cd $(PKGNAME);\
	Rscript -e 'roxygen2::roxygenise(".")'

build: rd
	R CMD build $(PKGSRC)

install: build
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check:
	Rscript -e 'devtools::check("STITCH")'

check2: build
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz

clean:
	$(RM) -r $(PKGNAME).Rcheck/

