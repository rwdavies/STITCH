PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" STITCH/DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" STITCH/DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: rd check clean

rd:
	cd $(PKGNAME);\
	Rscript -e 'roxygen2::roxygenise(".")'

readme:
	Rscript -e 'rmarkdown::render("README.Rmd", "html_document")'

readme2:
	Rscript -e 'rmarkdown::render("README.Rmd")'
	Rscript -e 'pkgdown::build_site()'
	Rscript -e 'pkgdown::build_articles()'

build: rd
	R CMD build $(PKGSRC)

install: build
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check:
	# Rscript -e 'rcmdcheck::rcmdcheck("$(PKGNAME)_$(PKGVERS).tar.gz", args="--as-cran")'
	Rscript -e 'devtools::check()'

check2: build
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz

clean:
	$(RM) -r $(PKGNAME).Rcheck/

