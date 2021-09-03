# these examples rely on the R/qtl package, www.rqtl.org
# install R/qtl package, within R, by typing install.packages("qtl")

all: Congruence_Classes.pdf

R_OPTS=--no-save --no-restore --no-init-file --no-site-file # vanilla, but with --environ

Congruence_Classes.html: Congruence_Classes.Rnw
	R ${R_OPTS} -e "rmarkdown::render('Congruence_Classes.Rnw')"

Congruence_Classes.pdf: Congruence_Classes.Rnw
	R ${R_OPTS} -e "library(knitr);Sweave2knitr('Congruence_Classes.Rnw');knit('Congruence_Classes-knitr.Rnw')"
	pdflatex Congruence_Classes-knitr.tex
	bibtex Congruence_Classes-knitr
	pdflatex Congruence_Classes-knitr.tex

clean:
	rm -f Congruence_Classes.tex Congruence_Classes.log Congruence_Classes.out Congruence_Classes.aux Congruence_Classes.toc texput.log Congruence_Classes.bbl Congruence_Classes.blg   
