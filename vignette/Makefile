# these examples rely on the R/qtl package, www.rqtl.org
# install R/qtl package, within R, by typing install.packages("qtl")

all: Congruence_Classes.html

R_OPTS=--no-save --no-restore --no-init-file --no-site-file # vanilla, but with --environ

Congruence_Classes.html: Congruence_Classes.Rmd
	R ${R_OPTS} -e "rmarkdown::render('Congruence_Classes.Rmd', output_format = 'bookdown::gitbook', encoding = 'UTF-8')"

Congruence_Classes.pdf: Congruence_Classes.Rmd
	R ${R_OPTS} -e "rmarkdown::render('Congruence_Classes.Rmd', output_format = 'bookdown::pdf_document2', encoding = 'UTF-8')"

clean:
	rm -f Congruence_Classes.tex Congruence_Classes.log Congruence_Classes.out Congruence_Classes.aux Congruence_Classes.toc texput.log Congruence_Classes.bbl Congruence_Classes.blg   
