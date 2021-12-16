# these examples rely on the R/qtl package, www.rqtl.org
# install R/qtl package, within R, by typing install.packages("qtl")

all: overview.html

R_OPTS=--no-save --no-restore --no-init-file --no-site-file # vanilla, but with --environ

overview.html: Congruence_Classes.Rmd
	R ${R_OPTS} -e "rmarkdown::render('Congruence_Classes.Rmd', output_format = 'bookdown::gitbook', encoding = 'UTF-8')"
	mv *.html ../docs/
	cp -r docs ../docs/
	cp -r figures ../docs/

clean:
	rm -f overview.html   
