cleanbook: clean gitbook

gitbook:
	Rscript --quiet _render.R "bookdown::gitbook" 
book:
	Rscript --quiet _render.R "bookdown::gitbook" 
pdf:
	Rscript --quiet _render.R "bookdown::pdf_book" 
all:
	Rscript --quiet _render.R 

clean:
	rm -rf docs
