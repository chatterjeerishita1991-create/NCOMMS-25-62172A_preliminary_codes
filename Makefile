.PHONY: expression pan-cancer ovarian all

expression:
	Rscript scripts/01_fetch_expression.R

pan-cancer:
	Rscript scripts/02_pan_cancer_figures.R

ovarian:
	Rscript scripts/03_ovarian_de_gsea_figures.R

all: expression pan-cancer ovarian
