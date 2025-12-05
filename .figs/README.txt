1) rulegraph.svg originally obtained with snakemake --rulegraph, the output goes to logs/weavepop.log
You need to remove the lines corresponding to the config that are added by Sanakefile 
 dot -Tsvg weavepop.log > rulegraph.svg
Then edit with Inkscape and saved as SVG, PNG, and PDF

2) graphical_abstract.svg created with Inkscape and saved as SVG, PNG, and PDF

3) pipeline_diagram and cnv_calling_diagram created in Google Slides in magwenelab/Pipeline diagram.pptx
and saved as SVG, PNG, and PDF

4) database created with Duckerd, edited in Inkscape and saved as SVG, PNF and PDF
```
conda activate snakemake
mamba install conda-forge::nodejs
npm install -g duckerd
duckerd -d test/results/02.Dataset/database.db -f pdf -t neutral -o .figs/database.pdf
```
5) shiny.png is a screenshot of the Shiny app using the test dataset

