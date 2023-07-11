#!/usr/bin/bash

R -e "rmarkdown::render('input_output_relationships.Rmd', output_format='pdf_document')"
