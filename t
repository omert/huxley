#!/bin/csh -f

latex $1
bibtex $1
latex $1
latex $1
latex $1
latex $1
dvips $1 -o
ps2pdf $1.ps
evince $1.pdf  &
rm $1.aux
rm $1.bbl
rm $1.blg
rm $1.dvi
rm $1.log
rm $1.ps
