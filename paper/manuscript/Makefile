LATEX=pdflatex
TEXFLAGS=

TARGET=fastZKP.pdf 

all: $(TARGET)

$(TARGET): $(TARGET:.pdf=.tex) *.tex *.bib Makefile
	$(LATEX) $(TEXFLAGS) $<
	$(LATEX) $(TEXFLAGS) $<
	bibtex $(basename $<)
	bibtex $(basename $<)
	$(LATEX) $(TEXFLAGS) $<
	$(LATEX) $(TEXFLAGS) $<

test: test-fonts test-info

test-fonts: $(TARGET)
	pdffonts $<

test-info: $(TARGET)
	pdfinfo $<

print: $(TARGET)
	pdftops $< - | lpr -P lws510

view: $(TARGET)
	kpdf $<

clean:
	rm -f $(TARGET) *.out *.log *.aux *.bbl *.blg *.bak .*.xtx.cache parsetab.py \#*\# *~