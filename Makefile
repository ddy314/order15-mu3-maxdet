.PHONY: verify paper clean

verify:
	python verify_all.py

paper:
	cd paper && pdflatex -interaction=nonstopmode -halt-on-error main.tex
	cd paper && pdflatex -interaction=nonstopmode -halt-on-error main.tex
	cd paper && pdflatex -interaction=nonstopmode -halt-on-error main.tex

clean:
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.toc paper/*.synctex.gz paper/*.pdf
