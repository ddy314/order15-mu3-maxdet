.PHONY: verify paper arxiv clean

verify:
	python verify_all.py

paper:
	cd paper && pdflatex -interaction=nonstopmode -halt-on-error main.tex
	cd paper && pdflatex -interaction=nonstopmode -halt-on-error main.tex
	cd paper && pdflatex -interaction=nonstopmode -halt-on-error main.tex

arxiv: paper
	python3 tools/package_arxiv.py

clean:
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.toc paper/*.synctex.gz paper/*.pdf
