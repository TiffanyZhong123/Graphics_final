al: al.mdl lex.py main.py matrix.py mdl.py display.py draw.py gmath.py yacc.py
	python main.py skyscraper.mdl
clean:
	rm *pyc *out parsetab.py
	rm ./anim/*.png

clear:
	rm *pyc *out parsetab.py *ppm
