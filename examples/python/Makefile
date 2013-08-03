NULL = 

EXTRA_DIST = example.py

SOURCE_FILES = \
	mpsolve.pyx \
	mpsolve.pxd \
	$(NULL)

mpsolve.so: $(SOURCE_FILES)
	python setup.py build_ext
	cp build/lib*/mpsolve.so .

clean:
	rm -r mpsolve.so
	rm -r build 
	rm -r *.c
