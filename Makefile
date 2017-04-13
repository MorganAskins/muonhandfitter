all:
	gcc -shared -Wl,-soname,libmuviewer.so -o libmuviewer.so -fPIC muviewer.cpp
