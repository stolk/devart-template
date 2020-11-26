run8: tile8
	./tile8 > out.svg
	chromium-browser out.svg

tile8: tile8.c hsv.h prt.h
	clang -g -O3 -mavx tile8.c -o tile8 -lm


