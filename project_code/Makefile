run: tile
	./tile > out.svg
	chromium-browser out.svg

tile: tile.c hsv.h prt.h
	clang -g -O3 -mavx tile.c -o tile -lm


