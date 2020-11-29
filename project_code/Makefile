CC=clang-10

CFLAGS=\
	-g \
	-O3 \
	-march=skylake-avx512

runfb: tilefb
	./tilefb > fb.svg
	inkscape --export-png=fb.png --export-width=1280 fb.svg
	display fb.png


run16: tile16
	./tile16 > out.svg
	inkscape --export-png=out.png --export-width=1280 out.svg
	display out.png

run8: tile8
	./tile8 > out.svg
	chromium-browser out.svg

tile8: tile8.o
	$(CC) -o tile8 tile8.o -lm

tile16: tile16.o
	$(CC) -o tile16 tile16.o -lm

tilefb: tilefb.o
	$(CC) -o tilefb tilefb.o -lm


