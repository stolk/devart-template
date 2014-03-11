# Building

To build the tile program, it's easiest to use the Makefile on a UNIX machine.
It invokes the clang compiler, but gcc will probably work as well.
You do require AVX support in both compiler and the CPU you are using.
Nearly all recent PCs will have AVX support.

# Running

When running the tile program, redirect the output to a svg file, like so:
	$ ./tile > out.svg

# Viewing

You can view the output with pretty much any browser.
Or use Inkscape (http://inkscape.org) which even lets you edit it.
I use Inkscape to create png files (export as bitmap) by first setting the background alpha to 1.

