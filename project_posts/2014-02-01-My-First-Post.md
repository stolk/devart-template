Insert a description to describe your update to document your creative process. The Commissioned Interactive Artists will be doing the same so you can reference their Project Posts.

![Example Image](../project_images/cover.jpg?raw=true "Example Image")

Technology fascinates me. 
In particular, the evolution in micro processors are interesting.
The latest crop from Intel and AMD feature wide Single-Instruction-Multiple-Data registers through the AVX and AVX2 extensions.

I decided to take the random space filling algorithm as described by David Bourke (http://paulbourke.net/texture_colour/randomtile/), and map it onto the AVX technology inside the latest processors.

As the AVX registers are 256 bit wide, it is ideally suited to hold eight 32-bit floating point values.
This means that for fast tiling, I use eight sided polygons.
