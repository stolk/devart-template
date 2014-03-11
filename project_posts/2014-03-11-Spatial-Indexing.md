# Spatial Indexing

When finding overlapping shapes, it is naive to test against all shapes in the plane.
It is better to search the local neighbourhood only.
What is required, is a way of spatially indexing the database of shapes.

To do this, there are some options.
A grid subdivision is easiest, where the bounded plane is divided into grid cells, and for each grid cell we list all occupants (partial or full.)
If we know which cells are occupied by the new shape (to be placed) we can then only test the new shape against the inhabitents of these cells.

I think it is OK to only do this optimization at the latter stages of the algorithm, when shape sizes get really small.
If the bounding sphere of the shape is r, with r much smaller than the cell size, we can then have the cell sizes extend by r, and overlap their neighbours.
This way, we only have to test a single cell when attempting to place a new shape.

A more complex spatial subdivision scheme would be recursive like the quad tree approach.

# More samples

I created another sample image.
This uses a single convex shape, with a fixed hue.
Value and saturation fall off from the centre.
The packing constant c=1.3 was used.

![Example Image](../project_images/sample6.png?raw=true "c=1.3")
