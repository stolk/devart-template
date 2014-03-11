# Edge intersections

To determine whether two polygon overlap, we need to:

- See if any vertices of polygon A lies in polygon B.
- See if any of the edges of polygon A intersect any of the edges of polygon B.

To address the latter, we need a edge-vs-edge test.

I've implemented both a scalar version, and a SIMD vector version.
The latter is able to test 8 edges against a single other edge in one go (and no branching.)

It took some time to iron out the bugs, but having the scalar code as a baseline helps.

