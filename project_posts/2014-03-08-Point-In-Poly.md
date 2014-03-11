# Point in Poly problem.

To efficiently determine whether a point lies in a polygon, two approaches are possible.
- Counting ray intersections.
- Determining the winding order.

I decided to go with the latter.
I've managed to map this algorithm to an efficient AVX implementation.
It can test eight edges in a single go, without branching.
