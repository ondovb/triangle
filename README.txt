To compile:

   gcc triangle.cpp -lstdc++ -o triangle

To run:

   ./triangle region.bmp <radius>
   
   Produces:

      out-<radius>-orig.svg (original triangulated quadtree)
      out-<radius>-grow.svg (expanded quadtree)
      out-<radius>-comb.svg (overlay of original and expanded)

To run time trials:

   ./triangle region.bmp

   Produces:
   
      Same as above, with various radii.
      Writes times for 10,000 iterations to STDOUT.

Custom images:

   Images must be in 24-bit bitmap format, and should not be large (~100
   pixels per side). All non-white pixels will be considered black.
