
from UtilityFuncs import *

def OutlinePoints(x_coords, y_coords, expansion_percent):

    """
    Given a set of x-coordinates and a set of y-coordinates
    describing a polygon, the algorithm will return a list of 
    x and y coordinates that produce a shape that outlines the 
    original polygon and is larger by the percent specified by
    the expansion percent. The expansion percent can be negative.
    In this case that the new shape will sit inside the original one.

    OutlinePoitns: (listof Float) (listof Float) Float -> (listof (tupleof Float Float))
    """

    expansion_factor = 1 + expansion_percent / 100

    poly_coords = list(zip(x_coords, y_coords))

    scale_factor = ScaleFactor(poly_coords, expansion_factor)

    new_points = []

    for i in range(len(poly_coords)):

        i_next = (i + 1) % len(poly_coords)

        i_prev = (i - 1) % len(poly_coords)
        
        new_points.append(TriangulateNewCoord(poly_coords[i_prev], poly_coords[i], poly_coords[i_next], scale_factor))


    return new_points