import numpy as np
import numpy.linalg as linalg
from math import tan, acos, pi, sqrt

"""
Documentation for functions is mostly ommitted.
Many just preform common mathematical operations
and the ones that don't have some comments and a
name that specifies their use.
"""

def QuadraticEquation(a, b, c):

    plus_answer = (-1 * b + sqrt(b**2 - 4 * a * c)) / (2 * a)

    negative_answer = (-1 * b - sqrt(b**2 - 4 * a * c)) / (2 * a)
     
    return plus_answer, negative_answer

def CreateVector(p1, p2):
    
    x_component = p2[0] - p1[0]

    y_component = p2[1] - p1[1]
    
    return np.array([x_component, y_component, 0])

def Point2Vec(p):
    return np.array([p[0], p[1], 0])

def UnitNormalVector(vec):

    if abs(vec[1]) > 1e-5: # points are not on the same horizontal line, could have done != 0.0 but float comparison scares me 

        norm_slope = -1 * vec[0] / vec[1]

        norm_line = lambda x: x * norm_slope

        norm_vec = np.array([1, norm_line(1), 0])

    else: # Case of horizontal line, 1/m -> undeffined, vector is just straight line along y-axis

        norm_vec = np.array([0, 1, 0])

    if np.cross(norm_vec, vec)[2] > 0: # make sure it is facing away from the shape
            
        norm_vec *= -1

    return norm_vec / linalg.norm(norm_vec)


def PolygonArea(vertices):
    n = len(vertices)
    area = 0

    for i in range(n):
        j = (i + 1) % n
        area += vertices[i][0] * vertices[j][1]
        area -= vertices[j][0] * vertices[i][1]

    area = abs(area) / 2.0
    return area

def ScaleFactor(poly_contour, expansion_factor):

    """
    Calculates the scale factor that must be applied
    to get the desired area expansion.
    """
    
    num_of_points = len(poly_contour)

    enum_poly_contour = list(enumerate(poly_contour))
    
    enum_poly_contour_with_end = enum_poly_contour + [(num_of_points, poly_contour[0])] # adds the start so there is a full loop

    # Create the list of vectors
    vec_create = lambda ind, point: CreateVector(point, enum_poly_contour_with_end[ind + 1][1])

    vec_list = list(enumerate(map(lambda item: vec_create(*item), enum_poly_contour)))

    vec_list_with_end = vec_list + [(num_of_points, vec_list[0][1])] # Again creates a full loop of vectors

    ## Getting the sum of the lengths of the vectors composing the perimeter of the polygon
    vector_sum = sum(map(lambda item: linalg.norm(item[1]), vec_list))

    ## Get each angle in the polygon
    ## Calulating the interior angle via the dot product
    angle_func = lambda ind, vec: (pi - acos((np.dot(vec, vec_list_with_end[ind + 1][1])) / (linalg.norm(vec) * linalg.norm(vec_list_with_end[ind + 1][1]))),
                                   vec,
                                   vec_list_with_end[ind + 1][1])

    ## Get the angles between each vector
    anlge_vec_list = list(map(lambda item: angle_func(*item), vec_list))

    ## Subtract or add the triangle areas depending on whether turn is ccw or cw
    tan_func = lambda theta, v, u: tan((pi - theta) / 2) if (np.cross(v, u)[2] < 0) else (-1 * tan((pi - theta) / 2))

    tan_angle_sum = sum(map(lambda item: tan_func(*item), anlge_vec_list))

    poly_area = PolygonArea(poly_contour)
    
    desired_area = poly_area * expansion_factor - poly_area # Negative result is ok, it will result in a negative scaling factor meaning inward scale

    return QuadraticEquation(tan_angle_sum, vector_sum, -1 * desired_area)[0] # use the plus result

def slope_yint(v, u):

            m = (v[1] - u[1]) / (v[0] - u[0])

            b = v[1] - m * v[0]

            return (m, b)


def TriangulateNewCoord(p1, p2, p3, scale_factor):

        p1_vec = Point2Vec(p1)

        p2_vec = Point2Vec(p2)

        p3_vec = Point2Vec(p3)

        vec1 = CreateVector(p1, p2)

        vec2 = CreateVector(p2, p3)

        norm_1 = UnitNormalVector(vec1) * scale_factor

        norm_2 = UnitNormalVector(vec2) * scale_factor

        if abs(vec1[0]) < 1e-5: # This is the case when vec1 is vertical

            line2_params = slope_yint(p2_vec + norm_2, p3_vec + norm_2)

            new_x = (p1_vec + norm_1)[0]

            new_y = line2_params[0] * new_x + line2_params[1]

        elif abs(vec2[0]) < 1e-5: # Case when vec2 is vertical

            line1_params = slope_yint(p1_vec + norm_1, p2_vec + norm_1)

            new_x = (p2_vec + norm_2)[0]

            new_y = line1_params[0] * new_x  + line1_params[1]

        else:
            
            line1_params = slope_yint(p1_vec + norm_1, p2_vec + norm_1)

            line2_params = slope_yint(p2_vec + norm_2, p3_vec + norm_2)

            new_x = (line2_params[1] - line1_params[1]) / (line1_params[0] - line2_params[0]) 

            new_y = line2_params[0] * new_x + line2_params[1]

        return new_x, new_y
