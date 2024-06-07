# Global Moran Index based on GIS Algorithms, Ch9, p197-198, by Ningchuan Xiao, publ. 2016
# revised by Xiao Cui, 05.05.2024
# applicable for point patterns
# Classes and methods for geospatial algorithms based on points with x, y coordinates

# for GEO877 Spatial Algorithms (points part)
# - requires Python 3.6 or later (or replace f-string with older print format)
# - Point distance measures return results in metres (or assume metres)
####### Point #######

from numpy import sqrt, radians, arcsin, sin, cos

# class and methods for a geometric point
# =======================================
from numpy import sqrt

class Point():
    # initialise
    def __init__(self, x=None, y=None):
        self.x = x
        self.y = y
    
    # representation
    def __repr__(self):
        return f'Point(x={self.x:.2f}, y={self.y:.2f})'

        # Test for equality between Points
    def __eq__(self, other): 
        if not isinstance(other, Point):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return self.x == other.x and self.y == other.y

    # We need this method so that the class will behave sensibly in sets and dictionaries
    def __hash__(self):
        return hash((self.x, self.y))
    
    # calculate Euclidean distance between two points
    def distEuclidean(self, other):
        return sqrt((self.x-other.x)**2 + (self.y-other.y)**2)
    
    # calculate Manhattan distance between two points
    def distManhattan(self, other):
        return abs(self.x-other.x) + abs(self.y-other.y)

    # Haversine distance between two points on a sphere - requires lat/lng converted to radians
    def distHaversine(self, other):
        r = 6371000  # Earth's radius in metres (will return result in metres)
        phi1 = radians(self.y) # latitudes
        phi2 = radians(other.y)
        lam1 = radians(self.x) # longitudes
        lam2 = radians(other.x)

        d = 2 * r * arcsin(sqrt(sin((phi2 - phi1)/2)**2 + 
                                      cos(phi1) * cos(phi2) * sin((lam2 - lam1)/2)**2))
        return d   

    
    # Calculate determinant with respect to three points. Note the order matters here - we use it to work out left/ right in the next method
    def __det(self, p1, p2):
        det = (self.x-p1.x)*(p2.y-p1.y)-(p2.x-p1.x)*(self.y-p1.y)       
        return det

    def leftRight(self, p1, p2):
        # based on GIS Algorithms, Ch2, p11-12, by Ningchuan Xiao, publ. 2016
        # -ve: this point is on the left side of a line connecting p1 and p2
        #   0: this point is collinear
        # +ve: this point is on the right side of the line
        side = int(self.__det(p1, p2))
        if side != 0:
            side = side/abs(side)  # will return 0 if collinear, -1 for left, 1 for right
        return side

###### Moran
class PointSF(Point):
    def __init__(self, x=None, y=None, z=None):
        super().__init__(x, y)
        self.z = z
    
    def __repr__(self):
        return f'PointSF(x={self.x:.2f}, y={self.y:.2f}, z={self.z:.2f})'

def generate_weight_list(points, threshold):
    """
    Generate a weight list based on the distance threshold between points.
    
    Parameters:
        points (list): List of Point objects.
        threshold (float): Maximum distance threshold for considering points as neighbors.
    
    Returns:
        list: Weight list containing pairs of indices representing neighbor relationships.
    """
    wlist = []
    n = len(points)
    for i in range(n):
        for j in range(i + 1, n):
            # Euclidean distance
            distance = ((points[i].x - points[j].x) ** 2 + (points[i].y - points[j].y) ** 2) ** 0.5
            if distance <= threshold:
                wlist.append((i, j))
                wlist.append((j, i))  # Ensure symmetry for Moran's I
    return wlist

def moransi2(points, wlist):
    """
    Compute Moran's I measure for a list of points and a weight list.
    
    Parameters:
        points (list): List of PointSF objects containing attribute values.
        wlist (list): Weight list containing pairs of indices representing neighbor relationships.
    
    Returns:
        float: Moran's I measure.
    """
    # Extract attribute values from points
    z = [point.z for point in points]
    
    n = len(z)
    d = 0.0
    var = 0.0
    mean = float(sum(z)) / n
    for i in range(n):
        var += (z[i] - mean) ** 2
    S0 = len(wlist)
    for e in wlist:
        d += (z[e[0]] - mean) * (z[e[1]] - mean)
    I = n * d / S0 / var
    return I

# Given points with x, y coordinates and z attribute values
points = [PointSF(39, -81, 10), PointSF(40, -82, 20), PointSF(41, -83, 15), PointSF(43, -82, 25), PointSF(41, -84, 30)]

# Define weight list (you can generate it using the previously discussed method)
threshold = 5
wlist = generate_weight_list(points, threshold)

# Compute Moran's I measure using moransi2 function and the weight list
I = moransi2(points, wlist)
print("Moran's I Measure:", I)
