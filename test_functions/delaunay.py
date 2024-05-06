# 2024-05-05, v1.0: Add delaunay triangulation and Voronoi
"""
Simple structured Delaunay triangulation in 2D with Bowyer-Watson algorithm.
Based on code from Ayron Catteau. Published at http://github.com/ayron/delaunay
"""
# author:  Xiao Cui
# version: 1.0

# --------------------------
# Geometry operations are mainly from:
# Classes and methods for geospatial algorithms based on points with x, y coordinates

# for GEO877 Spatial Algorithms
# - requires Python 3.6 or later (or replace f-string with older print format)
# - Point distance measures return results in metres (or assume metres)

# release history:
# 2021-03-27, v1.0: updated version of points.py (v1.0) for solution 1 (comparing distance measures)
# author:  Sharon Richardson
# version: 1.0

# 2022-03-25, v2.1: Modified Sharon's original polygon code slightly
# author:  Ross Purves


# 2022-04-02, v2.2: Added Segment class and determinants to Point class
# author:  Ross Purves
# version: 2.2
# --------------------------


from numpy import sqrt, radians, arcsin, sin, cos

# class and methods for a geometric point
# =======================================
from numpy import sqrt


####### Point #######
def sort_points(points):
    """
    Sorts a given list of Point objects and removes duplicates
    """
    unique_points = list(set(points))
    return sorted(unique_points, key=lambda p: (p.x, p.y))

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

####### Segment #######
class Segment():    
    # initialise
    def __init__(self, p0, p1):
        self.start = p0
        self.end = p1
        self.length = p0.distEuclidean(p1)
    
    # representation
    def __repr__(self):
        return f'Segment with start {self.start} and end {self.end}.' 

    # Test for equality between Segments - we treat segments going in opposite directions as equal here
    def __eq__(self, other): 
        if (self.start == other.start or self.start == other.end) and (self.end == other.end or self.end == other.start):
            return True
        else:
            return False
            # We need this method so that the class will behave sensibly in sets and dictionaries
    
    def __hash__(self):
        return hash((start.start.x, self.start.y,self.end.x, self.end.y))    
        
    # determine if intersects with another segment (using Point method leftRight)
    # - should we incorporate testing for identical segments and non-zero lengths?
    def intersects(self, other):        
        # create bounding boxes for each segment, using Bbox class
        self_bbox = Bbox(self)
        other_bbox = Bbox(other)

        # test if bounding boxes overlap (using Bbox method testOverlap)
        bbox_overlap = self_bbox.intersects(other_bbox)
               
        # if the two bboxes do not overlap, there can be no intersection        
        if bbox_overlap == False:
            return False
        
        else: # bboxes overlap, test if lines intersect or are collinear
            
            # find side for each point of each line in relation to the points of the other line
            apq = self.start.leftRight(other.start, other.end)
            bpq = self.end.leftRight(other.start, other.end)
            pab = other.start.leftRight(self.start, self.end)
            qab = other.end.leftRight(self.start, self.end)           
            
            # leftRight sums to 0 for both segments if they intersect (+1/-1) or are collinear (0/0)
            if (apq + bpq == 0 and pab + qab == 0): 
                return True
            else:
                return False      
    
####### Bounding Box #######    

# define bounding box for 2 or more points (initialised as PointGroup, Polygon or Segment)
class Bbox():
    
    # initialise
    def __init__(self, data):
    # using built-in `isinstance` to test what class has been used to initialise the object   

        # for Segment objects 
        if isinstance(data, Segment) == True:
            x = [data.start.x, data.end.x]
            y = [data.start.y, data.end.y]
        
        # for PointGroup objects (including Polygon)
        else:      
            x = [i.x for i in data]   # extract all x coords as a list
            y = [i.y for i in data]   # extract all y coords as a list

        # determine corners, calculate centre and area
        self.ll = Point(min(x), min(y))    # lower-left corner (min x, min y)
        self.ur = Point(max(x), max(y))    # upper-right corner (max x, max y)
        self.ctr = Point((max(x)-min(x))/2, (max(y)-min(y))/2)   # centre of box
        self.area = (abs(max(x)-min(x)))*abs((max(y)-min(y)))    # area of box
           
    # representation
    def __repr__(self):
        return f'Bounding box with lower-left {self.ll} and upper-right {self.ur}' 
    
    def __eq__(self, other): 
        if (self.ll == other.ll and self.ur == other.ur):
            return True
        else:
            return False
            # We need this method so that the class will behave sensibly in sets and dictionaries
    
    def __hash__(self):
        return hash(self.ll, other.ur)  
        
    # test for overlap between two bounding boxes
    def intersects(self, other):       
        # test if bboxes overlap (touching is not enough to be compatible with the approach to segments)
        if (self.ur.x > other.ll.x and other.ur.x > self.ll.x and
            self.ur.y > other.ll.y and other.ur.y > self.ll.y):
            return True
        else:
            return False
        
    # Check if bounding box is contained by this one     
    def contains(self, other):
        if (self.ur.x >= other.ur.x and self.ll.x <= other.ll.x and
            self.ur.y >= other.ur.y and self.ll.y <= other.ll.y):
            return True
        else:
            return False
        
    def intersectsRegion(self, other):
        if self == other:
            return self
        elif self.contains(other):
            return other
        elif self.intersects(other):
            llx = max(self.ll.x, other.ll.x)
            lly = max(self.ll.y, other.ll.y)
            urx = min(self.ur.x, other.ur.x)
            ury = min(self.ur.y, other.ur.y)
            return Bbox([Point(llx, lly), Point(urx, ury)])
        else:
            return None
                        
    # contains includes points on boundaries, otherwise we have problems when points are used to define
    # box but not in it
    def containsPoint(self, p):
        if (self.ur.x >= p.x and p.x >= self.ll.x and
            self.ur.y >= p.y and p.y >= self.ll.y):
            return True
        return False
        

# FOR GROUPS OF POINTS

# data provided should be as an array of points with x, y coordinates.

# class for a group of Points, assumes initial data is unsorted, spatially
class PointGroup(): 
    # initialise
    def __init__(self, data=None, xcol=None, ycol=None):
        self.points = []
        self.size = len(data)
        for d in data:
            self.points.append(Point(d[xcol], d[ycol]))
    
    # representation
    def __repr__(self):
        return f'PointGroup containing {self.size} points' 
 
    # create index of points in group for referencing
    def __getitem__(self, key):
        return self.points[key]
    

# Polygon class for polygons, assumes initial data is in a spatially sorted order
class Polygon(PointGroup):  
    # initialise
    def __init__(self, data=None, xcol=None, ycol=None):
        self.points = []
        self.size = len(data)
        for d in data:
            self.points.append(Point(d[xcol], d[ycol]))
        self.bbox = Bbox(self)
        
    # representation
    def __repr__(self):
        return f'Polygon PointGroup containing {self.size} points' 
  
    # test if polygon is closed: first and last point should be identical
    def isClosed(self):
        start = self.points[0]
        end = self.points[-1]
        return start == end

    def removeDuplicates(self):
        oldn = len(self.points)
        self.points = list(dict.fromkeys(self.points)) # Get rid of the duplicates
        self.points.append(self.points[0]) # Our polygon must have one duplicate - we put it back now
        n = len(self.points)
        print(f'The old polygon had {oldn} points, now we only have {n}.')
        
        # find area and centre of the polygon
    # - based on GIS Algorithms, Ch.2 p9-10, by Ningchuan Xiao, publ. 2016 
    
    def __signedArea(self):  # used for both area and centre calculations - this is a private method (only used within the class)
        a = 0
        xmean = 0
        ymean = 0
        for i in range(0, self.size-1):
            ai = self[i].x * self[i+1].y - self[i+1].x * self[i].y
            a += ai
            xmean += (self[i+1].x + self[i].x) * ai
            ymean += (self[i+1].y + self[i].y) * ai

        a = a/2.0   # signed area of polygon (can be a negative)
    
        return a, xmean, ymean
    
    def area(self):
        a, xmean, ymean = self.__signedArea()
        area = abs(a)   # absolute area of polygon
        
        return area

    def centre(self):
        a, xmean, ymean = self.__signedArea() # note we use the signed area here
        centre = Point(xmean/(6*a), ymean/(6*a)) # centre of polygon 
        return centre

    def containsPoint(self, p):
        if (self.bbox.containsPoint(p) == False):
            return False
        
        # Solution, as discussed in lecture, added here
        ray = Segment(p, Point(self.bbox.ur.x+1, p.y))
        count = 0
        
        for i in range(0, self.size-1):
            start = self[i]
            end = self[i+1]
            if (start.y != end.y): 
                if ((p.x < start.x) and (p.y == start.y)): 
                    count = count + 1
                else:
                    s = Segment(start, end)
                    if s.intersects(ray):
                        count = count + 1
              
        if (count%2 == 0):
            return False           

        return True



class Line:
	def __init__(self, m, b):
		self.m = m
		self.b = b


def midpoint(p1, p2):

	m_p = ((p1[0] + p2[0])/2, (p1[1] + p2[1])/2)

	return m_p


def find_b(p, m):
	return p[1] - m * p[0]


def slope(p1, p2):

	# If the slope is vertical
	if p2[0] == p1[0]:
		return 'undefined'
	return (p2[1] - p1[1])/(p2[0] - p1[0])


def slope_perpendicular_bisector(p1, p2):
		
	# If the line is horizontal
	if p2[0] == p1[0]:
		return 0
	
	if p2[1] == p1[1]:
		return 'undefined'
		
	return -((p2[0] - p1[0])/(p2[1] - p1[1]))


def perpendicular_bisector(p1, p2, xmin, xmax, ymin, ymax):


	m = slope_perpendicular_bisector(p1, p2)
	
	m_p = midpoint(p1, p2)

	if m == 'undefined':
		p_r1 = (m_p[0], ymin-2 )
		p_r2 = (m_p[0], ymax+2 )
			
	elif m == 0:
		p_r1 = (xmin-2, m_p[1])
		p_r2 = (xmax+2, m_p[1])
	
	else:

		b = find_b(m_p, m)

		x_small = xmin - 2
		x_large = xmax + 2

		# Find points that the bisector intersect the rectangle
		# Add more cases
		y_r1 = m*(x_small - 1) + b
		p_r1 = (x_small - 1, y_r1)
		y_r2 = m*(x_large+1) + b
		p_r2 = (x_large + 1, y_r2)

	return p_r1, p_r2
	
def intersection(p1, q1, p2, q2, shift=None):
		

	m1 = slope(p1, q1)
	if shift:
		m1 -= shift

	m2 = slope(p2, q2)

	if m1 != 'undefined':
		b1 = find_b(p1, m1)
	if m2 != 'undefined':
		b2 = find_b(p2, m2)

	# If the line is vertical, we take the x of the point to be b

	if m1 == 'undefined' and m2 == 0:
		x = p1[0]
		y = p2[1]
		
	elif m2 == 'undefined' and m1 == 0:
		x = p2[0]
		y = p1[1]

	elif m1 == 'undefined':
		x = p1[0]
		y = m2*x+b2

	elif m2 == 'undefined':
		x = p2[0]
		y = m1*x+b1
		

	else:
		x = (b2-b1)/(m1-m2)
	
		y = m1*x+b1

	return (x, y)

#### triangle
class Triangle(object):
    def __init__(self, point_1, point_2, point_3):
        self.point_1 = point_1
        self.point_2 = point_2
        self.point_3 = point_3

        self.line_1 = Segment(self.point_1, self.point_2)
        self.line_2 = Segment(self.point_2, self.point_3)
        self.line_3 = Segment(self.point_3, self.point_1)

        self.clockwise = self.is_clockwise()

    def __str__(self):
        return 'Triangle({}, {}, {})'.format(self.point_1,
                                             self.point_2,
                                             self.point_3)

    def __repr__(self):
        return 'Triangle({}, {}, {})'.format(self.point_1,
                                             self.point_2,
                                             self.point_3)

    def __eq__(self, other):
        """
        Determines if the two triangles have the same points even if they are
        in a different order.
        """
        if isinstance(other, Triangle):
            point_1_eq = (self.point_1 == other.point_1 or
                          self.point_1 == other.point_2 or
                          self.point_1 == other.point_3)
            point_2_eq = (self.point_2 == other.point_1 or
                          self.point_2 == other.point_2 or
                          self.point_2 == other.point_3)
            point_3_eq = (self.point_3 == other.point_1 or
                          self.point_3 == other.point_2 or
                          self.point_3 == other.point_3)
            return point_1_eq and point_2_eq and point_3_eq
        return False


    def is_clockwise(self):
        """
        Determines whether the triangle is clockwise or counter clockwise
        """
        orientation = ((self.point_2.y - self.point_1.y) *
                       (self.point_3.x - self.point_2.x) -
                       (self.point_2.x - self.point_1.x) *
                       (self.point_3.y - self.point_2.y))

        if orientation == 0:
            # Points are on the same line
            return None
        if orientation > 0:
            # CW
            return True
        else:
            # CCW
            return False

    def is_collinear(self):
        return self.clockwise is None

    @classmethod
    def get_circumcircle(cls, line_1, line_2):
        """
        The three points of a triangle form a circle. The center of
        the circle can be found via the point where the two perpendicular
        bisectors meet.

        returns :: center, radius of the circle
        """
        midpoint_1 = line_1.get_midpoint()
        perpendicular_line_1 = line_1.get_perpendicular_bisector(midpoint_1)

        midpoint_2 = line_2.get_midpoint()
        perpendicular_line_2 = line_2.get_perpendicular_bisector(midpoint_2)

        center = Line.intersect(perpendicular_line_1,
                                          perpendicular_line_2)
        radius = center.distance(line_1.point_1)
        return center, radius


class ConvexHull(object):
    """
    A polygon that wraps around a given set of points.
    """
    def __init__(self, points):
        # Sort the points first on x, then on y
        sorted_points = sort_points(points)
        self.points = sorted_points
        self.hull_points = []
        self.build()

    def build(self):
        """
        Loops through self.points, building a Triangle ABC. If the triangle is
        Counter Clockwise, point C is more left of point B and is a better
        candidate for being a hull point
        """
        n = len(self.points)
        start_point = 0
        A = 0

        while True:
            self.hull_points.append(self.points[A])
            B = (A + 1) % n  # Loops back around when we get to the last item
                             # in self.points

            for C in range(n):
                clockwise = Triangle(self.points[A], self.points[B], self.points[C]).clockwise

                if clockwise == False:
                    B = C

            # Now B is the point most left with respect to A and is a new
            # hull point.
            A = B

            if A == start_point:
                # We've looped back around and are finished.
                break

    def add_point(self, point):
        """
        Given a point, add it to the list of points and rebuild the hull
        """
        self.points.append(point)
        # The point might already be in the hull
        self.points = sort_points(self.points)

        self.hull_points = []
        self.build()
        
##
def find_difference(old_points, new_points):
    """
    Determines what elements in old_points are missing in new_points
    Note: Does not include new elements from new_points
    """
    return set(old_points) - set(new_points)


class Delaunay(object):
    """
    Creates the triangulation using a sweep line approach.
    """
    def __init__(self, points):
        self.points = sort_points(points)
        self.triangles = []
        self.convex_hull = None
        self.degenerate = False

        self.build_initial()
        self.triangulate()

    def build_initial(self):
        """
        Creates the starter Triangle object(s). Built to handle collinearity in
        the starting points. In such a case, the produced triangles will be
        degenerate.
        """
        # Form a triangle and convex hull from the first 3 points
        starter_hull = [self.points.pop(0) for i in range(3)]
        self.convex_hull = ConvexHull(starter_hull)
        first_triangle = Triangle(*starter_hull)
        if first_triangle.is_collinear():
            self.degenerate = True
            self.get_collinear_triangles()
        else:
            self.triangles.append(Triangle(*starter_hull))

    def get_collinear_triangles(self):
        """
        In the case of the initial points being on the same line, continue
        adding points to the convex hull until we get the first non-collinear
        point.
        """
        while True:
            if not self.points:
                # All the points are collinear
                raise ValueError("All points provided are collinear and a ",
                                 "triangulation could not be formed.")

            new_point = self.points.pop(0)
            new_triangle = Triangle(new_point, *self.convex_hull.hull_points[-2:])
            # We add the point to the convex hull regardless

            if not new_triangle.is_collinear():
                self.triangles = self.build_triangles(new_point, self.convex_hull.hull_points)
                break
            self.convex_hull.add_point(new_point)
        self.convex_hull.add_point(new_point)

    def triangulate(self):
        """
        Builds a list of Triangle objects based on the triangles created in
        build_initial.
        """
        for point in self.points:
            # Copy convex hull point list
            old_convex_hull = copy(self.convex_hull.hull_points)
            # Add new point to convex hull
            self.convex_hull.add_point(point)

            points_to_connect = []

            points_to_connect += self.get_displaced_points(old_convex_hull)

            points_to_connect += self.get_neighboring_points(point)

            # Sort points_to_connect on x, y
            points_to_connect = sort_points(points_to_connect)

            # Create triangles and add to triangle_list
            new_triangles = self.build_triangles(point, points_to_connect)
            self.triangles += new_triangles

    def build_triangles(self, point, points_to_connect):
        """
        Creates a list of triangles using point as a corner in every triangle.
        """
        triangle_list = []
        n = len(points_to_connect)
        for i in range(n-1):
            triangle_list.append(
                Triangle(point, points_to_connect[i], points_to_connect[i+1])
            )
        return triangle_list

    def get_displaced_points(self, old_convex_hull):
        """
        Find the points that are no longer part of the convex hull. These will
        form triangles with the newly added point.
        """
        return find_difference(old_convex_hull, self.convex_hull.hull_points)

    def get_neighboring_points(self, point):
        """
        Find the two hull points next to the newly added point in the convex hull
        """
        # Add the two points next to new_point in convex hull new_point
        new_point_index = self.convex_hull.hull_points.index(point)

        left_point = self.convex_hull.hull_points[new_point_index-1]

        if len(self.convex_hull.hull_points) == new_point_index + 1:
            # The newly added point neighbors the first hull_point
            right_point = self.convex_hull.hull_points[0]
        else:
            right_point = self.convex_hull.hull_points[new_point_index+1]

        return [left_point, right_point]
