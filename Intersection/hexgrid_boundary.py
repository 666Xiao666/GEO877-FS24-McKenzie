from read_in_boundary import columbus_bd
from mean_distance import mean_distance
import math

# Step 4: Generate hexagonal grid over the boundary
def calculate_polygons(startx, starty, endx, endy, radius):
    sl = (2 * radius) * math.tan(math.pi / 6)
    p = sl * 0.5
    b = sl * math.cos(math.radians(30))
    w = b * 2
    h = 2 * sl

    startx = startx - w
    starty = starty - h
    endx = endx + w
    endy = endy + h

    origx = startx
    origy = starty

    xoffset = b
    yoffset = 3 * p

    polygons = []
    row = 1

    while starty < endy:
        if row % 2 == 0:
            startx = origx + xoffset
        else:
            startx = origx
        while startx < endx:
            p1x = startx
            p1y = starty + p
            p2x = startx
            p2y = starty + (3 * p)
            p3x = startx + b
            p3y = starty + h
            p4x = startx + w
            p4y = starty + (3 * p)
            p5x = startx + w
            p5y = starty + p
            p6x = startx + b
            p6y = starty
            poly = [
                (p1x, p1y),
                (p2x, p2y),
                (p3x, p3y),
                (p4x, p4y),
                (p5x, p5y),
                (p6x, p6y),
                (p1x, p1y)]
            polygons.append(poly)
            startx += w
        starty += yoffset
        row += 1
    return polygons

# Get the bounding box of the converted boundary in the NAD83 / California Albers projection
min_x = min(x for x, y in columbus_bd)
min_y = min(y for x, y in columbus_bd)
max_x = max(x for x, y in columbus_bd)
max_y = max(y for x, y in columbus_bd)

# Define the radius for the hexagons
hex_radius = mean_distance  # mean distance of trips

# Generate the hexagonal grid using the converted coordinates
hexagons = calculate_polygons(min_x, min_y, max_x, max_y, hex_radius)

# Function to check if a point is inside a polygon
def is_point_in_polygon(x, y, polygon):
    n = len(polygon)
    inside = False
    p1x, p1y = polygon[0]
    for i in range(n + 1):
        p2x, p2y = polygon[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y
    return inside

# Filter hexagons to keep only those that overlap with the boundary
columbus_hexagons = []
for hexagon in hexagons:
    for vertex in hexagon[:-1]:  # Check each vertex except the duplicate last one
        if is_point_in_polygon(vertex[0], vertex[1], columbus_bd):
            columbus_hexagons.append(hexagon)
            break

# Step 6: Store hexagons in a list of dictionaries
columbus_hexagons_data = []

for i, hexagon in enumerate(columbus_hexagons):
    columbus_hexagons_data.append({
        'hexagon_id': i,
        'polygon': hexagon
    })