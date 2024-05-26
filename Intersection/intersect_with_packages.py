import pandas as pd
from shapely.geometry import LineString, Polygon
from hexgrid_boundary import calculate_polygons, min_x, min_y, max_x, max_y, mean_distance
import matplotlib.pyplot as plt

# Update this path to the correct location of your CSV file on your system
csv_file_path = 'mandist.csv'

# Load the hexagons data
hexagons = calculate_polygons(min_x, min_y, max_x, max_y, mean_distance)

# Convert hexagons to Shapely Polygon objects
hexagon_polygons = [Polygon(hex) for hex in hexagons]

# Print the number of hexagons loaded
print(f"Number of hexagons loaded: {len(hexagon_polygons)}")

# Load the CSV data
csv_data = pd.read_csv(csv_file_path)

# Parse lines and convert to Shapely LineString objects
def parse_line(line):
    points = line.split('->')
    coordinates = [tuple(map(float, point.strip()[1:-1].split(','))) for point in points]
    return LineString(coordinates)

lines = [parse_line(line) for line in csv_data['line']]

# Print the number of lines loaded
print(f"Number of lines loaded: {len(lines)}")

# Find intersections
intersections = []

for line in lines:
    for hexagon in hexagon_polygons:
        if line.intersects(hexagon):
            intersections.append((line, hexagon))
            # Print debug information for each intersection found
            #print(f"Found intersection: Line {line} intersects with Hexagon {hexagon}")
print(f"Number of lines loaded: {len(intersections)}")
# Check if no intersections are found
if not intersections:
    print("No intersections found.")

# Function to plot hexagons and lines
def plot_hexagons_and_lines(hexagons, lines, intersections):
    fig, ax = plt.subplots()
    
    # Plot hexagons
    for hexagon in hexagons:
        x, y = hexagon.exterior.xy
        ax.plot(x, y, 'blue')
    
    # Plot lines
    for line in lines:
        x, y = line.xy
        ax.plot(x, y, 'green')  # Swap x and y
    
    # Highlight intersections
    for line, hexagon in intersections:
        x, y = line.xy
        ax.plot(x, y, 'red', linewidth=2)  # Swap x and y
    
    ax.set_aspect('equal')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Hexagons and Lines with Intersections Highlighted')
    plt.show()

# Plot the hexagons and lines
plot_hexagons_and_lines(hexagon_polygons, lines, intersections)
