from NewCabMetric import DistanceCalculation as DC
import csv
import matplotlib.pyplot as plt
import re
from read_in_boundary import columbus_bd
from hexgrid_boundary import columbus_hexagons

# Create an instance of DistanceCalculation
dc = DC("/Users/benedikt/Documents/GitHub/GEO877-FS24-McKenzie/CoGo_Bikerental_Colorado_US/cleaned_data/cleaned202007-cogo-tripdata.csv")

# Process the data and save the results to a CSV file. This functions are imported from NewCabMetric.py. They calculate the manhattan distances and safe them to a csv-file.
dc.process_data()
dc.save_results_to_csv("/Users/benedikt/Documents/GitHub/GEO877-FS24-McKenzie/Intersection/mandist.csv")

# Read the generated CSV file
csv_file = "/Users/benedikt/Documents/GitHub/GEO877-FS24-McKenzie/Intersection/mandist.csv"

# Create a figure
#plt.figure(figsize=(10, 10))

# Plot the boundary polygon
#xs, ys = zip(*columbus_bd)
#plt.plot(xs, ys, 'b-')
#plt.fill(xs, ys, 'skyblue', alpha=0.5)

# Plot the hexagons
for hexagon in columbus_hexagons:
    hexagon_xs = [coord[0] for coord in hexagon]
    hexagon_ys = [coord[1] for coord in hexagon]
    #plt.plot(hexagon_xs, hexagon_ys, 'k-')



lines = []
# Read the CSV file and plot each line on the map
with open(csv_file, 'r') as file:
    csv_reader = csv.DictReader(file)
    for line in csv_reader:
        # Extract the line coordinates from the 'line' column
        line_coords_str = line['line'].strip('()').split(') -> (')
        
        line_coords = []
        for coord in line_coords_str:
            coord = coord.strip('() ').replace(') -> (', ' ')
            parts = coord.split(', ')
            if len(parts) != 2:
                print(f"Unexpected number of parts in coordinate: {coord}")
                continue
            x, y = parts
            try:
                line_coords.append([float(x), float(y)])
            except ValueError as e:
                print(f"Error converting coordinates to float: {coord} - {e}")
                continue

        # Extract x and y coordinates separately
        x_coords = [coord[0] for coord in line_coords]
        y_coords = [coord[1] for coord in line_coords]

        # Create a dictionary representing the line
        line_dict = {
            'line_coords': line_coords,
            'x_coords': x_coords,
            'y_coords': y_coords
        }
        
        # Append the line dictionary to the list of lines
        lines.append(line_dict)


        # Plot the line
        #plt.plot(x_coords, y_coords, color='red', linewidth=1)

# Function to check if two line segments intersect
def segment_intersect(p1, q1, p2, q2):
    def other_dir(p1, p2, p3):
        return (p3[1] - p1[1]) * (p2[0] - p1[0]) >= (p2[1] - p1[1]) * (p3[0] - p1[0])

    def on_segment(p, q, r):
        return min(p[0], q[0]) <= r[0] <= max(p[0], q[0]) and min(p[1], q[1]) <= r[1] <= max(p[1], q[1])

    p1_q1 = other_dir(p1, q1, p2) != other_dir(p1, q1, q2)
    p2_q2 = other_dir(p2, q2, p1) != other_dir(p2, q2, q1)

    if p1_q1 and p2_q2:
        x1 = p1[0]
        y1 = p1[1]
        x2 = q1[0]
        y2 = q1[1]
        x3 = p2[0]
        y3 = p2[1]
        x4 = q2[0]
        y4 = q2[1]

        d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
        if d == 0:
            return False

        t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / d
        u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / d

        return 0 <= t <= 1 and 0 <= u <= 1 and on_segment(p1, q1, (x1 + t * (x2 - x1), y1 + t * (y2 - y1))) and on_segment(p2, q2, (x3 + u * (x4 - x3), y3 + u * (y4 - y3)))

    return False

# Find segment intersections between lines and hexagons
intersections = []
line_segments = []

# Read the CSV file and process each line
with open(csv_file, 'r') as file:
    csv_reader = csv.DictReader(file)
    for line in csv_reader:
        # Extract the line coordinates from the 'line' column
        line_coords_str = line['line'].strip('()').split(') -> (')
        line_coords = []
        for coord in line_coords_str:
            coord = coord.strip('() ')
            x, y = coord.split(', ')
            line_coords.append((float(x), float(y)))

        # Check for intersections with each hexagon
        for hexagon in columbus_hexagons:
            hexagon_coords = hexagon

            # Check for intersections between each line segment and hexagon edge
            for i in range(len(line_coords) - 1):
                p1 = line_coords[i]
                q1 = line_coords[i + 1]
                line_segments.append((p1, q1))

                for j in range(len(hexagon_coords) - 1):
                    p2 = hexagon_coords[j]
                    q2 = hexagon_coords[j + 1]

                    if segment_intersect(p1, q1, p2, q2):
                        intersections.append((p1, q1, p2, q2))


# Print the intersections
#print("Intersections:")
#for intersection in intersections:
#    print(intersection)

# Set the title and axis labels
#plt.title('Columbus Boundary, Hexagons, and Lines in NAD83 / California Albers Projection')
#plt.xlabel('X (meters)')
#plt.ylabel('Y (meters)')

# Display the plot
#plt.tight_layout()
#plt.show()