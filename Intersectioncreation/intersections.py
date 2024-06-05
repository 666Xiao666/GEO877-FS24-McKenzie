from NewCabMetric import DistanceCalculation as DC, results
import csv
import matplotlib.pyplot as plt
import re
from read_in_boundary import columbus_bd
from hexgrid_boundary import columbus_hexagons

class createIntersections:
    def __init__(self, csv_file, columbus_bd, columbus_hexagons):
        self.csv_file = csv_file
        self.columbus_bd = columbus_bd
        self.columbus_hexagons = columbus_hexagons
        self.lines = []
        self.intersections = []
        self.line_segments = []
        

    def read_lines(self):
        linedata = [entry['line'] for entry in results]
        for line in linedata:
            line_coords_str = line.strip('()').split(') -> (')
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
            x_coords = [coord[0] for coord in line_coords]
            y_coords = [coord[1] for coord in line_coords]

            line_dict = {
                'line_coords': line_coords,
                'x_coords': x_coords,
                'y_coords': y_coords
            }
            

            self.lines.append(line_dict)



    def segment_intersect(self, p1, q1, p2, q2):
        def other_dir(p1, p2, p3):
            return (p3[1] - p1[1]) * (p2[0] - p1[0]) >= (p2[1] - p1[1]) * (p3[0] - p1[0])

        def on_segment(p, q, r):
            return min(p[0], q[0]) <= r[0] <= max(p[0], q[0]) and min(p[1], q[1]) <= r[1] <= max(p[1], q[1])

        p1_q1 = other_dir(p1, q1, p2) != other_dir(p1, q1, q2)
        p2_q2 = other_dir(p2, q2, p1) != other_dir(p2, q2, q1)

        if p1_q1 and p2_q2:
            x1, y1 = p1
            x2, y2 = q1
            x3, y3 = p2
            x4, y4 = q2

            d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
            if d == 0:
                return False

            t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / d
            u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / d

            return 0 <= t <= 1 and 0 <= u <= 1 and on_segment(p1, q1, (x1 + t * (x2 - x1), y1 + t * (y2 - y1))) and on_segment(p2, q2, (x3 + u * (x4 - x3), y3 + u * (y4 - y3)))

        return False

    def calculate_intersections(self):
        self.read_lines()

        for line in self.lines:
            line_coords = line['line_coords']

            for hexagon in self.columbus_hexagons:
                hexagon_coords = hexagon

                for i in range(len(line_coords) - 1):
                    p1 = line_coords[i]
                    q1 = line_coords[i + 1]
                    self.line_segments.append((p1, q1))

                    for j in range(len(hexagon_coords) - 1):
                        p2 = hexagon_coords[j]
                        q2 = hexagon_coords[j + 1]

                        if self.segment_intersect(p1, q1, p2, q2):
                            self.intersections.append((p1, q1, p2, q2))

    def plot_map(self):
        plt.figure(figsize=(10, 10))

        xs, ys = zip(*self.columbus_bd)
        plt.plot(xs, ys, 'b-')
        plt.fill(xs, ys, 'skyblue', alpha=0.5)

        for hexagon in self.columbus_hexagons:
            hexagon_xs = [coord[0] for coord in hexagon]
            hexagon_ys = [coord[1] for coord in hexagon]
            plt.plot(hexagon_xs, hexagon_ys, 'k-')

        for line in self.lines:
            plt.plot(line['x_coords'], line['y_coords'], color='red', linewidth=1)

        plt.title('Columbus Boundary, Hexagons, and Lines in NAD83 / California Albers Projection')
        plt.xlabel('X (meters)')
        plt.ylabel('Y (meters)')

        plt.tight_layout()
        plt.show()

    def get_intersections(self):
        return self.intersections

    def get_lines(self):
        return self.lines

    def get_line_segments(self):
        return self.line_segments

# Usage example
mandist = "/Users/benedikt/Documents/GitHub/GEO877-FS24-McKenzie/Intersection/mandist.csv"

intersect = createIntersections(mandist, columbus_bd, columbus_hexagons)
intersect.calculate_intersections()

intersections = intersect.get_intersections()
lines = intersect.get_lines()
line_segments = intersect.get_line_segments()


intersect.plot_map()