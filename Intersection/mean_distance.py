import csv
import math
from albers_conversion import convert_to_albers

# Step 1: Read the CSV file
file_path = 'CoGo_Bikerental_Colorado_US/cleaned_data/cleaned202007-cogo-tripdata.csv'
data = []

with open(file_path, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        data.append(row)



for row in data:
    row['start_x'], row['start_y'] = convert_to_albers(float(row['start_lng']), float(row['start_lat']))
    row['end_x'], row['end_y'] = convert_to_albers(float(row['end_lng']), float(row['end_lat']))
    #row['start_x'], row['start_y'] = float(row['start_lng']), float(row['start_lat'])
    #row['end_x'], row['end_y'] = float(row['end_lng']), float(row['end_lat'])

# Step 3: Define the bounding box based on the provided coordinates
bbox_lon_lat = (-83.2101797, 39.8086936, -82.7713119, 40.1573082)
min_x, min_y = convert_to_albers(bbox_lon_lat[0], bbox_lon_lat[1])
max_x, max_y = convert_to_albers(bbox_lon_lat[2], bbox_lon_lat[3])
#min_x, min_y = bbox_lon_lat[0], bbox_lon_lat[1]
#max_x, max_y = bbox_lon_lat[2], bbox_lon_lat[3]
bbox = (min_x, min_y, max_x, max_y)

def calculate_distance(x1, y1, x2, y2):
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

distances = [calculate_distance(row['start_x'], row['start_y'], row['end_x'], row['end_y']) for row in data]
mean_distance = sum(distances) / len(distances)