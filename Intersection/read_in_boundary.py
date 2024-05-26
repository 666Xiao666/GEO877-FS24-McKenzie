import json
from albers_conversion import convert_to_albers
# Step 1: Read the GeoJSON file
geojson_file_path = 'CoGo_Bikerental_Colorado_US/Corporate_Boundary.geojson'

with open(geojson_file_path, 'r') as f:
    geojson_data = json.load(f)

# Step 2: Extract the single polygon coordinates from the GeoJSON data
# GeoJSON file contains a single Polygon
polygon = []

for feature in geojson_data['features']:
    if feature['geometry']['type'] == 'Polygon':
        polygon = feature['geometry']['coordinates'][0]  # Assuming it's a single polygon
        break
    elif feature['geometry']['type'] == 'MultiPolygon':
        polygon = feature['geometry']['coordinates'][0][0]  # First polygon of the first multipolygon
        break



# Step 4: Convert the coordinates of the polygon to NAD83 / California Albers projection (EPSG:3310)
columbus_bd = [convert_to_albers(lon, lat) for lon, lat in polygon]
#columbus_bd = polygon
