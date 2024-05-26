#This is the file containing the clas DistanceCalculation. This class is for the calculation of the Manhattan Distance and returns the following: 
# -start latitude and start longitude
# -end latitude and end longitude
# -distance (Manhattan)
# -the coordinates of the line (should always be 3 coordinates -> one is new the other 2 are start and end coordinate)

import csv
from albers_conversion import convert_to_albers, convert_from_albers
# Quick addition: I left the correction for the slight offset northward alignment. Looking more closely at the city's road system, I could see that the offset could vary a bit in certain parts of the city.
file_path = 'CoGo_Bikerental_Colorado_US/cleaned_data/cleaned202007-cogo-tripdata.csv'
data = []

with open(file_path, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        data.append(row)

class DistanceCalculation:
    def __init__(self, file_path):
        self.file_path = file_path
        self.results = []

    def manhattan_distance_with_least_turns(self, x1, y1, x2, y2):
        #if the start coordinates from start to end match each other it will return 0 (no Manhattan calculation possible)
        if x1 == x2 and y1 == y2:
            return None, None

        # this code line is responsible for the calculation of the Manhattan distance
        distance = abs(x1 - x2) + abs(y1 - y2)

        #The next code line finds the Manhattan route with the least turning points (our assumption)
        #Due to the assumption, the turning point will always be where either x or y coordinates match each other (total route will be like a right-angled triangle)
        #Then it will directely return the line representation (the coordinate in the middle is the connection point between the 2 stations)
        # Generate the line representation minimizing turns
        if x1 == x2 or y1 == y2:
            line_representation = f"({x1}, {y1}) -> ({x2}, {y2})"
        else:
            line_representation = f"({x1}, {y1}) -> ({x1}, {y2}) -> ({x2}, {y2})"

        return distance, line_representation

    def process_data(self):
        with open(self.file_path, mode='r') as file:
            lines = file.readlines()
            # This code determines the headers (lines[0] is always the 1st row) and separates it with a comma
            # Read the header and determine the index of relevant columns
            header = lines[0].strip().split(',')
            #the next 4 lines of the code have all the same function, it is for finding the position of the corresponding index and gives it the name (with _idx end)
            start_lat_idx = header.index('start_lat')
            start_lng_idx = header.index('start_lng')
            end_lat_idx = header.index('end_lat')
            end_lng_idx = header.index('end_lng')
            # Process each line of data
            for line in lines[1:]:
                #here we first define what a row. This splits the line into a list of values
                row = line.strip().split(',')
                #definition of the start and end coordinates with the aid of the before defined row
                start_lat = float(row[start_lat_idx])
                start_lng = float(row[start_lng_idx])
                end_lat = float(row[end_lat_idx])
                end_lng = float(row[end_lng_idx])

                #this line converts the coordinates to the albers projection.
                start_lng, start_lat = convert_to_albers(start_lng, start_lat)
                end_lng, end_lat = convert_to_albers(end_lng, end_lat)


                #this line returns a tuple containing the distance and the line representation (coordinates) by calling the Manhatten distance calculation (line 6-23)
                distance, line_representation = self.manhattan_distance_with_least_turns(start_lng, start_lat, end_lng, end_lat)

                #this part is executed if the start station and the end station is not the same (otherwise they will be None)
                if distance is not None and line_representation is not None:
                    self.results.append({
                        'start_lat': start_lat,
                        'start_lng': start_lng,
                        'end_lat': end_lat,
                        'end_lng': end_lng,
                        'distance': distance,
                        'line': line_representation
                    })


    def get_results(self):
        #this returns all the results we later use in the display_results part 
        return self.results

    def display_results(self):
        #returns all the results for the dataset
        for result in self.results:
            print(result)
    #Saves the results to a csv
    def save_results_to_csv(self, output_file_path):
        with open(output_file_path, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=['start_lng', 'start_lat', 'end_lng', 'end_lat', 'distance', 'line'])
            writer.writeheader()
            #for testing purposes a lower number can be added after self.results (for example self.results[:20])
            for result in self.results:
                writer.writerow(result)



#this is the file path -> change it if you do it on your device
file_path = r'CoGo_Bikerental_Colorado_US/cleaned_data/cleaned202007-cogo-tripdata.csv'

#creating an instance of the class "DistanceCalculation" (do not forget the change the path, otherwise it will not work :) 
distance_calculation = DistanceCalculation(file_path)

#processing the data input
distance_calculation.process_data()

#this is for getting the results properly displayed
results = distance_calculation.get_results()
#for result in results:
#    print(result)

#distance_calculation.save_results_to_csv("mandist.csv")


