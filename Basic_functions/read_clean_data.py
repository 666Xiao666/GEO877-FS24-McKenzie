import os
import csv
from datetime import datetime

class Cogo_Tripdata:
    def __init__(self, input_path, output_path):
        self.input_path = input_path
        self.output_path = output_path
        self.count_original = 0 
        self.count_clean = 0
        self.filelist = os.listdir(input_path)

    # Read the original csv-bike trip data and add each row into a dictionary, wich then are stored within a list.
    def read_csv(self, file_path):
        data = []
        with open(file_path, 'r') as file:
            csv_reader = csv.DictReader(file)
            for row in csv_reader:
                data.append(row)
        return data
    
    # List of rows is converted back into csv with the given attributes as headers.
    def write_csv(self, file_path, data, fieldnames):
        with open(file_path, 'w', newline='') as file:
            csv_writer = csv.DictWriter(file, fieldnames=fieldnames)
            csv_writer.writeheader()
            csv_writer.writerows(data)

    # Starting and ending trip timestamps aare stored as floats. Those are converted into date format to calculate the trip durations.
    # Trip durations are calculated in minutes and added as new attribute 'trip_duration [min]'
    def calculate_trip_duration(self, row):
        start = datetime.strptime(row['started_at'], '%Y-%m-%d %H:%M:%S')
        end = datetime.strptime(row['ended_at'], '%Y-%m-%d %H:%M:%S')  
        trip_duration = int((end - start).total_seconds() / 60) 
        row['trip_duration [min]'] = trip_duration 
        return row

    # The original data consists of corrupt data. To get reliable results:
    # I set a maximum trip duration of 100 min (based on histogram plots of all trip durations).
    # I delete trips that lasts less than 10 minutes AND end at the same coordinates as they start
    # I delete trips that have no coordinate values
    def clean_data(self, data):
        max_duration = 100
        min_duration = 10
        cleaned_data = []

        for row in data:
            if ((row['trip_duration [min]'] < max_duration and
                row['start_lat'] and row['start_lng'] and row['end_lat'] and row['end_lng']) 
                and not
                (row['trip_duration [min]'] < min_duration and row['start_lat'] == row['end_lat'] and row['start_lng'] == row['end_lng'])):
                cleaned_data.append(row)

        return cleaned_data

    # I read in, write, calculate trip distances and clean data for each of the original csv-files. 
    # Then I store the output files in a new folder, which we use for further processing.
    # To have an idea of how many data entries have been removed from the original data sets, I subtract the sum of the final entries from the sum of initial entries. 
    def process_files(self):
        for f in self.filelist:
            source = os.path.join(self.input_path, f)
            data = self.read_csv(source)
            self.count_original += len(data)

            data_with_trip_duration = [self.calculate_trip_duration(row) for row in data]

            cleaned_data = self.clean_data(data_with_trip_duration)
            self.count_clean += len(cleaned_data)

            output_file = os.path.join(self.output_path, "cleaned" + f)
            self.write_csv(output_file, cleaned_data, fieldnames=data[0].keys())

        removed = self.count_original - self.count_clean
        share_removed = round((removed / self.count_original) * 100, 3)

        print(f"From a total number of {self.count_original} entries, {removed} entries ({share_removed} %) have been removed.")

input_path = "CoGo_Bikerental_colorado_US/original_data/"
output_path = "CoGo_Bikerental_colorado_US/cleaned_data/"
cogo_datacleaning = Cogo_Tripdata(input_path, output_path)
cogo_datacleaning.process_files()
