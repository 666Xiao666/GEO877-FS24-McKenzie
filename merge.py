import pandas as pd

# List of file names
files = ['/Users/benedikt/Documents/GitHub/GEO877-FS24-McKenzie/data/cleaned_data/cleaned202007-cogo-tripdata.csv', '/Users/benedikt/Documents/GitHub/GEO877-FS24-McKenzie/data/cleaned_data/cleaned202107-cogo-tripdata.csv', '/Users/benedikt/Documents/GitHub/GEO877-FS24-McKenzie/data/cleaned_data/cleaned202207-cogo-tripdata.csv', '/Users/benedikt/Documents/GitHub/GEO877-FS24-McKenzie/data/cleaned_data/cleaned202307-cogo-tripdata.csv']

# Initialize an empty list to hold dataframes
file_list = []

# Read each file into a dataframe and append to the list
for path in files:
    file = pd.read_csv(path)
    file_list.append(file)

# Concatenate all dataframes in the list
merged_file = pd.concat(file_list)

# Save the concatenated dataframe to a new CSV file
merged_file.to_csv('/Users/benedikt/Documents/GitHub/GEO877-FS24-McKenzie/data/cleaned_data/cleaned_July-cogo-tripdata.csv', index=False)
