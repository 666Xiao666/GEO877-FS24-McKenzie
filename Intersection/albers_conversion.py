import math

def convert_to_albers(lon, lat):
    # Constants for the Albers equal-area projection (EPSG:3310)
    a = 6378137.0  # Semi-major axis
    e = 0.081819191  # Eccentricity
    phi0 = math.radians(37.66666666666666)  # Latitude of the origin
    phi1 = math.radians(38.73333333333333)  # First standard parallel
    phi2 = math.radians(40.03333333333333)  # Second standard parallel
    lambda0 = math.radians(-82.5)  # Central meridian
    fe = 600000.0  # False easting
    fn = 0.0  # False northing

    phi = math.radians(lat)
    lambda_ = math.radians(lon)
    
    n = (math.sin(phi1) + math.sin(phi2)) / 2
    C = math.cos(phi1)**2 + 2 * n * math.sin(phi1)
    rho0 = a * math.sqrt(C - 2 * n * math.sin(phi0)) / n
    rho = a * math.sqrt(C - 2 * n * math.sin(phi)) / n
    theta = n * (lambda_ - lambda0)
    
    x = fe + rho * math.sin(theta)
    y = fn + rho0 - rho * math.cos(theta)
    
    return x, y



def convert_from_albers(x, y):
    # Constants for the Albers equal-area projection (EPSG:3310)
    a = 6378137.0  # Semi-major axis
    e = 0.081819191  # Eccentricity
    phi0 = math.radians(37.66666666666666)  # Latitude of the origin
    phi1 = math.radians(38.73333333333333)  # First standard parallel
    phi2 = math.radians(40.03333333333333)  # Second standard parallel
    lambda0 = math.radians(-82.5)  # Central meridian
    fe = 600000.0  # False easting
    fn = 0.0  # False northing

    n = (math.sin(phi1) + math.sin(phi2)) / 2
    C = math.cos(phi1)**2 + 2 * n * math.sin(phi1)
    rho0 = a * math.sqrt(C - 2 * n * math.sin(phi0)) / n

    rho = math.sqrt((x - fe)**2 + (rho0 - y + fn)**2)
    theta = math.atan2(x - fe, rho0 - y + fn)

    phi = math.asin((C - (rho * n / a)**2) / (2 * n))
    lambda_ = lambda0 + theta / n

    lat = math.degrees(phi)
    lon = math.degrees(lambda_)

    return lon, lat
