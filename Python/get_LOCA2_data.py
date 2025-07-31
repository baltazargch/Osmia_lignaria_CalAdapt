import climakitae as ck                           # Import the base climakitae package
import climakitaegui as ckg                       # Import the climakitaegui package

selections = ckg.Select()                         # Initialize selections object
selections.show()                                 # Pull up selections GUI to make data settings
data = sel.retrieve()                             # Retrieve the data from the AWS catalog
data = ck.load(data)                              # Read the data into memory
ckg.view(data)                                    # Generate a basic visualization of the data

from climakitae.core.data_interface import get_data, get_data_options

dataopts = get_data_options() 

# Retrieve temperature data for California
data = get_data(
    variable="Maximum air temperature at 2m",
    downscaling_method="Dynamical", 
    resolution="3 km",
    timescale="daily",
    scenario="SSP 3-7.0",
    cached_area="CA"
)

# Data is returned as an xarray Dataset
print(data)
ck.export(data, filename="example1", format="NetCDF") 

dataopts = get_data_options(
    downscaling_method = "Statistical", 
    resolution = "3 km", 
    timescale = 'daily'
)
