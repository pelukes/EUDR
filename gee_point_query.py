import ee
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from tkinter import Tk, Button, filedialog, Label
import requests
from PIL import Image
from io import BytesIO
import os
from datetime import datetime  # For handling timestamps

# Initialize Google Earth Engine API
ee.Initialize()

# Function to select a shapefile through a GUI dialog
def select_shapefile():
    shapefile_path = filedialog.askopenfilename(
        title="Select Shapefile with points",
        filetypes=[("Shapefiles", "*.shp")]
    )
    return shapefile_path

# Function to get JRC Global Forest Cover 2020 and Hansen GFC logging data from 2020 onwards
def get_forest_cover_and_hansen_image_url(geometry):
    # Get the JRC Global Forest Cover for 2020
    jrc_forest_cover = ee.ImageCollection('JRC/GFC2020/V1').filterDate('2020-12-31').first();
    
    # Get the Hansen Global Forest Change data (2023 version)
    hansen_forest_change = ee.Image('UMD/hansen/global_forest_change_2023_v1_11')
    
    # Select forest loss from 2020 onwards using the 'lossyear' band
    forest_loss_2020_onwards = hansen_forest_change.select('lossyear').updateMask(hansen_forest_change.select('lossyear').gt(20))
    
    # Buffer the geometry
    buffered_geometry = geometry.buffer(1000).bounds()  # 500 meters buffer
    
    # Visualize JRC Global Forest Cover 2020
    jrc_visualized = jrc_forest_cover.visualize(
        bands=['Map'],  # Assuming 'Map' is the correct band for forest cover
        min=0,
        max=100,
        palette=['006400']  # Green for forest
    )
    
    # Visualize Hansen forest loss from 2020 onwards
    hansen_visualized = forest_loss_2020_onwards.visualize(
        bands=['lossyear'],  # Visualize loss year band
        min=21,
        max=23,
        palette=['0000FF','00FF00', 'FF0000']  # Red for forest loss
    )
    
    # Combine JRC forest cover with Hansen forest loss
    combined_forest_cover = ee.ImageCollection([jrc_visualized, hansen_visualized]).mosaic()
    #combined_forest_cover = ee.ImageCollection([jrc_forest_cover, forest_loss_2020_onwards]).mosaic()
    
    # Define thumbnail parameters for the combined image
    thumbnail_params = {
        'region': buffered_geometry.getInfo(),
        'dimensions': [512,512],  # Image size
        'format': 'png'
    }
    
    # Get the URL for the combined forest cover and forest loss thumbnail image
    url = combined_forest_cover.getThumbURL(thumbnail_params)
    
    return url


# Function to calculate NDVI and extract NDVI values at a point
def calculate_ndvi(geometry, start_date, end_date):
    # Get Sentinel-2 harmonized collection
    s2_collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED') \
                    .filterBounds(geometry) \
                    .filterDate(start_date, end_date) \
                    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 10)) \
                    .map(lambda image: image.normalizedDifference(['B8A', 'B4'])
                         .rename('NDVI')
                         .set('system:time_start', image.get('system:time_start')))  # Keep time metadata
    
    # Extract NDVI values at the given point
    def extract_ndvi(image):
        ndvi_value = image.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=geometry,
            scale=500
        ).get('NDVI')
        
        return ee.Feature(None, {
            'NDVI': ndvi_value,
            'system:time_start': image.get('system:time_start')
        })
    
    # Map the extraction over the collection
    ndvi_features = s2_collection.map(extract_ndvi)
    
    return ndvi_features

# Convert Unix timestamp to human-readable datetime
def convert_timestamp_to_date(timestamp):
    return datetime.utcfromtimestamp(timestamp / 1000)  # Convert from milliseconds to seconds

# Function to create NDVI time series chart and include RGB images + JRC Forest Cover map
# Function to create NDVI time series chart and include RGB images + JRC Forest Cover map
def create_ndvi_chart_with_rgb_and_forest(point_geom, plot_id, pdf_writer):
    geometry = ee.Geometry.Point(point_geom.x, point_geom.y)
    
    # Get NDVI time series for 2019-2024
    ndvi_data = calculate_ndvi(geometry, '2017-01-01', '2024-12-31').getInfo()
    
    # Handle cases where NDVI data or time_start is missing
    times = []
    ndvi_values = []
    
    for feature in ndvi_data['features']:
        time = feature['properties'].get('system:time_start', None)
        ndvi = feature['properties'].get('NDVI', None)
        
        if time is not None and ndvi is not None:
            # Convert the time from Unix to human-readable format
            readable_time = convert_timestamp_to_date(time)
            times.append(readable_time)
            ndvi_values.append(ndvi)
   
    # Get Sentinel-2 RGB image URLs for 2020 and 2024
    reference_rgb_url = get_rgb_image(geometry, '2020-06-01', '2020-08-31')
    current_rgb_url = get_rgb_image(geometry, '2024-06-01', '2024-08-31')
    
    # Get JRC and Hansen combined forest cover and logging data for 2020 onwards
    jrc_hansen_url = get_forest_cover_and_hansen_image_url(geometry)
    
    # Fetch images using requests
    reference_rgb = Image.open(BytesIO(requests.get(reference_rgb_url).content))
    current_rgb = Image.open(BytesIO(requests.get(current_rgb_url).content))
    jrc_hansen = Image.open(BytesIO(requests.get(jrc_hansen_url).content))
    
    # Set figure size to A4 landscape (11.69 x 8.27 inches)
    fig, axs = plt.subplots(1, 4, figsize=(20, 5))
    
    # Ensure all subplots have the same size and square aspect ratio
    #plt.subplots_adjust(left=0.05, right=0.95, top=0.85, bottom=0.15, wspace=0.4)
    
    # NDVI Time Series Plot (set equal aspect ratio)
    if len(times) > 0:
        axs[0].plot(times, ndvi_values, marker='o', linestyle='-')
        axs[0].set_title(f'NDVI Time Series for Plot {plot_id}')
        axs[0].set_xlabel('Time')
        axs[0].set_ylabel('NDVI')
    else:
        axs[0].text(0.5, 0.5, 'No NDVI data available', horizontalalignment='center', verticalalignment='center', transform=axs[0].transAxes)
    
    
    # Reference RGB Image (2020) with equal aspect ratio
    axs[1].imshow(reference_rgb)
    axs[1].set_title('S2 Reference (2020)')
    axs[1].axis('off')
    axs[1].set_aspect('equal', 'box')  # Ensure square shape
    
    # Current RGB Image (2024) with equal aspect ratio
    axs[2].imshow(current_rgb)
    axs[2].set_title('S2 Current (2024)')
    axs[2].axis('off')
    axs[2].set_aspect('equal', 'box')  # Ensure square shape
    
    # JRC Global Forest Cover + Hansen Logging (2020 onwards) with equal aspect ratio
    axs[3].imshow(jrc_hansen)
    axs[3].set_title('JRC Forest Cover 2020 + Hansen Logging (2021-2023)')
    axs[3].axis('off')
    axs[3].set_aspect('equal', 'box')  # Ensure square shape
    
    # Save the figure to the PDF
    pdf_writer.savefig(fig)
    plt.close()

# Function to get RGB image with lowest cloud cover and stretch histogram locally
def get_rgb_image(geometry, start_date, end_date):
    # Get Sentinel-2 harmonized collection and filter
    s2_collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED') \
                    .filterBounds(geometry) \
                    .filterDate(start_date, end_date) \
                    .sort('CLOUDY_PIXEL_PERCENTAGE') \
                    .first()
    
    # Compute the min/max pixel values for B4, B3, B2 bands within the geometry
    region_stats = s2_collection.reduceRegion(
        reducer=ee.Reducer.percentile([2, 98]),  # 2nd and 98th percentiles for stretching
        geometry=geometry,
        scale=10,
        bestEffort=True
    )
    
    # Get the min/max values for each band
    min_b4 = region_stats.get('B4_p2').getInfo()
    max_b4 = region_stats.get('B4_p98').getInfo()
    min_b3 = region_stats.get('B3_p2').getInfo()
    max_b3 = region_stats.get('B3_p98').getInfo()
    min_b2 = region_stats.get('B2_p2').getInfo()
    max_b2 = region_stats.get('B2_p98').getInfo()
    
    # Visualize the image with local stretching
    rgb_image = s2_collection.visualize(
        bands=['B11', 'B8A', 'B4'],
        #min=[min_b4, min_b3, min_b2],
        #max=[max_b4, max_b3, max_b2]
        min=[0,0,0],
        max=[4000,4000,1000]
    )
    
    # Create a buffered region around the point
    buffered_geometry = geometry.buffer(1000).bounds()  # 1000 meters buffer
    
    # Define thumbnail parameters
    thumbnail_params = {
        'region': buffered_geometry.getInfo(),
        'dimensions': 512,
        'format': 'png'
    }
    
    # Get URL for the RGB thumbnail
    url = rgb_image.getThumbURL(thumbnail_params)
    
    return url

# Main function to process all points and save charts with RGB images as PDF
def process_shapefile_to_ndvi_time_series(shapefile_gdf):
    # Use PdfPages to save the PDF in A4 landscape format
    with PdfPages('ndvi_time_series_with_rgb_and_forest_cover.pdf') as pdf_writer:
        for idx, row in shapefile_gdf.iterrows():
            point_geom = row.geometry
            plot_id = row['id']  # Assuming 'id' is a column in the shapefile
            create_ndvi_chart_with_rgb_and_forest(point_geom, plot_id, pdf_writer)

# Function to handle button click
def on_button_click():
    shapefile_path = select_shapefile()
    if shapefile_path:
        gdf = gpd.read_file(shapefile_path)
        process_shapefile_to_ndvi_time_series(gdf)
        label.config(text="Process Completed! PDF Saved.")
    else:
        label.config(text="No shapefile selected.")

# GUI Setup
root = Tk()
root.title("NDVI and RGB Image Time Series Generator")

# Button to trigger shapefile selection
button = Button(root, text="Select Shapefile with Points", command=on_button_click)
button.pack(pady=20)

# Label to show the status of processing
label = Label(root, text="")
label.pack(pady=10)

# Start the GUI
root.mainloop()
