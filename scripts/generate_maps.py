import argparse
import os
import warnings

import fiona
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch
from matplotlib_scalebar.scalebar import ScaleBar
from shapely.geometry import LineString, MultiPoint, Point, Polygon, box

warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')

def calculate_distance(geometry):
    if isinstance(geometry, MultiPoint):
        points = [point for point in geometry.geoms]
        line = LineString(points)
        return round(line.length / 1000, 2)  # Convert meters to kilometers and round to the nearest hundredth
    elif isinstance(geometry, LineString):
        return round(geometry.length / 1000, 2)  # Convert meters to kilometers and round to the nearest hundredth
    elif isinstance(geometry, Point):
        print(f"Warning: Found a Point geometry. Cannot calculate distance.")
        return None
    else:
        print(f"Warning: Unexpected geometry type: {type(geometry)}. Cannot calculate distance.")
        return None

def is_available(row):
    return row['availability'] == 's'

def plot_data_on_basemap(basemap, gdf, institution, filename, output_folder):
    # Define the colors for the categories
    category_colors = {
        'Ocean': '#a3bdd1',  # Correct blue color for ocean
        'Ice shelf': '#cfe1eb',
        'Land': '#f0f0f0',
        'Sub-antarctic_G': 'lightgreen',
        'Sub-antarctic_L': 'lightblue',
        'Ice tongue': 'lightgrey',
        'Rumple': 'yellow',
    }

    fig, ax = plt.subplots(figsize=(10, 10))

    # Set background color to white
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    # Plot basemap with custom colors
    basemap['color'] = basemap['Category'].map(category_colors)

    # Split basemap into inside and outside the circle
    center = Point(0, 0)
    radius = 3e6
    circle_poly = center.buffer(radius)

    inside_circle = basemap[basemap.geometry.intersects(circle_poly)]
    outside_circle = basemap[~basemap.geometry.intersects(circle_poly)]

    # Plot the 'Ocean' category differently inside and outside the circle
    for category, color in category_colors.items():
        inside = inside_circle[inside_circle['Category'] == category]
        outside = outside_circle[outside_circle['Category'] == category]

        # Inside circle: plot with original color
        if not inside.empty:
            if category == 'Ocean':
                inside.plot(ax=ax, color=color, edgecolor='none')
            else:
                inside.plot(ax=ax, color=color, edgecolor='none')

        # Outside circle: plot with white color for 'Ocean'
        if not outside.empty:
            if category == 'Ocean':
                outside.plot(ax=ax, color='white', edgecolor='none')
            else:
                outside.plot(ax=ax, color=color, edgecolor='none')

    # Debug: print basemap information
    print("Basemap categories and colors:")
    print(basemap[['Category', 'color']].drop_duplicates())

    # Define colors and labels based on availability
    availability_colors = {'u': '#fb9a99', 's': '#1f78bc', 'a': 'grey'}
    availability_labels = {'u': 'Unavailable', 's': 'Available', 'a': 'Unsupported'}

    # Plot the data for the institution
    if not gdf.empty:
        for availability, color in availability_colors.items():
            subset = gdf[gdf['availability'] == availability]
            if not subset.empty:
                label = availability_labels.get(availability, 'Other')
                subset.plot(ax=ax, color=color, markersize=0.55, linewidth=0.25, label=label, zorder=2 if availability == 's' else 1)

    # Debug: print GeoDataFrame information
    print("Institution GeoDataFrame:")
    print(gdf.head())

    # Set limits to zoom in on Antarctica
    ax.set_xlim(-3e6, 3e6)
    ax.set_ylim(-3e6, 3e6)

    # Set equal aspect ratio for circular plot and remove axes
    ax.set_aspect('equal')
    ax.axis('off')

    # Add a custom scale bar
    scalebar_length_km = 1000  # Length of scale bar in kilometers
    scalebar_length_m = scalebar_length_km * 1000  # Convert to meters

    # Calculate the position and draw the scale bar
    scalebar_x = -2.5e6
    scalebar_y = -2.8e6
    ax.plot([scalebar_x, scalebar_x + scalebar_length_m], [scalebar_y, scalebar_y], color='black', lw=2)
    ax.text(scalebar_x + scalebar_length_m / 2, scalebar_y - 2e5, f'{scalebar_length_km} km', ha='center', va='top')

    plt.title(f'{institution} Data Availability', fontsize=14)

    # Handle legend
    legend_patches = [Patch(color=color, label=availability_labels[availability]) for availability, color in availability_colors.items()]
    ax.legend(handles=legend_patches, loc='upper right', fontsize=8, title='Availability')

    # Draw a circular boundary
    bbox_size = 3.5e6  # Adjusted to cover entire plotting area
    bbox = box(-bbox_size, -bbox_size, bbox_size, bbox_size)

    # Create the mask polygon with the bounding box as the outer boundary and the circle as a hole
    mask = Polygon(bbox.exterior.coords, [circle_poly.exterior.coords])

    # Convert the mask to a matplotlib patch
    mask_patch = plt.Polygon(list(mask.exterior.coords) + list(mask.interiors[0].coords), closed=True, facecolor='white', edgecolor='none')
    ax.add_patch(mask_patch)

    # Draw the circle outline again
    circle = plt.Circle((0, 0), radius, facecolor='none', edgecolor='black', linewidth=1)
    ax.add_artist(circle)

    # Save the plot to the data folder
    output_path = os.path.join(output_folder, filename)
    plt.savefig(output_path, bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.close(fig)
    print(f"Saved map for {institution} to {output_path}")

def main(input_file_path, map_path, output_folder, statistics_folder):
    print(f"Processing GeoPackage: {input_file_path}")
    print(f"Using map file: {map_path}")
    print(f"Output folder: {output_folder}")
    print(f"Statistics folder: {statistics_folder}")

    # Ensure the output folders exist
    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(statistics_folder, exist_ok=True)

    # Check if input files exist
    if not os.path.exists(input_file_path):
        raise FileNotFoundError(f"Input file not found: {input_file_path}")
    if not os.path.exists(map_path):
        raise FileNotFoundError(f"Map file not found: {map_path}")

    # Read the GeoPackage file into a GeoDataFrame
    fp = gpd.read_file(input_file_path)
    if fp.empty:
        print("The GeoPackage file is empty.")
        return

    # Read the map shapefile
    mp = gpd.read_file(map_path)
    if mp.empty:
        print("The map file is empty.")
        return

    # Check if the 'Category' column exists
    if 'Category' not in mp.columns:
        print("The 'Category' column does not exist in the GeoDataFrame.")
        return

    # Define the colors for the categories
    category_colors = {
        'Ocean': '#a3bdd1',
        'Ice shelf': '#cfe1eb',
        'Land': '#f0f0f0',
        'Sub-antarctic_G': 'lightgreen',
        'Sub-antarctic_L': 'lightblue',
        'Ice tongue': 'lightgrey',
        'Rumple': 'yellow',
    }

    # Create a color column based on the 'Category' column
    mp['color'] = mp['Category'].map(category_colors).fillna('grey')  # Default to grey for undefined categories

    try:
        layers = fiona.listlayers(input_file_path)
        if not layers:
            raise ValueError("No layers found in the GeoPackage.")
        print("Layers in the GeoPackage:")
        print(layers)
    except Exception as e:
        print(f"Error listing layers in the GeoPackage: {e}")
        return

    # Create a dictionary to map institutions to layers
    institution_layers = {}

    # Populate the dictionary
    for layer in layers:
        try:
            gdf = gpd.read_file(input_file_path, layer=layer)  # Read the full layer
            print(f"Successfully read layer: {layer}")
            print(f"Layer shape: {gdf.shape}")
            print(f"Columns: {gdf.columns}")
            if 'institution' in gdf.columns:
                institutions = gdf['institution'].unique()
                for institution in institutions:
                    if institution not in institution_layers:
                        institution_layers[institution] = []
                    institution_layers[institution].append(layer)
        except Exception as e:
            print(f"Error reading layer {layer}: {str(e)}")

    print("Institution layers mapping:")
    print(institution_layers)

    # Create an overview GeoDataFrame
    overview_gdf = gpd.GeoDataFrame()

    # Iterate through each institution and create maps
    for institution, layers in institution_layers.items():
        institution_data = []

        print(f"\nProcessing institution: {institution}")
        print(f"Layers for this institution: {layers}")

        for layer in layers:
            try:
                gdf = gpd.read_file(input_file_path, layer=layer)
                print(f"Successfully read layer '{layer}' for institution '{institution}'")
                print(f"Layer shape: {gdf.shape}")
                institution_subset = gdf[gdf['institution'] == institution]
                print(f"Subset shape for institution: {institution_subset.shape}")
                institution_data.append(institution_subset)
            except Exception as e:
                print(f"Error processing layer {layer} for institution {institution}: {str(e)}")

        # Combine all dataframes for the institution
        if institution_data:
            institution_gdf = gpd.GeoDataFrame(pd.concat(institution_data, ignore_index=True))
        else:
            institution_gdf = gpd.GeoDataFrame()

        # Ensure the CRS matches between the basemap and the GeoDataFrame
        if not institution_gdf.empty:
            institution_gdf = institution_gdf.to_crs(mp.crs)

        print(f"Basemap CRS: {mp.crs}")
        print(f"Institution data CRS: {institution_gdf.crs}")

        # Calculate distances
        institution_gdf['distance'] = institution_gdf['geometry'].apply(calculate_distance)
        invalid_geometries = institution_gdf[institution_gdf['distance'].isnull()]
        if not invalid_geometries.empty:
            print(f"Warning: Found {len(invalid_geometries)} geometries that don't support distance calculation.")
            print(invalid_geometries[['institution', 'campaign']])  # Adjust columns as needed
        institution_gdf = institution_gdf.dropna(subset=['distance'])

        # Add to the overview GeoDataFrame
        if not institution_gdf.empty:
            overview_gdf = pd.concat([overview_gdf, institution_gdf])

        # Plot and save the aggregated data for the institution
        try:
            plot_data_on_basemap(mp, institution_gdf, institution, f'Antarctica_coverage_{institution}.png', output_folder)
        except Exception as e:
            print(f"Error plotting data for {institution}: {str(e)}")

    # Ensure the CRS matches between the basemap and the overview GeoDataFrame
    if not overview_gdf.empty:
        overview_gdf = gpd.GeoDataFrame(overview_gdf).to_crs(mp.crs)

    print(f"Creating overview map")
    print(f"Overview GeoDataFrame shape: {overview_gdf.shape}")
    print(f"Unique institutions in overview: {overview_gdf['institution'].unique()}")

    # Plot and save the overview map
    try:
        plot_data_on_basemap(mp, overview_gdf, "Overview", 'Antarctica_coverage_overview.png', output_folder)
    except Exception as e:
        print(f"Error plotting overview map: {str(e)}")

    print("All maps have been created and saved to the output folder.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate maps and statistics from GeoPackage data.')
    parser.add_argument('input_file_path', type=str, help='Path to the GeoPackage input file.')
    parser.add_argument('map_path', type=str, help='Path to the map shapefile.')
    parser.add_argument('output_folder', type=str, help='Folder to save the generated maps.')
    parser.add_argument('statistics_folder', type=str, help='Folder to save the generated statistics.')

    args = parser.parse_args()
    main(args.input_file_path, args.map_path, args.output_folder, args.statistics_folder)

# Run the script with example input
# python3 generate_maps.py /Users/nathanbekele/Documents/Coding Projects/Research/antarctic_index.gpkg /Users/nathanbekele/Downloads/Quantarctica3/Miscellaneous/SimpleBasemap/ADD_DerivedLowresBasemap.shp scripts/data1 scripts/data_statistics