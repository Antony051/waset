import ee
import geemap.foliumap as geemap
import geopandas as gpd
import pandas as pd
import streamlit as st
from matplotlib import pyplot as plt
import folium.plugins as plugins
import os
import zipfile
import tempfile
import plotly.express as px

# email = "greppo-app@jak-streamlit.iam.gserviceaccount.com"
# key_file = "private-key.json"

# credentials = ee.ServiceAccountCredentials(email=email, key_file=key_file)
# ee.Initialize(credentials)
geemap.ee_initialize()
st.set_page_config(layout="wide")

st.title("Water Spread Area Estimation Tool (WASET)")

st.sidebar.info(
    "This web app was developed by Antony Kishoare J, a Research Scholar under Dr. E. Arunbabu at Centre for Water Resources, Anna University")

st.write("This app is developed to calculate the water spread area of any water body.")

st.info(
    "There are two ways to use this app. If you have the shapefile of the water body, you can directly upload the file and the output will be generated.")
st.info("The second method allows to draw a polygon around any water body and the output will be generated.")

Map = geemap.Map(locate_control=True, plugin_ScaleBar=True,plugin_Draw=False,  plugin_FullScreen=True, center=[13, 80],zoom=10)
plugins.Draw(export=True).add_to(Map)
# taluk = ee.FeatureCollection("projects/jak-streamlit/assets/taluk")
basemaps = list(geemap.basemaps.keys())
basemaps = basemaps[:5]
col_1, col_2 = st.columns(2)

width = 950
height = 600

def save_uploaded_file(file_content, file_name):
    """
    Save the uploaded file to a temporary directory
    """
    import tempfile
    import os
    import uuid

    _, file_extension = os.path.splitext(file_name)
    file_id = str(uuid.uuid4())
    file_path = os.path.join(tempfile.gettempdir(),
                             f"{file_id}{file_extension}")

    with open(file_path, "wb") as file:
        file.write(file_content.getbuffer())

    return file_path

with col_2:
    data = st.file_uploader("Upload a shapefile of the water body", type=["zip", "geojson"])

with col_1:
    option = st.radio("Do you have the shapefile of the water body?", ["Yes", "No"])
    basemap = st.selectbox(
        "Select a basemap",
        basemaps,
        index=basemaps.index("HYBRID")
    )
    Map.add_basemap(basemap)

if option == "No":

    if data is None:
        Map.to_streamlit(width=width, height=height)

    if data is not None:
        uploaded_file = save_uploaded_file(data, data.name)
        input_gdf = gpd.read_file(uploaded_file)
        input_ee = geemap.gdf_to_ee(input_gdf)
        jrc = ee.Image("JRC/GSW1_3/GlobalSurfaceWater").select("max_extent")
        jrc_clip = jrc.clip(input_ee)
        extent = jrc_clip.eq(1).selfMask().rename("extent")

        object_id = extent.connectedComponents(
            connectedness=ee.Kernel.plus(1), maxSize=1023)

        all_tanks = extent.addBands(object_id).reduceToVectors(
            geometry=input_ee,
            crs=jrc.projection(),
            scale=30,
            geometryType='polygon',
            eightConnected=True,
            maxPixels=1e9,
            reducer=ee.Reducer.mean(),
        )


        def addarea(feature):
            return feature.set({'area': feature.geometry().area(1).divide(100 * 100)})


        area = all_tanks.map(addarea)

        all_tanks_gdf = geemap.ee_to_gdf(area)

        threshold = all_tanks_gdf["area"].max() * 0.9

        tank = area.filter(ee.Filter.gte('area', threshold))

        Map.addLayer(tank, {}, "Water Body")
        Map.centerObject(tank, 12)
        Map.to_streamlit(width=width, height=height)

if option == "Yes":

    if data is None:
        Map.to_streamlit(width=width, height=height)

    if data is not None:
        uploaded_file = save_uploaded_file(data, data.name)
        input_gdf = gpd.read_file(uploaded_file)
        tank = geemap.gdf_to_ee(input_gdf)

        Map.addLayer(tank, {}, "Water Body")
        Map.centerObject(tank, 12)
        Map.to_streamlit()

col1, col2 = st.columns(2)

with col1:
    if data is not None:
        uploaded_file = save_uploaded_file(data, data.name)
        gdf = geemap.ee_to_gdf(tank)

        # plot gdf using matplotlib
        fig, ax = plt.subplots(figsize=(10, 10))
        gdf.plot(ax=ax, color='red')
        # plot df box plot
        with st.expander(" Expand to view the water body"):
            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(111)
            gdf.plot(ax=ax, color='#015ce0', alpha=0.5)
            ax.set_title("Water Body")
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")
            st.pyplot(fig)


with col2:
    if data is not None:
        uploaded_file = save_uploaded_file(data, data.name)
        gdf = geemap.ee_to_gdf(tank)
        lat = round(gdf.centroid.iloc[0].y, 2)
        lon = round(gdf.centroid.iloc[0].x, 2)
        db = gpd.read_file("data/tn_taluks.shp")
        # find the taluk intersecting with tank
        st.subheader("Water Body Details")
        geoloc = db.loc[db.intersects(gdf.geometry.iloc[0])]
        try:
            states = geoloc.state.unique()
            districts = geoloc.district.unique()
            taluks = geoloc.taluk.unique()
            st.info("The water body is in the state {}".format(states))
            st.info("It is in {} district".format(districts))
            st.info("The water body lies on {} taluks ".format(taluks))
        except Exception as e:
            st.info("This water body lies outside India")
        tank_area = round(tank.geometry().area(10).divide(100 * 100).getInfo(), 2)
        st.info("The maximum water spread area of the water body is {} Ha".format(tank_area))
        st.info("The water body is located at latitude {} and longitude {}".format(lat, lon))


if data is not None:

    cloudBitmask = ee.Number(2).pow(10).int()
    cirrusBitmask = ee.Number(2).pow(11).int()


    def masks2clouds(image):
        qa = image.select('QA60')
        mask = qa.bitwiseAnd(cloudBitmask).eq(0).And(qa.bitwiseAnd(cirrusBitmask).eq(0))
        return image.updateMask(mask)


    def ndwi(image):

        ndwi = image.normalizedDifference(['B3', 'B8'])
        return image.addBands(ndwi.rename('ndwi'))


    sentinel = ee.ImageCollection('COPERNICUS/S2') \
        .filterBounds(tank) \
        .map(masks2clouds) \
        .map(ndwi)


    def area_fn(shp):
        shp = ee.Feature(shp)
        area = shp.area(1)
        return shp.set("Area", area)

    def wsarea(in_date, out_date, fc):
        img = sentinel.filterDate(in_date, out_date).select("ndwi").max().clip(fc)
        img1 = img.gt(0.1).selfMask()
        r_v = img1.reduceToVectors(
            geometry=fc,
            crs = img1.projection(),
            scale=10,
            eightConnected=False,
            maxPixels=5e12
        )
        try:
            vector = r_v.union(1)
            geometry = geemap.ee_to_gdf(vector)
            a1 = vector.map(area_fn)
            area = a1.aggregate_sum("Area").divide(1e4).getInfo()
            area = round(area, 4)
        except Exception as e:
            geometry = geemap.ee_to_gdf(fc).boundary
            area = 0

        return geometry, area


    year = ['2017', '2018', '2019', '2020', '2021']
    in_month = ['01-01', '02-01', '03-01', '04-01', '05-01', '06-01', '07-01', '08-01', '09-01', '10-01', '11-01',
                '12-01']
    out_month = ['01-31', '02-28', '03-31', '04-30', '05-31', '06-30', '07-31', '08-31', '09-30', '10-31', '11-30',
                 '12-31']

    # create array combining year and month

    in_date = []
    for i in range(len(year)):
        for j in range(len(in_month)):
            in_date.append(year[i] + '-' + in_month[j])

    out_date = []

    for i in range(len(year)):
        for j in range(len(in_month)):
            out_date.append(year[i] + '-' + out_month[j])

    for i in range(len(in_date)):
        in_date[i] = ee.Date(in_date[i])
        out_date[i] = ee.Date(out_date[i])

    month_list = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    index_month = []
    for i in range(len(year)):
        for j in range(len(month_list)):
            index_month.append(month_list[j] + year[i])

    file_name = pd.DataFrame()
    for i in range(len(index_month)):
        file_name = file_name.append([index_month[i] + ".shp"], ignore_index=True)

    feat_temp = gpd.GeoDataFrame(columns=["geometry"])  # to store temp input files
    shape = gpd.GeoDataFrame()  # to store dissolved shapes
    area_df = pd.DataFrame(columns=["Area(Ha)"])  # to store area of shapes
    bounds = gpd.GeoDataFrame(gdf.boundary, columns=["geometry"])

    temp_shp = tempfile.mkdtemp()
    temp_max = tempfile.mkdtemp()
    temp_plt = tempfile.mkdtemp()
    temp_out = tempfile.mkdtemp()

    gdf.to_file(temp_max + "/max_ws.shp")
    with zipfile.ZipFile(temp_max + "/max_ws.zip", 'w') as myzip:
        for file in os.listdir(temp_max):
            if file.endswith("shp"):
                myzip.write(os.path.join(temp_max, file), file)
            if file.endswith("cpg"):
                myzip.write(os.path.join(temp_max, file), file)
            if file.endswith("shx"):
                myzip.write(os.path.join(temp_max, file), file)
            if file.endswith("dbf"):
                myzip.write(os.path.join(temp_max, file), file)

    with st.spinner("Calculating water area"):
        for i in range(len(index_month)):
            # get gdf and area using function
            feat_temp["geometry"], area_df.loc[index_month[i]] = wsarea(in_date[i], out_date[i], tank)
            # set gdf as geometry
            feature = feat_temp.set_geometry(col='geometry')
            # save feature as shapefile
            feature.to_file(os.path.join(temp_shp, file_name.iloc[i, 0]))
            with zipfile.ZipFile(temp_shp + "/" + index_month[i] + ".zip", "w") as zip:
                for file in os.listdir(temp_shp):
                    if file.startswith(index_month[i]):
                        if file.endswith("shp"):
                            zip.write(os.path.join(temp_shp, file), file)
                        if file.endswith("cpg"):
                            zip.write(os.path.join(temp_shp, file), file)
                        if file.endswith("shx"):
                            zip.write(os.path.join(temp_shp, file), file)
                        if file.endswith("dbf"):
                            zip.write(os.path.join(temp_shp, file), file)
            # append boundary to feature
            shape = feat_temp.append(bounds)
            # dissolve to create one feature with boundary
            shp = shape.dissolve()
            # create plot of the feature
            fig, ax = plt.subplots(figsize=(10, 10))
            shp.plot(ax=ax)
            ax.set_title(index_month[i])
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")
            # save the plot
            plt.savefig(temp_plt + "/" + index_month[i] + ".png")
            print("Completed for" + index_month[i])
        area_csv = area_df.to_csv(temp_plt + "/" + "area_csv.csv")
        with zipfile.ZipFile(temp_plt + "/" + "plots.zip", "w") as zip:
            for file in os.listdir(temp_plt):
                if file.endswith(".zip"):
                    zip.write(os.path.join(temp_plt, file), file)

        with zipfile.ZipFile(temp_out +"/" +"output.zip", "w") as zip:
            for file in os.listdir(temp_plt):
                if file.endswith(".png"):
                    zip.write(os.path.join(temp_plt, file), file)
            for file in os.listdir(temp_shp):
                if file.endswith(".zip"):
                    zip.write(os.path.join(temp_shp, file), file)
            for file in os.listdir(temp_max):
                if file.endswith("zip"):
                    zip.write(os.path.join(temp_max, file), file)
            zip.write(temp_plt + "/" + "area_csv.csv", "area_csv.csv")
            zip.write(temp_max + "/" + "max_ws.zip", "max_ws.zip")
    st.success("Water area calculation completed")

    st.subheader("Water Area")
    fig = px.bar(area_df,labels={"x":"Month","y":"Water spread Area (Ha)"}, title="Water spread area in Ha")
    st.plotly_chart(fig, use_container_width=True)

    with open(temp_out +"/" +"output.zip", 'rb') as f:
        st.download_button('Download Zip', f, file_name='output.zip')
    st.stop()
