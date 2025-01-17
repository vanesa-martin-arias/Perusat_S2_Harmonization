import ee
ee.Initialize()

# Visualization parameters
vis_params = {
    'bands': ['b3', 'b2', 'b1'],  # RGB natural
    'min': 0,
    'max': 3000
}

# Load Sentinel-2 image
s2_ = ee.Image('projects/servir-sco-assets/assets/PeruSat-1_L8_L9_S2_Harmonization/s2img_2019')
s2bands_og = ['B2', 'B3', 'B4', 'B8']
s2bands = ['b1', 'b2', 'b3', 'b4']
s2 = s2_.select(s2bands_og, s2bands)

# Load PeruSat-1 images
vis_params_ps = {
    'bands': ['b3', 'b2', 'b1'],
    'min': 0,
    'max': 3000
}

perusat_2019_1 = ee.Image('projects/ee-perudeforestation/assets/ATCOR_ORTO_DIM_PER1_20190915151221_SEN_MS_000054')
perusat_2019_2 = ee.Image('projects/ee-perudeforestation/assets/ATCOR_ORTO_DIM_PER1_20190915151221_SEN_MS_000372')
psbands_og = ['b1', 'b2', 'b3', 'b4']
ps1 = perusat_2019_1.select(psbands_og)
ps2 = perusat_2019_2.select(psbands_og)

# Create image collection and mosaic
collection = ee.ImageCollection([ps1, ps2])
mosaic_ps1 = collection.mosaic()

# Reproject and combine bands
s2_19 = s2.reproject('EPSG:4326', None, 10)
ps_ = mosaic_ps1.reproject('EPSG:4326', None, 10)
comb_ = ee.Image.cat([s2_19, ps_])
comb = comb_.select(['b1', 'b2', 'b3', 'b4', 'b1_1', 'b2_1', 'b3_1', 'b4_1']) \
           .rename(['b1', 'b2', 'b3', 'b4', 'b1_PS', 'b2_PS', 'b3_PS', 'b4_PS'])

# Chart style
chart_style = {
    'title': 'Frequency (number of pixels with a band value in the bucket)',
    'colors': ['#0000FF', '#89CFF0', '#2E8B57', '#93C572', '#D2042D', '#C04000', '#BF40BF', '#CCCCFF']
}

# Create histogram chart
chart = ui.Chart.image.histogram({
    'image': comb,
    'scale': 10,
    'maxPixels': 1e13
})
chart.setOptions(chart_style)

# Histogram matching functions
bands = ['b1', 'b2', 'b3', 'b4']
image1 = s2_19
image2 = ps_
image_geometry = comb.geometry()

def get_fc(image, band):
    """Get FeatureCollection representation of CDF from histogram"""
    histo = image.reduceRegion(
        reducer=ee.Reducer.histogram({
            'maxBuckets': pow(2, 12),
            'minBucketWidth': 16
        }),
        geometry=image_geometry,
        scale=2,
        maxPixels=1e12,
        tileScale=8
    )
    
    vals_list = ee.List(ee.Dictionary(histo.get(band)).get('bucketMeans'))
    freqs_list = ee.List(ee.Dictionary(histo.get(band)).get('histogram'))
    cdf_array = ee.Array(freqs_list).accum(0)
    total = cdf_array.get([-1])
    normalized_cdf = cdf_array.divide(total)
    array = ee.Array.cat([vals_list, normalized_cdf], 1)
    
    return ee.FeatureCollection(array.toList().map(lambda list_: ee.Feature(None, {
        'dn': ee.List(list_).get(0),
        'probability': ee.List(list_).get(1)
    })))

def equalize(image1, image2, band):
    """Equalize a given band between two images"""
    fc1 = get_fc(image1, band)
    fc2 = get_fc(image2, band)
    
    classifier1 = ee.Classifier.smileRandomForest(100) \
        .setOutputMode('REGRESSION') \
        .train(features=fc1, classProperty='dn', inputProperties=['probability'])
    
    classifier2 = ee.Classifier.smileRandomForest(100) \
        .setOutputMode('REGRESSION') \
        .train(features=fc2, classProperty='probability', inputProperties=['dn'])
    
    b = image2.select(band).rename('dn')
    return b.classify(classifier2, 'probability').classify(classifier1, band)

def match(image1, image2):
    """Map matching function over each band"""
    return ee.Image.cat([equalize(image1, image2, band) for band in bands])

# Perform matching
matched = match(image1, image2)

# Display results
vis_params_display = {'bands': ['b3', 'b2', 'b1'], 'min': 200, 'max': 4000}

Map.addLayer(image1, vis_params_display, 'img-default-s2')
Map.addLayer(image2, vis_params_display, 'to be matched-ps')
Map.addLayer(matched, vis_params_display, 'matched')

# Create final combined image
image_s2 = image1.reproject('EPSG:4326', None, 10)
ps_matched = matched.reproject('EPSG:4326', None, 10)
comb_m = ee.Image.cat([image_s2, ps_matched])
comb_mb = comb_m.select(['b1', 'b2', 'b3', 'b4', 'b1_1', 'b2_1', 'b3_1', 'b4_1']) \
                .rename(['b1', 'b2', 'b3', 'b4', 'b1_PS', 'b2_PS', 'b3_PS', 'b4_PS'])

# Create final histogram chart
chart_style1 = {
    'title': 'Frequency (number of pixels with a band value in the bucket) after matching',
    'colors': ['#0000FF', '#89CFF0', '#2E8B57', '#93C572', '#D2042D', '#C04000', '#BF40BF', '#CCCCFF']
}

chart1 = ui.Chart.image.histogram({
    'image': comb_mb,
    'region': image_geometry,
    'scale': 10,
    'maxPixels': 1e13
})
chart1.setOptions(chart_style1)
