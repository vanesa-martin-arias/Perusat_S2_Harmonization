import ee
ee.Initialize()

# Load images
perusat2020 = ee.Image('projects/servir-sco-assets/assets/PeruSat-1_L8_L9_S2_Harmonization/ATCOR_ORTO_DIM_PER1_20200804151918_SEN_MS_001290')
s2img2020 = ee.Image('projects/servir-sco-assets/assets/PeruSat-1_L8_L9_S2_Harmonization/s2img_2020')

# Set map center
Map.setCenter(-74.9484, -8.1118, 12)

# Define bands
s2bands_og = ['B2', 'B3', 'B4', 'B8']
s2bands = ['b1', 'b2', 'b3', 'b4']

# Get geometry and clip
perusat_geometry = perusat2020.geometry()
Map.addLayer(perusat_geometry, {}, 'perusat geometry')
s2_clip = s2img2020.clip(perusat_geometry)
Map.addLayer(s2_clip, {}, 's2 2020 clip', 1)

# Process PeruSat-1
ps_ = perusat2020.reproject('EPSG:4326', None, 10)
psbands_og = ['b1', 'b2', 'b3', 'b4']
ps = ps_.select(psbands_og)

# Combine bands
s2_ = s2_clip.reproject('EPSG:4326', None, 10)
s2_19 = s2_.select(s2bands_og, s2bands)
comb_ = ee.Image.cat([s2_19, ps])
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
    'region': perusat_geometry,
    'scale': 10,
    'maxPixels': 1e13
})
chart.setOptions(chart_style)

# Histogram matching
bands = ['b1', 'b2', 'b3', 'b4']
image1 = s2_19
image2 = ps

def get_fc(image, band):
    """Get FeatureCollection representation of CDF from histogram"""
    histo = image.reduceRegion({
        'reducer': ee.Reducer.histogram({
            'maxBuckets': pow(2, 12),
            'minBucketWidth': 16
        }),
        'geometry': image.geometry(),
        'scale': 2,
        'maxPixels': 1e12,
        'tileScale': 8
    })
    
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
Map.centerObject(image1)
vis_params = {'bands': ['b3', 'b2', 'b1'], 'min': 200, 'max': 4000}
Map.addLayer(image1, vis_params, 'img-default-s2')
Map.addLayer(image2, vis_params, 'to be matched-ps')
Map.addLayer(matched, vis_params, 'matched')

# Create chart after matching
ps_matched = matched.reproject('EPSG:4326', None, 10)
comb_m = ee.Image.cat([s2_19, ps_matched])
comb_mb = comb_m.select(['b1', 'b2', 'b3', 'b4', 'b1_1', 'b2_1', 'b3_1', 'b4_1']) \
                .rename(['b1', 'b2', 'b3', 'b4', 'b1_PS', 'b2_PS', 'b3_PS', 'b4_PS'])

chart_style1 = {
    'title': 'Frequency (number of pixels with a band value in the bucket) after matching',
    'colors': ['#0000FF', '#89CFF0', '#2E8B57', '#93C572', '#D2042D', '#C04000', '#BF40BF', '#CCCCFF']
}

chart1 = ui.Chart.image.histogram({
    'image': comb_mb,
    'region': perusat_geometry,
    'scale': 10,
    'maxPixels': 1e13
})
chart1.setOptions(chart_style1)
