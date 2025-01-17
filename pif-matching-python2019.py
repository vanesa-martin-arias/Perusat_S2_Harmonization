import ee
ee.Initialize()

# Define geometry
perusat_geometry = ee.Geometry.Polygon([
    [-75.02202399319583, -8.149775716074409],
    [-74.89259131497317, -8.176963230858934],
    [-74.86512549466067, -8.051885407633481],
    [-74.99490149563724, -8.024689443490146]
])

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

# Select bands and create mosaic
ps1 = perusat_2019_1.select(psbands_og)
ps2 = perusat_2019_2.select(psbands_og)
collection = ee.ImageCollection([ps1, ps2])
mosaic_ps1 = collection.mosaic()

# Reproject images
s2_1 = s2.reproject('EPSG:4326', None, 10)
ps_ = mosaic_ps1.reproject('EPSG:4326', None, 10)

# PIF Matching
bands = ee.List(['b1', 'b2', 'b3', 'b4'])
s2_pif = s2_1
ps_pif = ps_.register(s2_pif, 100)  # max displacement in meters

# Calculate spectral distance
distance = s2_pif.spectralDistance(ps_pif, 'SID')
geom = ps_pif.geometry()

# Calculate threshold and create PIF mask
threshold = distance.reduceRegion(
    reducer=ee.Reducer.percentile([5]),
    geometry=geom,
    scale=10,
    bestEffort=True,
    maxPixels=1e13
).getNumber('distance')

pif = distance.lt(threshold)

def match_band(band):
    """Match bands using PIF"""
    before_pif = s2_pif.select([band]).updateMask(pif)
    after_pif = ps_pif.select([band]).updateMask(pif)
    
    args = {
        'reducer': ee.Reducer.linearFit(),
        'geometry': geom,
        'scale': 10,
        'maxPixels': 1e13,
        'bestEffort': True
    }
    
    coeffs = ee.Image.cat([after_pif, before_pif]).reduceRegion(args)
    
    return ps_pif.select([band]) \
        .multiply(coeffs.getNumber('scale')) \
        .add(coeffs.getNumber('offset')) \
        .toUint16()

# Apply matching to all bands
matched_bands = bands.map(match_band)
matched = ee.ImageCollection(matched_bands).toBands().rename(bands)

# Visualization parameters
viz = {'min': 350, 'max': 1300, 'bands': ['b3', 'b2', 'b1']}

# Add layers to map
Map.addLayer(s2_pif, viz, 's2 Before', True)
Map.addLayer(matched, viz, 'ps After (Matched)', True)
Map.addLayer(ps_pif, viz, 'ps After (Original)', False)

# Create charts for before and after comparison
comb_ = ee.Image.cat([s2_1, ps_])
comb = comb_.select(
    ['b1', 'b2', 'b3', 'b4', 'b1_1', 'b2_1', 'b3_1', 'b4_1'],
    ['b1', 'b2', 'b3', 'b4', 'b1_PS', 'b2_PS', 'b3_PS', 'b4_PS']
)

chart_style1 = {
    'title': 'Frequency (number of pixels with a band value in the bucket) Before',
    'colors': ['#0000FF', '#89CFF0', '#2E8B57', '#93C572', '#D2042D', '#C04000', '#BF40BF', '#CCCCFF']
}

chart_before = ui.Chart.image.histogram({
    'image': comb,
    'region': geom,
    'scale': 10,
    'maxPixels': 1e13
})
chart_before.setOptions(chart_style1)

# Create chart for after matching
comb_after = ee.Image.cat([s2_pif, matched])
comb = comb_after.select(
    ['b1', 'b2', 'b3', 'b4', 'b1_1', 'b2_1', 'b3_1', 'b4_1'],
    ['b1', 'b2', 'b3', 'b4', 'b1_PS', 'b2_PS', 'b3_PS', 'b4_PS']
)

chart_style2 = {
    'title': 'Frequency (number of pixels with a band value in the bucket) After',
    'colors': ['#0000FF', '#89CFF0', '#2E8B57', '#93C572', '#D2042D', '#C04000', '#BF40BF', '#CCCCFF']
}

chart_matched = ui.Chart.image.histogram({
    'image': comb_after,
    'region': geom,
    'scale': 10,
    'maxPixels': 1e13
})
chart_matched.setOptions(chart_style2)
