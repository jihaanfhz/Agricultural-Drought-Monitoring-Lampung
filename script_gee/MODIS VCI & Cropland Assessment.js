// ============================================================
// PROJECT: AGRICULTURAL DROUGHT MONITORING - LAMPUNG
// PHASE 4: VEGETATION RESPONSE & GAP-FILLED VCI ASSESSMENT
// SCRIPT 07: MODIS VCI WITH GAP-FILLING (CROPLAND & ALL VEGETATION)
// DATASET: MOD13Q1 (NDVI 250m) & MCD12Q1 (Land Cover)
// ============================================================


// ------------------------------------------------------------
// 1. STUDY AREA SETUP
// ------------------------------------------------------------
var lampung = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filter(ee.Filter.eq('ADM0_NAME', 'Indonesia'))
  .filter(ee.Filter.eq('ADM1_NAME', 'Lampung'));

var roi = lampung.geometry();

Map.centerObject(lampung, 8);


// ------------------------------------------------------------
// 2. MODIS LAND COVER (CROPLAND MASK)
// ------------------------------------------------------------
// MCD12Q1 Band LC_Type1 (IGBP): Kelas 12 = Croplands, Kelas 14 = Cropland Mosaics
var landcover = ee.Image('MODIS/061/MCD12Q1/2020_01_01').select('LC_Type1');
var croplandMask = landcover.eq(12).or(landcover.eq(14)).clip(roi);


// ------------------------------------------------------------
// 3. LOAD & FILTER MODIS NDVI (MOD13Q1 250m)
// ------------------------------------------------------------
var startDate = '2015-01-01';
var endDate = '2026-01-01';

var filterQA = function(image) {
  var qa = image.select('SummaryQA');
  var mask = qa.lte(1); // 0 = Good Data, 1 = Marginal Data
  return image.updateMask(mask);
};

var modisNDVI = ee.ImageCollection('MODIS/061/MOD13Q1')
  .filterDate(startDate, endDate)
  .filterBounds(roi)
  .map(filterQA)
  .select('NDVI')
  .map(function(img) {
    var scaled = img.multiply(0.0001).rename('NDVI');
    return scaled.copyProperties(img, ['system:time_start']);
  });


// ------------------------------------------------------------
// 4. CLIMATOLOGICAL MONTHLY MEDIAN (FOR GAP-FILLING BACKUP)
// ------------------------------------------------------------
var months = ee.List.sequence(1, 12);
var years = ee.List.sequence(2015, 2025);

// Median historis per bulan (2015-2025) sebagai penambal gap awan
var monthlyClimatologyNDVI = ee.ImageCollection.fromImages(
  months.map(function(m) {
    m = ee.Number(m);
    var filtered = modisNDVI.filter(ee.Filter.calendarRange(m, m, 'month')).mean();
    return filtered.rename('NDVI_clim').set('month', m);
  })
);


// ------------------------------------------------------------
// 5. MONTHLY COMPOSITES WITH TEMPORAL & SPATIAL GAP-FILLING
// ------------------------------------------------------------
var monthlyNDVIList = years.map(function(year) {
  year = ee.Number(year);
  return months.map(function(month) {
    month = ee.Number(month);
    var start = ee.Date.fromYMD(year, month, 1);
    var end = start.advance(1, 'month');

    // 1. NDVI Aktual Bulan Tersebut
    var monthlyMean = modisNDVI.filterDate(start, end).mean();

    // 2. Backup Gap Filling 1: Median Historis Bulan yang Sama (2015-2025)
    var climBackup = monthlyClimatologyNDVI.filter(ee.Filter.eq('month', month)).first();

    // 3. Gabungkan: Jika piksel berawan/kosong, diisi oleh median historis
    var filledNDVI = monthlyMean.unmask(climBackup);

    // 4. Backup Gap Filling 2: Spatial Focal Smoothing (mengisi sisa piksel kosong di pinggiran)
    var finalNDVI = filledNDVI.unmask(filledNDVI.focal_median(500, 'circle', 'meters')).clip(roi);

    // Composite 1: Seluruh Vegetasi Lampung (Penuh Tanpa Celah Awan)
    var ndviAll = finalNDVI.rename('NDVI_All');

    // Composite 2: Khusus Lahan Pertanian
    var ndviCropland = finalNDVI.updateMask(croplandMask).rename('NDVI_Cropland');

    var yStr = year.format('%04d');
    var mStr = month.format('%02d');
    var bandName = ee.String('VCI_').cat(yStr).cat('_').cat(mStr);

    return ndviAll.addBands(ndviCropland).set({
      'year': year,
      'month': month,
      'band_name': bandName,
      'system:time_start': start.millis()
    });
  });
}).flatten();

var monthlyNDVI = ee.ImageCollection.fromImages(monthlyNDVIList);


// ------------------------------------------------------------
// 6. CALCULATE VCI (ALL VEGETATION & CROPLAND)
// ------------------------------------------------------------
var monthlyVCICollection = monthlyNDVIList.map(function(img) {
  var image = ee.Image(img);
  var m = image.get('month');

  var sameMonthColl = monthlyNDVI.filter(ee.Filter.eq('month', m));

  // A. VCI Seluruh Vegetasi Lampung (Gap-Filled)
  var ndviMinAll = sameMonthColl.select('NDVI_All').reduce(ee.Reducer.min());
  var ndviMaxAll = sameMonthColl.select('NDVI_All').reduce(ee.Reducer.max());
  var vciAll = image.select('NDVI_All').subtract(ndviMinAll)
    .divide(ndviMaxAll.subtract(ndviMinAll).where(ndviMaxAll.subtract(ndviMinAll).eq(0), 0.001))
    .multiply(100).rename('VCI_All').clamp(0, 100);

  // B. VCI Khusus Lahan Pertanian
  var ndviMinCrop = sameMonthColl.select('NDVI_Cropland').reduce(ee.Reducer.min());
  var ndviMaxCrop = sameMonthColl.select('NDVI_Cropland').reduce(ee.Reducer.max());
  var vciCrop = image.select('NDVI_Cropland').subtract(ndviMinCrop)
    .divide(ndviMaxCrop.subtract(ndviMinCrop).where(ndviMaxCrop.subtract(ndviMinCrop).eq(0), 0.001))
    .multiply(100).rename('VCI_Cropland').clamp(0, 100);

  return image.addBands(vciAll).addBands(vciCrop);
});

var vciCollection = ee.ImageCollection.fromImages(monthlyVCICollection);


// ------------------------------------------------------------
// 7. MAP VISUALIZATION (SEPTEMBER 2023 EL NINO SAMPLE)
// ------------------------------------------------------------
var sampleImage = ee.Image(
  vciCollection
    .filter(ee.Filter.eq('year', 2023))
    .filter(ee.Filter.eq('month', 9))
    .first()
);

var vciPalette = [
  '#730000', // VCI < 10: Sangat Kering
  '#e60000', // VCI 10-20: Kering Berat
  '#ffaa00', // VCI 20-35: Kering Sedang
  '#ffff00', // VCI 35-50: Kering Ringan
  '#a800e6', // VCI 50-70: Normal
  '#0070ff'  // VCI > 70: Sangat Sehat
];

// LAYER 1: VCI Seluruh Lampung (Penuh & Mulus / Gap-Filled)
Map.addLayer(
  sampleImage.select('VCI_All'),
  {min: 0, max: 100, palette: vciPalette},
  '1. VCI Seluruh Vegetasi Gap-Filled (Sept 2023)',
  true
);

// LAYER 2: VCI Khusus Lahan Pertanian (Terkikis Masking Sawah)
Map.addLayer(
  sampleImage.select('VCI_Cropland'),
  {min: 0, max: 100, palette: vciPalette},
  '2. VCI Khusus Pertanian (Sept 2023)',
  false
);

// LAYER 3: Batas Administrasi
Map.addLayer(
  lampung.style({color: 'black', fillColor: '00000000', width: 2}),
  {},
  'Batas Admin Lampung'
);


// ------------------------------------------------------------
// 8. TIME-SERIES CHART (COMPARISON)
// ------------------------------------------------------------
var comparisonChart = ui.Chart.image.series({
  imageCollection: vciCollection.select(['VCI_Cropland', 'VCI_All']),
  region: roi,
  reducer: ee.Reducer.mean(),
  scale: 500,
  xProperty: 'system:time_start'
}).setOptions({
  title: 'Perbandingan Tren VCI: Lahan Pertanian vs Seluruh Vegetasi Lampung (2015–2025)',
  hAxis: {title: 'Tahun', format: 'yyyy'},
  vAxis: {title: 'VCI (%)'},
  series: {
    0: {color: '#e60000', label: 'VCI Pertanian'},
    1: {color: '#0070ff', label: 'VCI Seluruh Vegetasi (Gap-Filled)'}
  },
  lineWidth: 1.5,
  pointSize: 1.5
});

print(comparisonChart);


// ------------------------------------------------------------
// 9. EXPORT COMPARATIVE CSV & GEOTIFF STACK TO DRIVE
// ------------------------------------------------------------
var comparativeStatisticsTable = ee.FeatureCollection(
  vciCollection.map(function(image) {
    var y = image.get('year');
    var m = image.get('month');

    var stats = image.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: roi,
      scale: 500,
      maxPixels: 1e13,
      tileScale: 4
    });

    return ee.Feature(null, {
      'year': y,
      'month': m,
      'VCI_All_mean': stats.get('VCI_All'),
      'VCI_Cropland_mean': stats.get('VCI_Cropland'),
      'NDVI_All_mean': stats.get('NDVI_All'),
      'NDVI_Cropland_mean': stats.get('NDVI_Cropland')
    });
  })
);

// A. Export CSV Statistik Pembanding (2015-2025)
Export.table.toDrive({
  collection: comparativeStatisticsTable,
  description: 'MODIS_VCI_Comparative_Statistics_2015_2025',
  folder: 'Agricultural_Drought_Lampung',
  fileNamePrefix: 'MODIS_VCI_Comparative_Statistics_2015_2025',
  fileFormat: 'CSV'
});

// B. Export GeoTIFF Stack VCI Seluruh Vegetasi Gap-Filled (2015-2025)
var vciBandNames = vciCollection.aggregate_array('band_name');
var vciAllStack = vciCollection.select('VCI_All').sort('system:time_start').toBands().rename(vciBandNames);

Export.image.toDrive({
  image: vciAllStack,
  description: 'MODIS_VCI_AllVegetation_Stack_2015_2025',
  folder: 'Agricultural_Drought_Lampung',
  fileNamePrefix: 'MODIS_VCI_AllVegetation_Stack_2015_2025',
  region: roi,
  scale: 250,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

print('========== EXPORT TASKS READY ==========');
print('1. File CSV: Statistik perbandingan VCI Pertanian vs Seluruh Vegetasi (2015-2025)');
print('2. File GeoTIFF Stack: Raster VCI Gap-Filled Seluruh Lampung (2015-2025)');