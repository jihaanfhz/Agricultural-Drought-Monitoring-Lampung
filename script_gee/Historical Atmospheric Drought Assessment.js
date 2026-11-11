// PHASE 2: HISTORICAL ATMOSPHERIC DROUGHT ASSESSMENT
// SCRIPT 02: CHIRPS MONTHLY RAINFALL & DEFICIT (1981-2025)
// ============================================================



//  STUDY AREA
// ------------------------------------------------------------

var admin = ee.FeatureCollection('FAO/GAUL/2015/level1');

var lampung = admin
  .filter(ee.Filter.eq('ADM0_NAME', 'Indonesia'))
  .filter(ee.Filter.eq('ADM1_NAME', 'Lampung'));

var roi = lampung.geometry();

print('========== STUDY AREA ==========');
print('Boundary feature count:', lampung.size());
print('Administrative name:', lampung.aggregate_array('ADM1_NAME'));



//  LOAD CHIRPS DAILY (1981-2025)
// ------------------------------------------------------------

var startDate = '1981-01-01';
var endDate = '2026-01-01';

var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .filterDate(startDate, endDate)
  .filterBounds(roi)
  .select('precipitation');

print('========== CHIRPS HISTORICAL DATA ==========');
print('Expected daily images:', 16436);
print('Actual daily images:', chirps.size());
print('First image date:', ee.Date(chirps.first().get('system:time_start')));
print('Last image date:', ee.Date(chirps.sort('system:time_start', false).first().get('system:time_start')));



//  CREATE MONTHLY RAINFALL ACCUMULATION (540 MONTHS)
// ------------------------------------------------------------

var years = ee.List.sequence(1981, 2025);
var months = ee.List.sequence(1, 12);

var monthlyRainfall = ee.ImageCollection.fromImages(
  years.map(function(year) {
    return months.map(function(month) {
      year = ee.Number(year);
      month = ee.Number(month);

      var start = ee.Date.fromYMD(year, month, 1);
      var end = start.advance(1, 'month');

      var monthly = chirps
        .filterDate(start, end)
        .sum()
        .rename('monthly_rainfall');

      return monthly
        .clip(roi)
        .set('system:time_start', start.millis())
        .set('year', year)
        .set('month', month);
    });
  }).flatten()
);



//  CALCULATE LONG-TERM CLIMATOLOGY BASELINE (1981-2025)
// ------------------------------------------------------------

// Hitung rerata dan standar deviasi historis bulanan (Jan-Des) selama 45 tahun
var monthlyClimatology = ee.ImageCollection.fromImages(
  months.map(function(m) {
    m = ee.Number(m);
    var filtered = monthlyRainfall.filter(ee.Filter.eq('month', m));

    var mean = filtered.select('monthly_rainfall').mean().rename('mean_rainfall');
    var stdDev = filtered.select('monthly_rainfall').reduce(ee.Reducer.stdDev()).rename('std_rainfall');

    return mean.addBands(stdDev).set('month', m);
  })
);



//  CALCULATE RAINFALL DEFICIT & SPI Z-SCORE (1981-2025)
// ------------------------------------------------------------

var monthlyDroughtCollection = monthlyRainfall.map(function(image) {
  var m = image.get('month');
  var clim = monthlyClimatology.filter(ee.Filter.eq('month', m)).first();
  
  var mean = ee.Image(clim).select('mean_rainfall');
  var stdDev = ee.Image(clim).select('std_rainfall');
  var actual = image.select('monthly_rainfall');

  // 1. Defisit Curah Hujan (mm) = Aktual - Rata-rata Normal
  var deficit = actual.subtract(mean).rename('rainfall_deficit');

  // 2. Persentase Normal (%) = (Aktual / Normal) * 100
  var percentNormal = actual.divide(mean).multiply(100).rename('percent_normal');

  // 3. Indeks Z-Score / Standar Deviasi (SPI-1 Proxy)
  var zScore = actual.subtract(mean)
    .divide(stdDev.where(stdDev.eq(0), 0.001))
    .rename('spi_zscore');

  return image
    .addBands(deficit)
    .addBands(percentNormal)
    .addBands(zScore);
});

print('========== MONTHLY DROUGHT COLLECTION ==========');
print('Expected monthly images:', 540);
print('Actual monthly images:', monthlyDroughtCollection.size());
print('Monthly band names:', ee.Image(monthlyDroughtCollection.first()).bandNames());



//  HISTORICAL TIME-SERIES CHARTS
// ------------------------------------------------------------

// 1. Curah Hujan Bulanan (1981-2025)
var monthlyChart = ui.Chart.image.series({
  imageCollection: monthlyDroughtCollection.select('monthly_rainfall'),
  region: roi,
  reducer: ee.Reducer.mean(),
  scale: 5566,
  xProperty: 'system:time_start'
}).setOptions({
  title: 'Historical Monthly Rainfall in Lampung (1981-2025)',
  hAxis: {title: 'Year', format: 'yyyy'},
  vAxis: {title: 'Rainfall (mm/month)'},
  lineWidth: 1,
  pointSize: 0
});

// 2. Defisit Curah Hujan (1981-2025)
var deficitChart = ui.Chart.image.series({
  imageCollection: monthlyDroughtCollection.select('rainfall_deficit'),
  region: roi,
  reducer: ee.Reducer.mean(),
  scale: 5566,
  xProperty: 'system:time_start'
}).setOptions({
  title: 'Historical Rainfall Deficit / Anomaly in Lampung (1981-2025)',
  hAxis: {title: 'Year', format: 'yyyy'},
  vAxis: {title: 'Rainfall Deficit (mm/month)'},
  colors: ['#d95f02'],
  lineWidth: 1,
  pointSize: 0
});

print('Historical monthly rainfall chart:', monthlyChart);
print('Historical rainfall deficit chart:', deficitChart);



//  SELECT MONTH FOR MAP VISUALIZATION
// ------------------------------------------------------------

var selectedYear = 2023; // Contoh: puncak kekeringan El Nino 2023
var selectedMonth = 9;   // September

var selectedDrought = ee.Image(
  monthlyDroughtCollection
    .filter(ee.Filter.eq('year', selectedYear))
    .filter(ee.Filter.eq('month', selectedMonth))
    .first()
);

print('========== SELECTED DROUGHT IMAGE ==========');
print('Selected year:', selectedYear);
print('Selected month:', selectedMonth);
print('Selected image:', selectedDrought);

var selectedStats = selectedDrought.reduceRegion({
  reducer: ee.Reducer.minMax().combine({
    reducer2: ee.Reducer.mean(),
    sharedInputs: true
  }),
  geometry: roi,
  scale: 5566,
  maxPixels: 1e9,
  tileScale: 4
});

print('Selected monthly drought statistics:', selectedStats);



//  MAP VISUALIZATION
// ------------------------------------------------------------

Map.centerObject(lampung, 8);

// 1. Curah Hujan Aktual
Map.addLayer(
  selectedDrought.select('monthly_rainfall'),
  {
    min: 0,
    max: 400,
    palette: ['fff7bc', 'fee391', 'fec44f', 'fe9929', 'ec7014', 'cc4c02', '993404']
  },
  'Monthly Rainfall ' + selectedYear + '-' + selectedMonth,
  false
);

// 2. Defisit Curah Hujan (Anomali mm)
Map.addLayer(
  selectedDrought.select('rainfall_deficit'),
  {
    min: -200,
    max: 200,
    palette: ['8c510a', 'd8b365', 'f6e8c3', 'f5f5f5', 'c7eae5', '5ab4ac', '01665e']
  },
  'Rainfall Deficit (mm) ' + selectedYear + '-' + selectedMonth,
  true
);

// 3. SPI Z-Score Kekeringan
Map.addLayer(
  selectedDrought.select('spi_zscore'),
  {
    min: -2.5,
    max: 2.5,
    palette: ['730000', 'e60000', 'ffaa00', 'ffff00', 'a800e6', '0070ff']
  },
  'SPI Z-Score ' + selectedYear + '-' + selectedMonth,
  true
);

// Layer Batas Administrasi
Map.addLayer(
  lampung.style({
    color: '000000',
    fillColor: '00000000',
    width: 2
  }),
  {},
  'Lampung Boundary',
  true
);



//  VALIDATION SUMMARY
// ------------------------------------------------------------

print('========== VALIDATION SUMMARY ==========');
print('Study area: Lampung Province');
print('Dataset: UCSB-CHG/CHIRPS/DAILY');
print('Historical period: 1981-2025 (45 Years Baseline)');
print('Actual daily images:', chirps.size());
print('Actual monthly images:', monthlyDroughtCollection.size());
print('Output bands: monthly_rainfall, rainfall_deficit, percent_normal, spi_zscore');
print('Nominal spatial resolution: approximately 5.5 km');
print('========================================');



//  EXPORT HISTORICAL DROUGHT STATISTICS TO CSV (1981-2025)
// ------------------------------------------------------------

var monthlyDroughtTable = ee.FeatureCollection(
  monthlyDroughtCollection.map(function(image) {
    var stats = image.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: roi,
      scale: 5566,
      maxPixels: 1e9,
      tileScale: 4
    });

    return ee.Feature(null, {
      'year': image.get('year'),
      'month': image.get('month'),
      'rainfall_actual_mm': stats.get('monthly_rainfall'),
      'rainfall_deficit_mm': stats.get('rainfall_deficit'),
      'percent_normal': stats.get('percent_normal'),
      'spi_zscore': stats.get('spi_zscore')
    });
  })
);

Export.table.toDrive({
  collection: monthlyDroughtTable,
  description: 'CHIRPS_Monthly_Drought_Indicators_1981_2025',
  folder: 'Agricultural_Drought_Lampung',
  fileNamePrefix: 'CHIRPS_Monthly_Drought_Indicators_1981_2025',
  fileFormat: 'CSV'
});



//  EXPORT SELECTED MONTH DROUGHT RASTER TO GEOTIFF
// ------------------------------------------------------------

Export.image.toDrive({
  image: selectedDrought.select(['monthly_rainfall', 'rainfall_deficit', 'spi_zscore']),
  description: 'CHIRPS_Drought_Raster_' + selectedYear + '_' + selectedMonth,
  folder: 'Agricultural_Drought_Lampung',
  fileNamePrefix: 'CHIRPS_Drought_Raster_' + selectedYear + '_' + selectedMonth,
  region: roi,
  scale: 5566,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});