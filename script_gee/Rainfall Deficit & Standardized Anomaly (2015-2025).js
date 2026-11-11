// PROJECT: AGRICULTURAL DROUGHT MONITORING - LAMPUNG
// SCRIPT: RAINFALL DEFICIT & STANDARDIZED ANOMALY (2015-2025)
// BASELINE CLIMATOLOGY: 1981-2025 (45 YEARS)
// ------------------------------------------------------------


//  STUDY AREA SETUP
// ------------------------------------------------------------
var admin = ee.FeatureCollection('FAO/GAUL/2015/level1');
var lampung = admin
  .filter(ee.Filter.eq('ADM0_NAME', 'Indonesia'))
  .filter(ee.Filter.eq('ADM1_NAME', 'Lampung'));
var roi = lampung.geometry();


//  LOAD & ACCUMULATE CHIRPS DAILY TO MONTHLY (1981-2025)
// ------------------------------------------------------------
var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .filterDate('1981-01-01', '2026-01-01')
  .filterBounds(roi)
  .select('precipitation');

var histYears = ee.List.sequence(1981, 2025);
var months = ee.List.sequence(1, 12);


var monthlyRainfallAll = ee.ImageCollection.fromImages(
  histYears.map(function(year) {
    return months.map(function(month) {
      year = ee.Number(year);
      month = ee.Number(month);
      var start = ee.Date.fromYMD(year, month, 1);
      var end = start.advance(1, 'month');

      var monthlySum = chirps
        .filterDate(start, end)
        .sum()
        .rename('monthly_rainfall');

      return monthlySum
        .clip(roi)
        .set('system:time_start', start.millis())
        .set('year', year)
        .set('month', month);
    });
  }).flatten()
);


//  CALCULATE LONG-TERM CLIMATOLOGY BASELINE (1981-2025)
// ------------------------------------------------------------

var monthlyClimatology = ee.ImageCollection.fromImages(
  months.map(function(m) {
    m = ee.Number(m);
    var filtered = monthlyRainfallAll.filter(ee.Filter.eq('month', m));
    
    var mean = filtered.select('monthly_rainfall').mean().rename('mean_rainfall');
    var stdDev = filtered.select('monthly_rainfall').reduce(ee.Reducer.stdDev()).rename('std_rainfall');
    
    return mean.addBands(stdDev).set('month', m);
  })
);


//  CALCULATE RAINFALL DEFICIT & SPI Z-SCORE (2015-2025)
// ------------------------------------------------------------
var analysisYears = ee.List.sequence(2015, 2025);

var droughtCollection = ee.ImageCollection.fromImages(
  analysisYears.map(function(year) {
    return months.map(function(month) {
      year = ee.Number(year);
      month = ee.Number(month);
      
      var start = ee.Date.fromYMD(year, month, 1);
      
      // curah hujan aktual bulan tersebut
      var actual = monthlyRainfallAll
        .filter(ee.Filter.eq('year', year))
        .filter(ee.Filter.eq('month', month))
        .first();
        
      // baseline klimatologi bulan tersebut
      var clim = monthlyClimatology
        .filter(ee.Filter.eq('month', month))
        .first();
        
      var mean = clim.select('mean_rainfall');
      var stdDev = clim.select('std_rainfall');
      var rainVal = ee.Image(actual).select('monthly_rainfall');
      
      // 1. Defisit Absolut (Anomali mm) = Aktual - Normal
      var deficit = rainVal.subtract(mean).rename('rainfall_deficit');
      
      // 2. Persentase Normal (%) = (Aktual / Normal) * 100
      var percentNormal = rainVal.divide(mean).multiply(100).rename('percent_normal');
      
      // 3. Z-Score (Standar Deviasi / Proxy SPI-1) = (Aktual - Mean) / StdDev
      // Menghindari pembagian dengan 0 menggunakan unmask
      var zScore = rainVal.subtract(mean)
        .divide(stdDev.where(stdDev.eq(0), 0.001))
        .rename('spi_zscore');
        
      return rainVal
        .addBands(deficit)
        .addBands(percentNormal)
        .addBands(zScore)
        .set('system:time_start', start.millis())
        .set('year', year)
        .set('month', month);
    });
  }).flatten()
);

print('========== DROUGHT COLLECTION (2015-2025) ==========');
print('Total monthly drought images (11 years x 12):', droughtCollection.size());
print('Bands available:', ee.Image(droughtCollection.first()).bandNames());


//  TIME-SERIES CHARTS (2015-2025)
// ------------------------------------------------------------

// 1. Defisit Curah Hujan Bulanan (mm)
var deficitChart = ui.Chart.image.series({
  imageCollection: droughtCollection.select('rainfall_deficit'),
  region: roi,
  reducer: ee.Reducer.mean(),
  scale: 5566,
  xProperty: 'system:time_start'
}).setOptions({
  title: 'Rata-rata Defisit Curah Hujan Bulanan Lampung (2015–2025)',
  hAxis: {title: 'Tahun', format: 'yyyy'},
  vAxis: {title: 'Defisit Curah Hujan (mm/bulan)'},
  colors: ['#d95f02'],
  lineWidth: 1.5,
  pointSize: 2
});

// 2. Indeks Z-Score / Standardized Anomaly
var zScoreChart = ui.Chart.image.series({
  imageCollection: droughtCollection.select('spi_zscore'),
  region: roi,
  reducer: ee.Reducer.mean(),
  scale: 5566,
  xProperty: 'system:time_start'
}).setOptions({
  title: 'Indeks Kekeringan Z-Score (CHIRPS) Lampung (2015–2025)',
  hAxis: {title: 'Tahun', format: 'yyyy'},
  vAxis: {title: 'Z-Score ( < -1 = Kering, > 1 = Basah)'},
  colors: ['#7570b3'],
  lineWidth: 1.5,
  pointSize: 2
});

print(deficitChart);
print(zScoreChart);


//  MAP VISUALIZATION (SAMPLE EXTREME DROUGHT: SEPT 2023 EL NINO)
// ------------------------------------------------------------
var selYear = 2023;
var selMonth = 9;

var selectedImg = ee.Image(
  droughtCollection
    .filter(ee.Filter.eq('year', selYear))
    .filter(ee.Filter.eq('month', selMonth))
    .first()
);

Map.centerObject(lampung, 8);

// 1. Curah Hujan Aktual
Map.addLayer(
  selectedImg.select('monthly_rainfall'),
  {min: 0, max: 300, palette: ['#fff7bc', '#fec44f', '#d95f0e', '#045a8d']},
  '1. Curah Hujan Aktual (mm) - ' + selYear + '/' + selMonth,
  false
);

// 2. Defisit Curah Hujan (Anomali mm)
Map.addLayer(
  selectedImg.select('rainfall_deficit'),
  {min: -200, max: 200, palette: ['#8c510a', '#d8b365', '#f6e8c3', '#f5f5f5', '#c7eae5', '#5ab4ac', '#01665e']},
  '2. Defisit Curah Hujan (mm) - ' + selYear + '/' + selMonth,
  true
);

// 3. Kelas Kekeringan Z-Score
Map.addLayer(
  selectedImg.select('spi_zscore'),
  {min: -2.5, max: 2.5, palette: ['#730000', '#e60000', '#ffaa00', '#ffff00', '#a800e6', '#0070ff']},
  '3. Z-Score Kekeringan - ' + selYear + '/' + selMonth,
  true
);

// Layer Boundary
Map.addLayer(
  lampung.style({color: '000000', fillColor: '00000000', width: 2}),
  {},
  'Batas Admin Lampung'
);


//  EXPORT DATA TO DRIVE
// ------------------------------------------------------------

// 1. Ekspor Ringkasan Statistik Bulanan (CSV)
var droughtTable = ee.FeatureCollection(
  droughtCollection.map(function(img) {
    var stats = img.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: roi,
      scale: 5566,
      maxPixels: 1e9,
      tileScale: 4
    });
    
    return ee.Feature(null, {
      'year': img.get('year'),
      'month': img.get('month'),
      'rainfall_actual_mm': stats.get('monthly_rainfall'),
      'rainfall_deficit_mm': stats.get('rainfall_deficit'),
      'percent_normal': stats.get('percent_normal'),
      'spi_zscore': stats.get('spi_zscore')
    });
  })
);

Export.table.toDrive({
  collection: droughtTable,
  description: 'CHIRPS_Drought_Indicators_Lampung_2015_2025',
  folder: 'Agricultural_Drought_Lampung',
  fileNamePrefix: 'CHIRPS_Drought_Indicators_Lampung_2015_2025',
  fileFormat: 'CSV'
});

// 2. Ekspor Raster Defisit September 2023 (GeoTIFF)
Export.image.toDrive({
  image: selectedImg.select(['rainfall_deficit', 'spi_zscore']),
  description: 'Raster_Deficit_Lampung_2023_09',
  folder: 'Agricultural_Drought_Lampung',
  fileNamePrefix: 'Raster_Deficit_Lampung_2023_09',
  region: roi,
  scale: 5566,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});