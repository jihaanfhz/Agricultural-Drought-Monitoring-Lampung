// PHASE 3: AGRICULTURAL DROUGHT ASSESSMENT (CHIRPS)
// SCRIPT 03: ROLLING 3-MONTH RAINFALL DEFICIT & SPI-3 (1981-2025)
// ============================================================



//  STUDY AREA SETUP
// ------------------------------------------------------------

var admin = ee.FeatureCollection('FAO/GAUL/2015/level1');

var lampung = admin
  .filter(ee.Filter.eq('ADM0_NAME', 'Indonesia'))
  .filter(ee.Filter.eq('ADM1_NAME', 'Lampung'));

var roi = lampung.geometry();



//  LOAD CHIRPS DAILY (1981-2025)
// ------------------------------------------------------------

var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .filterDate('1981-01-01', '2026-01-01')
  .filterBounds(roi)
  .select('precipitation');


// 3. CREATE ROLLING 3-MONTH RAINFALL ACCUMULATION
// ------------------------------------------------------------

var years = ee.List.sequence(1981, 2025);
var months = ee.List.sequence(1, 12);

var rainfall3MonthAll = ee.ImageCollection.fromImages(
  years.map(function(year) {
    return months.map(function(month) {
      year = ee.Number(year);
      month = ee.Number(month);

      var currentDate = ee.Date.fromYMD(year, month, 1);
      var windowStart = currentDate.advance(-2, 'month');
      var windowEnd = currentDate.advance(1, 'month');

      var rainfall3 = chirps
        .filterDate(windowStart, windowEnd)
        .sum()
        .rename('rainfall_3month');

      return rainfall3
        .clip(roi)
        .set('system:time_start', currentDate.millis())
        .set('year', year)
        .set('month', month)
        .set('window_start', windowStart.format('YYYY-MM-dd'))
        .set('window_end', windowEnd.advance(-1, 'day').format('YYYY-MM-dd'));
    });
  }).flatten()
).filter(ee.Filter.gte('system:time_start', ee.Date('1981-03-01').millis()));


//  CALCULATE LONG-TERM CLIMATOLOGY BASELINE (3-MONTH 1981-2025)
// ------------------------------------------------------------

// Hitung Rerata dan Standar Deviasi 3-bulanan historis (1981-2025)
var monthlyClimatology3M = ee.ImageCollection.fromImages(
  months.map(function(m) {
    m = ee.Number(m);
    var filtered = rainfall3MonthAll.filter(ee.Filter.eq('month', m));

    var mean = filtered.select('rainfall_3month').mean().rename('mean_rainfall_3m');
    var stdDev = filtered.select('rainfall_3month').reduce(ee.Reducer.stdDev()).rename('std_rainfall_3m');

    return mean.addBands(stdDev).set('month', m);
  })
);



//  CALCULATE 3-MONTH DEFICIT & SPI-3 Z-SCORE
// ------------------------------------------------------------

var drought3MonthCollection = rainfall3MonthAll.map(function(image) {
  var m = image.get('month');
  var clim = monthlyClimatology3M.filter(ee.Filter.eq('month', m)).first();

  var mean = ee.Image(clim).select('mean_rainfall_3m');
  var stdDev = ee.Image(clim).select('std_rainfall_3m');
  var actual = image.select('rainfall_3month');

  // 1. Defisit Curah Hujan 3 Bulanan (mm)
  var deficit = actual.subtract(mean).rename('rainfall_3month_deficit');

  // 2. Persentase Normal 3 Bulanan (%)
  var percentNormal = actual.divide(mean).multiply(100).rename('percent_normal_3m');

  // 3. Standar Deviasi / SPI-3 Z-Score (Indeks Kekeringan Pertanian)
  var spi3Zscore = actual.subtract(mean)
    .divide(stdDev.where(stdDev.eq(0), 0.001))
    .rename('spi3_zscore');

  return image
    .addBands(deficit)
    .addBands(percentNormal)
    .addBands(spi3Zscore);
});

print('========== ROLLING 3-MONTH DROUGHT COLLECTION ==========');
print('Expected images:', 538);
print('Actual images:', drought3MonthCollection.size());
print('Bands available:', ee.Image(drought3MonthCollection.first()).bandNames());



//  HISTORICAL TIME-SERIES CHARTS
// ------------------------------------------------------------

// 1. SPI-3 Z-Score (Indeks Utama Kekeringan Pertanian 1981-2025)
var spi3Chart = ui.Chart.image.series({
  imageCollection: drought3MonthCollection.select('spi3_zscore'),
  region: roi,
  reducer: ee.Reducer.mean(),
  scale: 5566,
  xProperty: 'system:time_start'
}).setOptions({
  title: 'Indeks SPI-3 (Agricultural Drought) Lampung (1981-2025)',
  hAxis: {title: 'Year', format: 'yyyy'},
  vAxis: {title: 'SPI-3 Z-Score ( < -1 = Kering Pertanian)'},
  colors: ['#e66101'],
  lineWidth: 1,
  pointSize: 0
});

print('Historical SPI-3 chart:', spi3Chart);



//  SELECT MONTH FOR MAP VISUALIZATION
// ------------------------------------------------------------

var selectedYear = 2023; // Contoh: Puncak El Nino 2023
var selectedMonth = 9;   // September (Akumulasi Juli-Agustus-September)

var selectedDrought3 = ee.Image(
  drought3MonthCollection
    .filter(ee.Filter.eq('year', selectedYear))
    .filter(ee.Filter.eq('month', selectedMonth))
    .first()
);

print('========== SELECTED 3-MONTH DROUGHT ==========');
print('Year:', selectedYear);
print('Ending Month:', selectedMonth);
print('Selected Image:', selectedDrought3);

var selectedStats = selectedDrought3.reduceRegion({
  reducer: ee.Reducer.minMax().combine({
    reducer2: ee.Reducer.mean(),
    sharedInputs: true
  }),
  geometry: roi,
  scale: 5566,
  maxPixels: 1e9,
  tileScale: 4
});

print('Selected 3-month drought statistics:', selectedStats);



// 8. MAP VISUALIZATION
// ------------------------------------------------------------

Map.centerObject(lampung, 8);

// 1. Hujan 3-Bulan Aktual (mm)
Map.addLayer(
  selectedDrought3.select('rainfall_3month'),
  {
    min: 0,
    max: 1200,
    palette: ['fff7bc', 'fee391', 'fec44f', 'fe9929', 'ec7014', 'cc4c02', '993404']
  },
  '1. Hujan 3-Bulan Aktual (mm) - ' + selectedYear + '/' + selectedMonth,
  false
);

// 2. Defisit Hujan 3-Bulan (mm)
Map.addLayer(
  selectedDrought3.select('rainfall_3month_deficit'),
  {
    min: -400,
    max: 400,
    palette: ['8c510a', 'd8b365', 'f6e8c3', 'f5f5f5', 'c7eae5', '5ab4ac', '01665e']
  },
  '2. Defisit Hujan 3-Bulan (mm) - ' + selectedYear + '/' + selectedMonth,
  true
);

// 3. SPI-3 Z-Score (Kekeringan Pertanian)
Map.addLayer(
  selectedDrought3.select('spi3_zscore'),
  {
    min: -2.5,
    max: 2.5,
    palette: ['730000', 'e60000', 'ffaa00', 'ffff00', 'a800e6', '0070ff']
  },
  '3. SPI-3 Z-Score (Kekeringan Pertanian) - ' + selectedYear + '/' + selectedMonth,
  true
);

// Layer Batas Administrasi
Map.addLayer(
  lampung.style({
    color: 'black',
    fillColor: '00000000',
    width: 2
  }),
  {},
  'Batas Admin Lampung',
  true
);


//  VALIDATION SUMMARY
// ------------------------------------------------------------

print('========== VALIDATION SUMMARY ==========');
print('Dataset: CHIRPS Daily');
print('Period: March 1981 - December 2025');
print('Expected images: 538');
print('Actual images:', drought3MonthCollection.size());
print('Output bands: rainfall_3month, rainfall_3month_deficit, percent_normal_3m, spi3_zscore');
print('========================================');


//  EXPORT 3-MONTH DROUGHT INDICATORS TO CSV (1981-2025)
// ------------------------------------------------------------

var drought3Table = ee.FeatureCollection(
  drought3MonthCollection.map(function(image) {
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
      'rainfall_3month_actual_mm': stats.get('rainfall_3month'),
      'rainfall_3month_deficit_mm': stats.get('rainfall_3month_deficit'),
      'percent_normal_3m': stats.get('percent_normal_3m'),
      'spi3_zscore': stats.get('spi3_zscore')
    });
  })
);

Export.table.toDrive({
  collection: drought3Table,
  description: 'CHIRPS_SPI3_Drought_Indicators_1981_2025',
  folder: 'Agricultural_Drought_Lampung',
  fileNamePrefix: 'CHIRPS_SPI3_Drought_Indicators_1981_2025',
  fileFormat: 'CSV'
});


//  EXPORT SELECTED MONTH SPI-3 RASTER TO GEOTIFF
// ------------------------------------------------------------

Export.image.toDrive({
  image: selectedDrought3.select(['rainfall_3month', 'rainfall_3month_deficit', 'spi3_zscore']),
  description: 'CHIRPS_SPI3_Raster_' + selectedYear + '_' + selectedMonth,
  folder: 'Agricultural_Drought_Lampung',
  fileNamePrefix: 'CHIRPS_SPI3_Raster_' + selectedYear + '_' + selectedMonth,
  region: roi,
  scale: 5566,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});