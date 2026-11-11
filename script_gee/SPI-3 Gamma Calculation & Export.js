// PHASE 4: SCIENTIFIC SPI-3 GAMMA DISTRIBUTION (CHIRPS)
// SCRIPT 04: SPI-3 GAMMA CALCULATION & EXPORT
// ============================================================



//  STUDY AREA SETUP
// ------------------------------------------------------------

var admin = ee.FeatureCollection('FAO/GAUL/2015/level1');

var lampung = admin
  .filter(ee.Filter.eq('ADM0_NAME', 'Indonesia'))
  .filter(ee.Filter.eq('ADM1_NAME', 'Lampung'));

var roi = lampung.geometry();



//  LOAD CHIRPS DAILY & ACCUMULATE TO MONTHLY (1981-2025)
// ------------------------------------------------------------

var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .filterDate('1981-01-01', '2026-01-01')
  .filterBounds(roi)
  .select('precipitation');

var years = ee.List.sequence(1981, 2025);
var months = ee.List.sequence(1, 12);

var monthlyRainfall = ee.ImageCollection.fromImages(
  years.map(function(year) {
    return months.map(function(month) {
      year = ee.Number(year);
      month = ee.Number(month);

      var start = ee.Date.fromYMD(year, month, 1);
      var end = start.advance(1, 'month');

      var sum = chirps
        .filterDate(start, end)
        .sum()
        .rename('monthly_rainfall');

      return sum
        .clip(roi)
        .set('year', year)
        .set('month', month)
        .set('system:time_start', start.millis());
    });
  }).flatten()
);



//  ROLLING 3-MONTH RAINFALL ACCUMULATION FUNCTION
// ------------------------------------------------------------

// Fungsi untuk membuat akumulasi 3 bulanan (Bulan T-2, T-1, T)
var get3MonthRainfall = function(targetYear, targetMonth) {
  var endDate = ee.Date.fromYMD(targetYear, targetMonth, 1).advance(1, 'month');
  var startDate = endDate.advance(-3, 'month');
  
  var sum3m = chirps
    .filterDate(startDate, endDate)
    .sum()
    .rename('rainfall_3month');
    
  return sum3m
    .clip(roi)
    .set('year', targetYear)
    .set('month', targetMonth);
};



//  CONFIGURATION: TARGET PERIOD & BASELINE
// ------------------------------------------------------------

var targetYear = 2025;
var targetMonth = 9; // September 2025 (Akumulasi Juli - Agustus - September)

var baselineStartYear = 1981;
var baselineEndYear = 2020;
var baselineYears = ee.List.sequence(baselineStartYear, baselineEndYear);

// Koleksi data historis 3-bulanan khusus bulan target (misal: semua September 1981-2020)
var baselineSep = ee.ImageCollection.fromImages(
  baselineYears.map(function(y) {
    return get3MonthRainfall(ee.Number(y), targetMonth);
  })
);

// Data 3-bulanan aktual yang akan dihitung SPI-3 nya
var rainfallTarget = get3MonthRainfall(targetYear, targetMonth);



//  GAMMA DISTRIBUTION PARAMETERS ESTIMATION (THOM'S MLE)
// ------------------------------------------------------------

// 1. Filter sampel positif (> 0 mm)
var positiveRainfall = baselineSep.map(function(img) {
  return img.updateMask(img.gt(0));
});

var positiveCount = positiveRainfall.count();
var totalCount = baselineSep.count();

// 2. Rata-rata aritmatika dan logaritma dari observasi positif
var meanRain = positiveRainfall.mean();
var meanLogRain = positiveRainfall.map(function(img) {
  return img.log().rename('log_rain');
}).mean();

// 3. Thom's Approximation untuk Parameter Alpha (Shape) dan Beta (Scale)
var A = meanRain.log().subtract(meanLogRain);

var alpha = ee.Image.constant(1)
  .add(
    ee.Image.constant(1)
      .add(A.multiply(4 / 3))
      .sqrt()
  )
  .divide(A.multiply(4))
  .rename('alpha');

var beta = meanRain.divide(alpha).rename('beta');



// PROBABILITY OF ZERO RAINFALL (q) & CUMULATIVE DENSITY (CDF)
// ------------------------------------------------------------

var zeroCount = baselineSep.map(function(img) {
  return img.eq(0);
}).sum();

var q = zeroCount.divide(totalCount).rename('q');
var standardizedRain = rainfallTarget.divide(beta);
var gammaCDF = alpha.gammainc(standardizedRain);
var probability = q.add(
  ee.Image.constant(1).subtract(q).multiply(gammaCDF)
).clamp(0.000001, 0.999999);



//  CONVERT TO STANDARDIZED PRECIPITATION INDEX (SPI-3)
// ------------------------------------------------------------

var validMask = A.gt(0)
  .and(alpha.gt(0))
  .and(beta.gt(0))
  .and(positiveCount.gte(30))
  .and(totalCount.gte(30));

var spi3 = probability
  .multiply(2)
  .subtract(1)
  .erfInv()
  .multiply(Math.sqrt(2))
  .rename('SPI3_Gamma')
  .updateMask(validMask)
  .clip(roi)
  .toFloat();


//  MAP VISUALIZATION
// ------------------------------------------------------------

Map.centerObject(lampung, 8);

var spiVis = {
  min: -2.5,
  max: 2.5,
  palette: [
    '#67001f', // Sangat Amat Kering (Extreme Drought)
    '#b2182b', // Kering Berat (Severe Drought)
    '#ef8a62', // Kering Sedang (Moderate Drought)
    '#fddbc7', // Agak Kering
    '#f7f7f7', // Normal
    '#d1e5f0', // Agak Basah
    '#67a9cf', // Basah Sedang
    '#2166ac', // Basah Berat
    '#053061'  // Sangat Basah
  ]
};

Map.addLayer(
  spi3,
  spiVis,
  'SPI-3 Gamma ' + targetYear + '-' + targetMonth
);

Map.addLayer(
  lampung.style({color: 'black', fillColor: '00000000', width: 2}),
  {},
  'Lampung Boundary'
);


// STATISTICAL REDUCTION & EXPORT
// ------------------------------------------------------------

var spiStats = spi3.reduceRegion({
  reducer: ee.Reducer.minMax().combine({
    reducer2: ee.Reducer.mean(),
    sharedInputs: true
  }),
  geometry: roi,
  scale: 5566,
  maxPixels: 1e9,
  tileScale: 4
});

print('========== SPI-3 GAMMA VALIDATION ==========');
print('Target Year:', targetYear);
print('Target Month:', targetMonth);
print('Baseline Period:', baselineStartYear + '-' + baselineEndYear);
print('SPI-3 Gamma Statistics:', spiStats);

// 1. Export GeoTIFF Raster
Export.image.toDrive({
  image: spi3,
  description: 'SPI3_Gamma_' + targetYear + '_' + (targetMonth < 10 ? '0' + targetMonth : targetMonth),
  folder: 'Agricultural_Drought_Lampung',
  fileNamePrefix: 'SPI3_Gamma_' + targetYear + '_' + (targetMonth < 10 ? '0' + targetMonth : targetMonth),
  region: roi,
  scale: 5566,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

// 2. Export CSV Statistics
var spiFeature = ee.Feature(null, {
  year: targetYear,
  month: targetMonth,
  baseline: baselineStartYear + '-' + baselineEndYear,
  SPI3_min: spiStats.get('SPI3_Gamma_min'),
  SPI3_max: spiStats.get('SPI3_Gamma_max'),
  SPI3_mean: spiStats.get('SPI3_Gamma_mean')
});

Export.table.toDrive({
  collection: ee.FeatureCollection([spiFeature]),
  description: 'SPI3_Gamma_Statistics_' + targetYear + '_' + (targetMonth < 10 ? '0' + targetMonth : targetMonth),
  folder: 'Agricultural_Drought_Lampung',
  fileNamePrefix: 'SPI3_Gamma_Statistics_' + targetYear + '_' + (targetMonth < 10 ? '0' + targetMonth : targetMonth),
  fileFormat: 'CSV'
});