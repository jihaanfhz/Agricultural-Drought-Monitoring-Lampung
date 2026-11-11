// PHASE 6: MASTER HISTORICAL SPI-3 PROCESSING & STACKING
// SCRIPT 06: AUTOMATED HISTORICAL SPI-3 GAMMA (1981-2025) & STACK (2015-2025)
// ============================================================


// STUDY AREA SETUP
// ------------------------------------------------------------
var lampung = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filter(ee.Filter.eq('ADM0_NAME', 'Indonesia'))
  .filter(ee.Filter.eq('ADM1_NAME', 'Lampung'));

var region = lampung.geometry();

Map.centerObject(lampung, 8);
Map.addLayer(lampung, {color: 'black'}, 'Lampung Boundary');

print('========== STUDY AREA ==========');
print('Boundary count:', lampung.size());


// CHIRPS DAILY & MONTHLY RAINFALL (1981-2025)
// ------------------------------------------------------------
var analysisStartYear = 1981;
var analysisEndYear = 2025;
var baselineStartYear = 1981;
var baselineEndYear = 2020;

var daily = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .filterDate(analysisStartYear + '-01-01', (analysisEndYear + 1) + '-01-01')
  .select('precipitation');

print('========== CONFIGURATION ==========');
print('Analysis Period:', analysisStartYear + '-' + analysisEndYear);
print('Baseline Period:', baselineStartYear + '-' + baselineEndYear);
print('CHIRPS daily images:', daily.size());

var years = ee.List.sequence(analysisStartYear, analysisEndYear);
var months = ee.List.sequence(1, 12);

var monthlyList = years.map(function(year) {
  year = ee.Number(year);
  return months.map(function(month) {
    month = ee.Number(month);
    var startDate = ee.Date.fromYMD(year, month, 1);
    var endDate = startDate.advance(1, 'month');

    var rainfall = daily
      .filterDate(startDate, endDate)
      .sum()
      .rename('monthly_rainfall');

    return rainfall.set({
      'year': year,
      'month': month,
      'system:time_start': startDate.millis()
    });
  });
}).flatten();

var monthly = ee.ImageCollection.fromImages(monthlyList);
print('Monthly rainfall count (540 images):', monthly.size());


// THREE-MONTH RAINFALL ACCUMULATION FUNCTION
// ------------------------------------------------------------
function rainfall3Month(year, month) {
  year = ee.Number(year);
  month = ee.Number(month);

  var endDate = ee.Date.fromYMD(year, month, 1);
  var startDate = endDate.advance(-2, 'month');

  var threeMonths = monthly.filterDate(startDate, endDate.advance(1, 'month'));
  var rainfall = threeMonths.sum().rename('rainfall_3month');
  var validThreeMonths = threeMonths.count().eq(3);

  rainfall = rainfall.updateMask(validThreeMonths);

  return rainfall.set({
    'year': year,
    'month': month,
    'system:time_start': endDate.millis()
  });
}

// SPI-3 GAMMA CALCULATION FUNCTION (FIXED PROPERTY & SYNTAX)
// ------------------------------------------------------------
function calculateSPI(year, month) {
  year = ee.Number(year);
  month = ee.Number(month);

  var firstBaselineYear = ee.Number(
    ee.Algorithms.If(month.lte(2), baselineStartYear + 1, baselineStartYear)
  );

  var baselineYears = ee.List.sequence(firstBaselineYear, baselineEndYear);

  var baseline = ee.ImageCollection.fromImages(
    baselineYears.map(function(y) {
      return rainfall3Month(y, month);
    })
  );

  var target = rainfall3Month(year, month);
  var totalCount = baseline.count();

  var positive = baseline.map(function(image) {
    return image.updateMask(image.gt(0));
  });

  var positiveCount = positive.count();
  var meanRainfall = positive.mean();
  var meanLog = positive.map(function(image) { return image.log(); }).mean();

  // Estimasi parameter Gamma (Thom approximation)
  var A = meanRainfall.log().subtract(meanLog);
  var alpha = ee.Image(1).add(ee.Image(1).add(A.multiply(4 / 3)).sqrt()).divide(A.multiply(4));
  var beta = meanRainfall.divide(alpha);

  // Probabilitas hujan nol
  var zeroCount = baseline.map(function(image) { return image.eq(0); }).sum();
  var q = zeroCount.divide(totalCount);

  // Sintaks GEE: alpha.gammainc(standardizedTarget)
  var standardizedTarget = target.divide(beta);
  var gammaCDF = alpha.gammainc(standardizedTarget);

  // Probabilitas Kumulatif Campuran H(x)
  var probability = q.add(ee.Image(1).subtract(q).multiply(gammaCDF)).clamp(0.000001, 0.999999);

  // Inversi Transformasi Normal Standar -> SPI3 Z-Score
  var spi = probability.multiply(2).subtract(1).erfInv().multiply(Math.sqrt(2));

  // Masking Validasi
  var validMask = A.gt(0)
    .and(alpha.gt(0))
    .and(beta.gt(0))
    .and(positiveCount.gte(30))
    .and(totalCount.gte(30));

  spi = spi.rename('SPI3')
    .updateMask(validMask)
    .clip(region);

  // Membuat properti nama band secara eksplisit
  var yStr = year.format('%04d');
  var mStr = month.format('%02d');
  var bandName = ee.String('SPI3_').cat(yStr).cat('_').cat(mStr);

  return spi.set({
    'year': year,
    'month': month,
    'band_name': bandName,
    'baseline_year_count': baselineYears.size(),
    'system:time_start': ee.Date.fromYMD(year, month, 1).millis()
  });
}

// AUTOMATED SPI-3 PROCESSING FOR ALL YEARS (1981-2025)
// ------------------------------------------------------------
var allHistoricalSPIList = years.map(function(year) {
  year = ee.Number(year);
  return months.map(function(month) {
    month = ee.Number(month);
    return calculateSPI(year, month);
  });
}).flatten();

var allHistoricalSPI = ee.ImageCollection.fromImages(allHistoricalSPIList);

print('========== HISTORICAL SPI-3 COLLECTION ==========');
print('Total monthly SPI-3 images (1981-2025):', allHistoricalSPI.size());


// MAP VISUALIZATION (SAMPLE: SEPT 2023 EL NINO EXTREME)
// ------------------------------------------------------------
var sampleSPI = ee.Image(
  allHistoricalSPI
    .filter(ee.Filter.eq('year', 2023))
    .filter(ee.Filter.eq('month', 9))
    .first()
);

var spiPalette = [
  '8c2d04', 'd7301f', 'fc8d59', 'fee08b',
  'ffffbf', 'd9ef8b', '91cf60', '1a9850', '006837'
];

Map.addLayer(
  sampleSPI,
  {min: -3, max: 3, palette: spiPalette},
  'SPI-3 Sample (September 2023)'
);


// TIME-SERIES CHART (HISTORICAL TREND 1981-2025)
// ------------------------------------------------------------
var historicalChart = ui.Chart.image.series({
  imageCollection: allHistoricalSPI,
  region: region,
  reducer: ee.Reducer.mean(),
  scale: 5566,
  xProperty: 'system:time_start'
}).setOptions({
  title: 'Tren SPI-3 Historis Provinsi Lampung (1981–2025)',
  hAxis: {title: 'Tahun', format: 'yyyy'},
  vAxis: {title: 'Nilai SPI-3 ( < -1 = Kering Pertanian)'},
  colors: ['#d7301f'],
  lineWidth: 1.2,
  pointSize: 0
});

print(historicalChart);


// AUTOMATED HISTORICAL STATISTICAL TABLE (1981-2025)
// ------------------------------------------------------------
var historicalStatisticsTable = ee.FeatureCollection(
  allHistoricalSPI.map(function(image) {
    var y = image.get('year');
    var m = image.get('month');

    var stats = image.reduceRegion({
      reducer: ee.Reducer.minMax()
        .combine({reducer2: ee.Reducer.mean(), sharedInputs: true})
        .combine({reducer2: ee.Reducer.count(), sharedInputs: true}),
      geometry: region,
      scale: 5566,
      maxPixels: 1e13,
      tileScale: 4
    });

    var validPixels = ee.Number(
      ee.Algorithms.If(stats.contains('SPI3_count'), stats.get('SPI3_count'), 0)
    );
    var hasValidData = validPixels.gt(0);

    return ee.Feature(null, {
      'year': y,
      'month': m,
      'baseline_year_count': image.get('baseline_year_count'),
      'valid_pixel_count': validPixels,
      'SPI3_min': ee.Algorithms.If(hasValidData, stats.get('SPI3_min'), -9999),
      'SPI3_max': ee.Algorithms.If(hasValidData, stats.get('SPI3_max'), -9999),
      'SPI3_mean': ee.Algorithms.If(hasValidData, stats.get('SPI3_mean'), -9999),
      'status': ee.Algorithms.If(hasValidData, 'VALID', 'NO_VALID_DATA')
    });
  })
);

print('========== STATISTICS TABLE SUMMARY ==========');
print('Total monthly rows for CSV (1981-2025):', historicalStatisticsTable.size());

// EXPORT 1: FULL HISTORICAL STATISTICS CSV (1981-2025)
// ------------------------------------------------------------
Export.table.toDrive({
  collection: historicalStatisticsTable,
  description: 'SPI3_Monthly_Statistics_Historical_1981_2025',
  folder: 'Agricultural_Drought_Lampung',
  fileNamePrefix: 'SPI3_Monthly_Statistics_Historical_1981_2025',
  fileFormat: 'CSV'
});


// EXPORT 2: TARGET PERIOD GEOTIFF RASTER STACK (2015-2025)
// ------------------------------------------------------------
var targetPeriodSPI = allHistoricalSPI
  .filter(ee.Filter.gte('year', 2015))
  .sort('system:time_start');

// Mengambil list nama band langsung dari properti metadata 'band_name'
var targetBandNames = targetPeriodSPI.aggregate_array('band_name');

// Mengubah ImageCollection ke multi-band Image lalu merename band secara presisi
var targetStackRenamed = targetPeriodSPI.toBands().rename(targetBandNames);

Export.image.toDrive({
  image: targetStackRenamed,
  description: 'SPI3_Monthly_Stack_2015_2025',
  folder: 'Agricultural_Drought_Lampung',
  fileNamePrefix: 'SPI3_Monthly_Stack_2015_2025',
  region: region,
  scale: 5566,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

print('========== EXPORT TASKS READY ==========');
print('1. File CSV: 540 baris data statistik SPI-3 (1981-2025)');
print('2. File GeoTIFF Stack: 132 band raster SPI-3 (2015-2025)');