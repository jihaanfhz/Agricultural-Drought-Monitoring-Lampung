// PHASE 5: FULL-YEAR SPI-3 PROCESSING & STACKING (2025)
// SCRIPT 05: MONTHLY SPI-3 GAMMA (JAN-DEC 2025) & 12-BAND GEOTIFF STACK
// ============================================================


// STUDY AREA
// ------------------------------------------------------------

var lampung = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filter(ee.Filter.eq('ADM0_NAME', 'Indonesia'))
  .filter(ee.Filter.eq('ADM1_NAME', 'Lampung'));

var region = lampung.geometry();

Map.centerObject(lampung, 8);
Map.addLayer(lampung, {color: 'black'}, 'Lampung Boundary');

print('========== STUDY AREA ==========');
print('Boundary count:', lampung.size());



// CHIRPS DAILY AND MONTHLY RAINFALL (1981-2025)
// ------------------------------------------------------------

var daily = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .filterDate('1981-01-01', '2026-01-01')
  .select('precipitation');

print('========== CHIRPS DATA ==========');
print('CHIRPS daily count:', daily.size());

var years = ee.List.sequence(1981, 2025);
var months = ee.List.sequence(1, 12);

var monthlyList = years.map(function(year) {
  year = ee.Number(year);
  return months.map(function(month) {
    month = ee.Number(month);
    var start = ee.Date.fromYMD(year, month, 1);
    var end = start.advance(1, 'month');

    var image = daily
      .filterDate(start, end)
      .sum()
      .rename('monthly_rainfall');

    return image.set({
      'year': year,
      'month': month,
      'system:time_start': start.millis()
    });
  });
}).flatten();

var monthly = ee.ImageCollection.fromImages(monthlyList);

print('Monthly rainfall count:', monthly.size());


// THREE-MONTH RAINFALL ACCUMULATION
// ------------------------------------------------------------

function rainfall3Month(year, month) {
  year = ee.Number(year);
  month = ee.Number(month);

  var endDate = ee.Date.fromYMD(year, month, 1);
  var startDate = endDate.advance(-2, 'month');

  var image = monthly
    .filterDate(startDate, endDate.advance(1, 'month'))
    .sum()
    .rename('rainfall_3month');

  return image.set({
    'year': year,
    'month': month,
    'system:time_start': endDate.millis()
  });
}


// SPI-3 GAMMA CALCULATION FUNCTION
// ------------------------------------------------------------

function calculateSPI(year, month) {
  year = ee.Number(year);
  month = ee.Number(month);

  // Jan-Feb 1981 diabaikan karena akumulasi 3 bulanan belum lengkap (butuh data 1980)
  var firstYear = ee.Number(
    ee.Algorithms.If(month.lte(2), 1982, 1981)
  );

  var baselineYears = ee.List.sequence(firstYear, 2020);

  // Curah hujan historis 3-bulanan pada bulan kalender yang sama
  var baseline = ee.ImageCollection.fromImages(
    baselineYears.map(function(y) {
      return rainfall3Month(y, month);
    })
  );

  // Akumulasi curah hujan target
  var target = rainfall3Month(year, month);

  // Observasi hujan positif (> 0 mm)
  var positive = baseline.map(function(image) {
    return image.updateMask(image.gt(0));
  });

  var positiveCount = positive.count();
  var totalCount = baseline.count();

  var meanRainfall = positive.mean();

  var meanLog = positive.map(function(image) {
    return image.log();
  }).mean();

  // Estimasi parameter Gamma (Alpha & Beta) menggunakan Thom Approximation
  var A = meanRainfall.log().subtract(meanLog);

  var alpha = ee.Image(1)
    .add(
      ee.Image(1).add(A.multiply(4 / 3)).sqrt()
    )
    .divide(A.multiply(4));

  var beta = meanRainfall.divide(alpha);

  // Probabilitas hujan nol
  var zeroCount = baseline.map(function(image) {
    return image.eq(0);
  }).sum();

  var q = zeroCount.divide(totalCount);

  // Sintaks GEE gammainc(a, x) -> alpha.gammainc(standardizedTarget)
  var standardizedTarget = target.divide(beta);
  var gammaCDF = alpha.gammainc(standardizedTarget);

  // Probabilitas kumulatif campuran H(x)
  var probability = q.add(
    ee.Image(1).subtract(q).multiply(gammaCDF)
  );

  probability = probability.clamp(0.000001, 0.999999);

  // Inversi Transformasi Normal Standar (SPI Z-Score)
  var spi = probability
    .multiply(2)
    .subtract(1)
    .erfInv()
    .multiply(Math.sqrt(2));

  // Masking area valid
  var validMask = A.gt(0)
    .and(alpha.gt(0))
    .and(beta.gt(0))
    .and(positiveCount.gte(30))
    .and(totalCount.gte(30));

  spi = spi
    .rename('SPI3')
    .updateMask(validMask)
    .clip(region);

  return spi.set({
    'year': year,
    'month': month,
    'baseline_year_count': baselineYears.size(),
    'system:time_start': ee.Date.fromYMD(year, month, 1).millis()
  });
}


// CALCULATE SPI-3 JANUARY-DECEMBER 2025
// ------------------------------------------------------------

var spiList = months.map(function(month) {
  return calculateSPI(2025, month);
});

var spiImages = ee.ImageCollection.fromImages(spiList);

print('========== SPI-3 FULL YEAR 2025 ==========');
print('SPI-3 monthly image count:', spiImages.size());


// MAP VISUALIZATION (SEPTEMBER 2025)
// ------------------------------------------------------------

var septemberSPI = ee.Image(
  spiImages.filter(ee.Filter.eq('month', 9)).first()
);

var spiPalette = [
  '8c2d04', // Sangat Kering Ekstrem
  'd7301f', // Kering Berat
  'fc8d59', // Kering Sedang
  'fee08b', // Agak Kering
  'ffffbf', // Normal
  'd9ef8b', // Agak Basah
  '91cf60', // Basah Sedang
  '1a9850', // Basah Berat
  '006837'  // Sangat Basah
];

Map.addLayer(
  septemberSPI,
  {
    min: -3,
    max: 3,
    palette: spiPalette
  },
  'SPI-3 September 2025'
);


// CREATE MONTHLY STATISTICS FEATURE COLLECTION (2025)
// ------------------------------------------------------------

var statistics = ee.FeatureCollection(
  months.map(function(month) {
    month = ee.Number(month);

    var image = ee.Image(
      spiImages.filter(ee.Filter.eq('month', month)).first()
    );

    var stats = image.reduceRegion({
      reducer: ee.Reducer.minMax().combine({
        reducer2: ee.Reducer.mean(),
        sharedInputs: true
      }),
      geometry: region,
      scale: 5566,
      maxPixels: 1e13,
      tileScale: 4
    });

    return ee.Feature(null, stats).set({
      'year': 2025,
      'month': month,
      'baseline_year_count': image.get('baseline_year_count')
    });
  })
);

print('========== STATISTICS 2025 ==========');
print('SPI-3 statistics row count:', statistics.size());
print('First monthly statistics row:', statistics.first());


// CREATE 12-BAND GEOTIFF STACK
// ------------------------------------------------------------

var bandNames = months.map(function(month) {
  return ee.String('SPI3_2025_')
    .cat(ee.Number(month).format('%02d'));
});

var spiStack = spiImages
  .sort('month')
  .toBands()
  .rename(bandNames);

print('========== 12-BAND RASTER STACK ==========');
print('SPI-3 stack:', spiStack);
print('SPI-3 band names:', spiStack.bandNames());


// TIME-SERIES CHART (2025 MONTHLY TREND)
// ------------------------------------------------------------

var chart2025 = ui.Chart.image.series({
  imageCollection: spiImages,
  region: region,
  reducer: ee.Reducer.mean(),
  scale: 5566,
  xProperty: 'system:time_start'
}).setOptions({
  title: 'Trend SPI-3 Rata-Rata Lampung Sepanjang Tahun 2025',
  hAxis: {title: 'Bulan', format: 'MMM YYYY'},
  vAxis: {title: 'Nilai SPI-3'},
  colors: ['#d7301f'],
  lineWidth: 2,
  pointSize: 4
});

print(chart2025);


// EXPORT 12-BAND GEOTIFF & CSV STACK TO DRIVE
// ------------------------------------------------------------

Export.image.toDrive({
  image: spiStack,
  description: 'SPI3_Monthly_2025_Stack',
  folder: 'Agricultural_Drought_Lampung',
  fileNamePrefix: 'SPI3_Monthly_2025_Stack',
  region: region,
  scale: 5566,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

Export.table.toDrive({
  collection: statistics,
  description: 'SPI3_Monthly_Statistics_2025',
  folder: 'Agricultural_Drought_Lampung',
  fileNamePrefix: 'SPI3_Monthly_Statistics_2025',
  fileFormat: 'CSV'
});

print('========== EXPORT STATUS ==========');
print('Task ekspor 12-band GeoTIFF Stack dan CSV statistik 2025 berhasil dibuat!');