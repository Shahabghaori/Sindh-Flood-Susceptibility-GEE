// ============================================================
// SINDH FLOOD SUSCEPTIBILITY MAPPING
// Pakistan | Google Earth Engine
// Water and River Priority Weighted Overlay
// Classes: Very Low, Low, Medium, High, Very High
// ============================================================


// 1. STUDY AREA
var studyArea = ee.FeatureCollection(
  'projects/causal-galaxy-482405-n9/assets/Sindh'
);

var roi = studyArea.geometry();

Map.centerObject(studyArea, 7);

Map.addLayer(
  studyArea.style({
    color: '000000',
    fillColor: '00000000',
    width: 2
  }),
  {},
  'Sindh Boundary'
);


// 2. SRTM ELEVATION AND SLOPE
var dem = ee.Image('USGS/SRTMGL1_003')
  .select('elevation')
  .clip(roi);

var slope = ee.Terrain.slope(dem).clip(roi);

// Lower elevation = higher susceptibility
var elevationRisk = ee.Image(1)
  .where(dem.lte(80), 2)
  .where(dem.lte(40), 3)
  .where(dem.lte(20), 4)
  .where(dem.lte(10), 5)
  .rename('Elevation_Risk');

// Flat slope = higher susceptibility
var slopeRisk = ee.Image(1)
  .where(slope.lte(10), 2)
  .where(slope.lte(5), 3)
  .where(slope.lte(2), 4)
  .where(slope.lte(1), 5)
  .rename('Slope_Risk');


// 3. MERIT HYDRO - HAND AND RIVER NETWORK
var merit = ee.Image('MERIT/Hydro/v1_0_1');

var hand = merit.select('hnd')
  .clip(roi);

var riverWidth = merit.select('wth')
  .clip(roi);

// HAND risk
var handRisk = ee.Image(1)
  .where(hand.lte(20), 2)
  .where(hand.lte(10), 3)
  .where(hand.lte(5), 4)
  .where(hand.lte(2), 5)
  .rename('HAND_Risk');

// River network
var riverMask = riverWidth.gt(0)
  .unmask(0)
  .clip(roi);

// Distance from rivers in metres
var riverDistance = riverMask
  .fastDistanceTransform(1024)
  .sqrt()
  .multiply(riverWidth.projection().nominalScale())
  .clip(roi);

// River proximity risk
var riverRisk = ee.Image(1)
  .where(riverDistance.lte(10000), 2)
  .where(riverDistance.lte(5000), 3)
  .where(riverDistance.lte(2000), 4)
  .where(riverDistance.lte(1000), 5)
  .rename('River_Risk');


// 4. JRC HISTORICAL WATER BODIES
var gsw = ee.Image('JRC/GSW1_4/GlobalSurfaceWater');

var occurrence = gsw.select('occurrence')
  .unmask(0)
  .clip(roi);

// Permanent water bodies
var permanentWater = occurrence.gte(90);

// Seasonal and intermittent water
var seasonalWater = occurrence.gte(10)
  .and(occurrence.lt(90));

// All historically observed water
var allWater = occurrence.gte(10);

// Historical water occurrence risk
var waterOccurrenceRisk = ee.Image(1)
  .where(occurrence.gte(10), 3)
  .where(occurrence.gte(50), 4)
  .where(occurrence.gte(90), 5)
  .rename('Water_Occurrence_Risk');

// Distance from historical water bodies
var waterDistance = allWater
  .unmask(0)
  .fastDistanceTransform(1024)
  .sqrt()
  .multiply(occurrence.projection().nominalScale())
  .clip(roi);

// Water body proximity risk
var waterProximityRisk = ee.Image(1)
  .where(waterDistance.lte(10000), 2)
  .where(waterDistance.lte(5000), 3)
  .where(waterDistance.lte(2000), 4)
  .where(waterDistance.lte(500), 5)
  .rename('Water_Proximity_Risk');


// 5. CHIRPS MONSOON RAINFALL
// June to September, 2001-2024

var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .select('precipitation');

var years = ee.List.sequence(2001, 2024);

var annualRainfall = ee.ImageCollection.fromImages(
  years.map(function(year) {
    year = ee.Number(year);

    return chirps
      .filter(ee.Filter.calendarRange(year, year, 'year'))
      .filter(ee.Filter.calendarRange(6, 9, 'month'))
      .sum()
      .set('year', year);
  })
);

var rainfall = annualRainfall
  .mean()
  .clip(roi)
  .rename('Monsoon_Rainfall');

// Calculate rainfall percentiles
var rainStats = rainfall.reduceRegion({
  reducer: ee.Reducer.percentile([20, 40, 60, 80]),
  geometry: roi,
  scale: 5500,
  maxPixels: 1e13,
  bestEffort: true,
  tileScale: 4
});

print('Rainfall thresholds', rainStats);

var p20 = ee.Number(rainStats.get('Monsoon_Rainfall_p20'));
var p40 = ee.Number(rainStats.get('Monsoon_Rainfall_p40'));
var p60 = ee.Number(rainStats.get('Monsoon_Rainfall_p60'));
var p80 = ee.Number(rainStats.get('Monsoon_Rainfall_p80'));

var rainfallRisk = ee.Image(1)
  .where(rainfall.gt(p20), 2)
  .where(rainfall.gt(p40), 3)
  .where(rainfall.gt(p60), 4)
  .where(rainfall.gt(p80), 5)
  .rename('Rainfall_Risk');


// 6. ESA WORLD COVER
var landcover = ee.ImageCollection('ESA/WorldCover/v200')
  .first()
  .select('Map')
  .clip(roi);

// LULC susceptibility scores
var lulcRisk = landcover.remap(
  [10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100],
  [ 1,  2,  2,  4,  4,  3,  1,  5,  5,  5,   1]
).rename('LULC_Risk');


// 7. WATER AND RIVER PRIORITY WEIGHTED OVERLAY
//
// Factor                       Weight
// River proximity               25%
// Historical water occurrence   20%
// HAND                           15%
// Elevation                      10%
// Monsoon rainfall               10%
// Slope                           5%
// LULC                            5%
// Water body proximity           10%
// Total                          100%

var floodScore = riverRisk.multiply(0.25)
  .add(waterOccurrenceRisk.multiply(0.20))
  .add(handRisk.multiply(0.15))
  .add(elevationRisk.multiply(0.10))
  .add(rainfallRisk.multiply(0.10))
  .add(slopeRisk.multiply(0.05))
  .add(lulcRisk.multiply(0.05))
  .add(waterProximityRisk.multiply(0.10))
  .rename('Flood_Susceptibility_Score');


// 8. FIVE FLOOD SUSCEPTIBILITY CLASSES

// 1 = Very Low
// 2 = Low
// 3 = Medium
// 4 = High
// 5 = Very High

var floodClass = ee.Image(1)
  .where(floodScore.gte(1.8), 2)
  .where(floodScore.gte(2.6), 3)
  .where(floodScore.gte(3.4), 4)
  .where(floodScore.gte(4.2), 5)
  .rename('Flood_Risk')
  .toByte()
  .clip(roi);


// 9. MAP COLORS

var riskPalette = [
  '006400',  // Very Low - Dark Green
  '7FFF00',  // Low - Light Green
  'FFFF00',  // Medium - Yellow
  'FF8C00',  // High - Orange
  'FF0000'   // Very High - Red
];

Map.addLayer(
  floodClass,
  {
    min: 1,
    max: 5,
    palette: riskPalette
  },
  'FINAL SINDH FLOOD SUSCEPTIBILITY MAP',
  true
);


// 10. WATER BODY AND RIVER MAPS

Map.addLayer(
  riverMask.selfMask(),
  {palette: ['000080']},
  'River and Drainage Network',
  true
);

Map.addLayer(
  permanentWater.selfMask(),
  {palette: ['0000FF']},
  'Permanent Water Bodies',
  false
);

Map.addLayer(
  seasonalWater.selfMask(),
  {palette: ['00FFFF']},
  'Seasonal Water Bodies',
  false
);

Map.addLayer(
  waterDistance,
  {
    min: 0,
    max: 10000,
    palette: ['0000FF', '00FFFF', 'FFFF00', 'FF0000']
  },
  'Distance from Historical Water',
  false
);


// 11. OTHER FACTOR MAPS

Map.addLayer(
  hand,
  {
    min: 0,
    max: 30,
    palette: ['red', 'orange', 'yellow', 'blue']
  },
  'HAND',
  false
);

Map.addLayer(
  dem,
  {
    min: 0,
    max: 1500,
    palette: ['blue', 'green', 'yellow', 'brown']
  },
  'Elevation',
  false
);

Map.addLayer(
  slope,
  {
    min: 0,
    max: 30,
    palette: ['blue', 'yellow', 'red']
  },
  'Slope',
  false
);

Map.addLayer(
  rainfall,
  {
    min: 0,
    max: 400,
    palette: ['white', 'yellow', 'orange', 'red']
  },
  'Mean Monsoon Rainfall',
  false
);

Map.addLayer(
  landcover,
  {
    min: 10,
    max: 100,
    palette: [
      '006400', 'FFBB22', 'FFFF4C',
      'F096FF', 'FA0000', 'B4B4B4',
      'F0F0F0', '0064C8', '0096A0',
      '00CF75', 'FAE6A0'
    ]
  },
  'ESA Land Cover',
  false
);


// 12. AREA CALCULATION IN SQUARE KILOMETRES

var areaImage = ee.Image.pixelArea()
  .divide(1000000)
  .rename('Area')
  .addBands(floodClass);

var areaStats = areaImage.reduceRegion({
  reducer: ee.Reducer.sum().group({
    groupField: 1,
    groupName: 'Risk_Class'
  }),
  geometry: roi,
  scale: 90,
  maxPixels: 1e13,
  tileScale: 4
});

print('Flood Risk Area Statistics (km2)', areaStats);


// 13. RISK CLASS AREA TABLE

var classNames = ee.Dictionary({
  '1': 'Very Low',
  '2': 'Low',
  '3': 'Medium',
  '4': 'High',
  '5': 'Very High'
});

var groups = ee.List(areaStats.get('groups'));

var areaTable = ee.FeatureCollection(
  groups.map(function(item) {
    item = ee.Dictionary(item);

    var id = ee.Number(item.get('Risk_Class')).toInt();

    return ee.Feature(null, {
      Risk_Class: id,
      Risk_Level: classNames.get(id.format()),
      Area_km2: item.get('sum')
    });
  })
);

print('FINAL RISK CLASSIFICATION TABLE', areaTable);


// 14. EXPORT FLOOD SUSCEPTIBILITY MAP

Export.image.toDrive({
  image: floodClass,
  description: 'Sindh_Flood_Susceptibility_Map',
  folder: 'GEE_Flood_Maps',
  fileNamePrefix: 'Sindh_Flood_Risk_5_Classes',
  region: roi,
  scale: 90,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});


// 15. EXPORT CONTINUOUS FLOOD SCORE

Export.image.toDrive({
  image: floodScore,
  description: 'Sindh_Flood_Susceptibility_Score',
  folder: 'GEE_Flood_Maps',
  fileNamePrefix: 'Sindh_Flood_Continuous_Score',
  region: roi,
  scale: 90,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});


// 16. EXPORT WATER NETWORK INVENTORY

var waterInventory = ee.Image(0)
  .where(riverMask.eq(1), 1)
  .where(seasonalWater, 2)
  .where(permanentWater, 3)
  .rename('Water_Features')
  .toByte()
  .clip(roi);

Export.image.toDrive({
  image: waterInventory,
  description: 'Sindh_River_Water_Inventory',
  folder: 'GEE_Flood_Maps',
  fileNamePrefix: 'Sindh_Rivers_Water_Bodies',
  region: roi,
  scale: 90,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});


// 17. EXPORT AREA STATISTICS

Export.table.toDrive({
  collection: areaTable,
  description: 'Sindh_Flood_Risk_Area_Statistics',
  folder: 'GEE_Flood_Maps',
  fileFormat: 'CSV'
});

print('======================================');
print('SINDH FLOOD MAPPING COMPLETED');
print('Five risk classes generated.');
print('Water and river priority included.');
print('Check Tasks tab for exports.');
print('======================================');
