// ============================================================================
// GEE EXTRACTION SCRIPT FOR BTEH (ELEPHANT HABITAT ANALYSIS)
// ============================================================================
// Extracts: Hansen Forest Cover, GEDI Canopy Height, Biomass, and Vegetation
// Output: 6 GeoTIFFs at 10m resolution
// 
// Period 1 (Dec 31): 2021, 2022, 2023 (December monthly mean)
// Period 2 (Jan 1): 2024, 2025 (January monthly mean)
// 
// USAGE: Copy this entire script into Google Earth Engine Code Editor
// URL: https://code.earthengine.google.com/
// ============================================================================

// ============================================================================
// 1. CONFIGURATION - UPDATE THESE VALUES
// ============================================================================

var STUDY_AREA_ASSET = 'projects/your-project-id/assets/StudyArea'; // UPDATE THIS!
var EXPORT_SCALE = 10; // 10 meter resolution
var EXPORT_TO_DRIVE = true; // true = Google Drive, false = Cloud Storage
var EXPORT_BUCKET = ''; // Leave empty if using Drive, else set bucket name

// ============================================================================
// 2. LOAD STUDY AREA
// ============================================================================

var studyArea = ee.FeatureCollection(STUDY_AREA_ASSET).geometry();
var bounds = studyArea.bounds();

// Center map and visualize study area
Map.centerObject(bounds, 8);
Map.addLayer(bounds, {color: 'FF0000', weight: 2}, 'Study Area Boundary');

print('Study Area loaded successfully');
print('Geometry bounds:', bounds.getInfo());

// ============================================================================
// 3. HANSEN GLOBAL FOREST CHANGE - DECEMBER 31 (2021, 2022, 2023)
// ============================================================================

print('\n=== EXTRACTING HANSEN FOREST COVER ===');

var hansen = ee.Image('UMD/hansen/global_forest_change_2024_v1_12');
var treecover2000 = hansen.select('treecover2000');
var loss = hansen.select('loss');
var lossyear = hansen.select('lossyear');

// Extract for December 31 of each year
var hansenYears = [2021, 2022, 2023];

hansenYears.forEach(function(year) {
  // Calculate loss year threshold
  // Loss year encoding: 1=2001, 2=2002, ..., 21=2021, 22=2022, 23=2023
  var lossYearThreshold = year - 2000;
  
  // Create cumulative loss mask up to end of year
  var lossMask = lossyear.lte(lossYearThreshold).and(loss.eq(1));
  
  // Calculate current forest cover
  // Formula: forest_cover = treecover2000 - (loss_mask * 100)
  var currentForest = treecover2000.subtract(lossMask.multiply(100));
  
  // Ensure values stay in valid range (0-100%)
  currentForest = currentForest.max(0).min(100);
  
  // Clip to study area
  var clipped = currentForest.clip(studyArea);
  
  // Add to map for visualization
  var visParams = {min: 0, max: 100, palette: ['brown', 'yellow', 'green']};
  Map.addLayer(clipped, visParams, 'Hansen Forest Cover - Dec 31 ' + year, false);
  
  // Export to GeoTIFF
  var description = 'HANSEN_ForestCover_Dec31_' + year;
  
  if (EXPORT_TO_DRIVE) {
    Export.image.toDrive({
      image: clipped.setDefaultProjection('EPSG:4326', null, EXPORT_SCALE),
      description: description,
      scale: EXPORT_SCALE,
      region: studyArea,
      fileFormat: 'GeoTIFF',
      maxPixels: 1e13,
      formatOptions: {
        cloudOptimized: true
      }
    });
  } else {
    Export.image.toCloudStorage({
      image: clipped.setDefaultProjection('EPSG:4326', null, EXPORT_SCALE),
      description: description,
      bucket: EXPORT_BUCKET,
      fileNamePrefix: 'bteh/' + description,
      scale: EXPORT_SCALE,
      region: studyArea,
      fileFormat: 'GeoTIFF',
      maxPixels: 1e13,
      formatOptions: {
        cloudOptimized: true
      }
    });
  }
  
  print('✓ Queued: ' + description);
});

// ============================================================================
// 4. GEDI L2A CANOPY HEIGHT - DECEMBER (2021, 2022, 2023)
// ============================================================================

print('\n=== EXTRACTING GEDI CANOPY HEIGHT - DECEMBER ===');

var gediCanopyCollection = ee.ImageCollection('LARSE/GEDI/GEDI02_A_002_MONTHLY')
  .filterBounds(studyArea);

print('Total GEDI L2A images available:', gediCanopyCollection.size().getInfo());

var decemberPeriods = [
  {year: 2021, start: '2021-12-01', end: '2021-12-31'},
  {year: 2022, start: '2022-12-01', end: '2022-12-31'},
  {year: 2023, start: '2023-12-01', end: '2023-12-31'}
];

decemberPeriods.forEach(function(period) {
  var decemberData = gediCanopyCollection
    .filterDate(period.start, period.end)
    .select('rh98');  // RH98 = Canopy top height (98th percentile)
  
  var count = decemberData.size().getInfo();
  
  if (count > 0) {
    // Calculate mean canopy height for December
    var canopyHeight = decemberData.mean();
    var clipped = canopyHeight.clip(studyArea);
    
    // Add to map
    var visParams = {min: 0, max: 50, palette: ['white', 'yellow', 'red']};
    Map.addLayer(clipped, visParams, 'GEDI Canopy Height - Dec ' + period.year, false);
    
    // Export
    var description = 'GEDI_CanopyHeight_Dec_' + period.year;
    
    if (EXPORT_TO_DRIVE) {
      Export.image.toDrive({
        image: clipped.setDefaultProjection('EPSG:4326', null, EXPORT_SCALE),
        description: description,
        scale: EXPORT_SCALE,
        region: studyArea,
        fileFormat: 'GeoTIFF',
        maxPixels: 1e13,
        formatOptions: {
          cloudOptimized: true
        }
      });
    } else {
      Export.image.toCloudStorage({
        image: clipped.setDefaultProjection('EPSG:4326', null, EXPORT_SCALE),
        description: description,
        bucket: EXPORT_BUCKET,
        fileNamePrefix: 'bteh/' + description,
        scale: EXPORT_SCALE,
        region: studyArea,
        fileFormat: 'GeoTIFF',
        maxPixels: 1e13,
        formatOptions: {
          cloudOptimized: true
        }
      });
    }
    
    print('✓ Queued: ' + description + ' (' + count + ' images)');
  } else {
    print('⚠ No GEDI data for December ' + period.year);
  }
});

// ============================================================================
// 5. GEDI L2A CANOPY HEIGHT - JANUARY (2024, 2025)
// ============================================================================

print('\n=== EXTRACTING GEDI CANOPY HEIGHT - JANUARY ===');

var januaryCanopyPeriods = [
  {year: 2024, start: '2024-01-01', end: '2024-01-31'},
  {year: 2025, start: '2025-01-01', end: '2025-01-31'}
];

januaryCanopyPeriods.forEach(function(period) {
  var januaryData = gediCanopyCollection
    .filterDate(period.start, period.end)
    .select('rh98');
  
  var count = januaryData.size().getInfo();
  
  if (count > 0) {
    var canopyHeight = januaryData.mean();
    var clipped = canopyHeight.clip(studyArea);
    
    var visParams = {min: 0, max: 50, palette: ['white', 'yellow', 'red']};
    Map.addLayer(clipped, visParams, 'GEDI Canopy Height - Jan ' + period.year, false);
    
    var description = 'GEDI_CanopyHeight_Jan_' + period.year;
    
    if (EXPORT_TO_DRIVE) {
      Export.image.toDrive({
        image: clipped.setDefaultProjection('EPSG:4326', null, EXPORT_SCALE),
        description: description,
        scale: EXPORT_SCALE,
        region: studyArea,
        fileFormat: 'GeoTIFF',
        maxPixels: 1e13,
        formatOptions: {
          cloudOptimized: true
        }
      });
    } else {
      Export.image.toCloudStorage({
        image: clipped.setDefaultProjection('EPSG:4326', null, EXPORT_SCALE),
        description: description,
        bucket: EXPORT_BUCKET,
        fileNamePrefix: 'bteh/' + description,
        scale: EXPORT_SCALE,
        region: studyArea,
        fileFormat: 'GeoTIFF',
        maxPixels: 1e13,
        formatOptions: {
          cloudOptimized: true
        }
      });
    }
    
    print('✓ Queued: ' + description + ' (' + count + ' images)');
  } else {
    print('⚠ No GEDI data for January ' + period.year);
  }
});

// ============================================================================
// 6. GEDI L4A BIOMASS DENSITY - DECEMBER (2021, 2022, 2023)
// ============================================================================

print('\n=== EXTRACTING GEDI BIOMASS - DECEMBER ===');

var gediBiomassCollection = ee.ImageCollection('LARSE/GEDI/GEDI04_A_002_MONTHLY')
  .filterBounds(studyArea);

print('Total GEDI L4A images available:', gediBiomassCollection.size().getInfo());

var decemberBiomassPeriods = [
  {year: 2021, start: '2021-12-01', end: '2021-12-31'},
  {year: 2022, start: '2022-12-01', end: '2022-12-31'},
  {year: 2023, start: '2023-12-01', end: '2023-12-31'}
];

decemberBiomassPeriods.forEach(function(period) {
  var decemberData = gediBiomassCollection
    .filterDate(period.start, period.end)
    .select('agbd');  // AGBD = Aboveground Biomass Density (Mg/ha)
  
  var count = decemberData.size().getInfo();
  
  if (count > 0) {
    // Calculate mean biomass for December
    var biomass = decemberData.mean();
    var clipped = biomass.clip(studyArea);
    
    // Add to map
    var visParams = {min: 0, max: 300, palette: ['blue', 'green', 'brown']};
    Map.addLayer(clipped, visParams, 'GEDI Biomass - Dec ' + period.year, false);
    
    // Export
    var description = 'GEDI_Biomass_Dec_' + period.year;
    
    if (EXPORT_TO_DRIVE) {
      Export.image.toDrive({
        image: clipped.setDefaultProjection('EPSG:4326', null, EXPORT_SCALE),
        description: description,
        scale: EXPORT_SCALE,
        region: studyArea,
        fileFormat: 'GeoTIFF',
        maxPixels: 1e13,
        formatOptions: {
          cloudOptimized: true
        }
      });
    } else {
      Export.image.toCloudStorage({
        image: clipped.setDefaultProjection('EPSG:4326', null, EXPORT_SCALE),
        description: description,
        bucket: EXPORT_BUCKET,
        fileNamePrefix: 'bteh/' + description,
        scale: EXPORT_SCALE,
        region: studyArea,
        fileFormat: 'GeoTIFF',
        maxPixels: 1e13,
        formatOptions: {
          cloudOptimized: true
        }
      });
    }
    
    print('✓ Queued: ' + description + ' (' + count + ' images)');
  } else {
    print('⚠ No GEDI biomass data for December ' + period.year);
  }
});

// ============================================================================
// 7. GEDI L4A BIOMASS DENSITY - JANUARY (2024, 2025)
// ============================================================================

print('\n=== EXTRACTING GEDI BIOMASS - JANUARY ===');

var januaryBiomassPeriods = [
  {year: 2024, start: '2024-01-01', end: '2024-01-31'},
  {year: 2025, start: '2025-01-01', end: '2025-01-31'}
];

januaryBiomassPeriods.forEach(function(period) {
  var januaryData = gediBiomassCollection
    .filterDate(period.start, period.end)
    .select('agbd');
  
  var count = januaryData.size().getInfo();
  
  if (count > 0) {
    var biomass = januaryData.mean();
    var clipped = biomass.clip(studyArea);
    
    var visParams = {min: 0, max: 300, palette: ['blue', 'green', 'brown']};
    Map.addLayer(clipped, visParams, 'GEDI Biomass - Jan ' + period.year, false);
    
    var description = 'GEDI_Biomass_Jan_' + period.year;
    
    if (EXPORT_TO_DRIVE) {
      Export.image.toDrive({
        image: clipped.setDefaultProjection('EPSG:4326', null, EXPORT_SCALE),
        description: description,
        scale: EXPORT_SCALE,
        region: studyArea,
        fileFormat: 'GeoTIFF',
        maxPixels: 1e13,
        formatOptions: {
          cloudOptimized: true
        }
      });
    } else {
      Export.image.toCloudStorage({
        image: clipped.setDefaultProjection('EPSG:4326', null, EXPORT_SCALE),
        description: description,
        bucket: EXPORT_BUCKET,
        fileNamePrefix: 'bteh/' + description,
        scale: EXPORT_SCALE,
        region: studyArea,
        fileFormat: 'GeoTIFF',
        maxPixels: 1e13,
        formatOptions: {
          cloudOptimized: true
        }
      });
    }
    
    print('✓ Queued: ' + description + ' (' + count + ' images)');
  } else {
    print('⚠ No GEDI data for January ' + period.year);
  }
});

// ============================================================================
// 8. OPTIONAL: GEDI GRIDDED VEGETATION METRICS (Full Mission: Apr 2019-Mar 2023)
// ============================================================================

print('\n=== OPTIONAL: GEDI VEGETATION STRUCTURE (Full Mission) ===');

var gediVeg = ee.Image('LARSE/GEDI/GRIDDEDVEG_002_V1/1KM/gediv002_rh98-mean_vf_20190417_20230316');
var vegStructure = gediVeg.select('mean').clip(studyArea);

var vegVisParams = {min: 0, max: 50, palette: ['white', 'yellow', 'orange', 'red']};
Map.addLayer(vegStructure, vegVisParams, 'GEDI Vegetation Structure (Full Mission)', false);

print('✓ GEDI Vegetation Structure added to map (Full mission: Apr 2019-Mar 2023)');

// ============================================================================
// 9. SUMMARY
// ============================================================================

print('\n' + '='.repeat(70));
print('✅ GEE EXTRACTION SCRIPT COMPLETE');
print('='.repeat(70));
print('\nSUMMARY OF EXPORTS:');
print('  - HANSEN_ForestCover_Dec31_2021');
print('  - HANSEN_ForestCover_Dec31_2022');
print('  - HANSEN_ForestCover_Dec31_2023');
print('  - GEDI_CanopyHeight_Dec_2021');
print('  - GEDI_CanopyHeight_Dec_2022');
print('  - GEDI_CanopyHeight_Dec_2023');
print('  - GEDI_CanopyHeight_Jan_2024');
print('  - GEDI_CanopyHeight_Jan_2025');
print('  - GEDI_Biomass_Dec_2021');
print('  - GEDI_Biomass_Dec_2022');
print('  - GEDI_Biomass_Dec_2023');
print('  - GEDI_Biomass_Jan_2024');
print('  - GEDI_Biomass_Jan_2025');
print('\nTotal exports queued: 13 GeoTIFFs at 10m resolution');
print('\nNEXT STEPS:');
print('  1. Go to Tasks tab (right panel)');
print('  2. Click "Run" on each export task to start download');
print('  3. GeoTIFFs will download to Google Drive /Earth Engine Exports/');
print('  4. Wait 10-60 minutes per file depending on size');
print('  5. Download to local machine and load into R/QGIS/ArcGIS');
print('\nAll outputs: Cloud-Optimized GeoTIFF, EPSG:4326, 10m resolution');
print('='.repeat(70));
