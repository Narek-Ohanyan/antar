"""Earth Engine exports: yearly vitality (kNDVI) composites, LandTrendr segmentation,
and disturbance ancillary layers, all on the master grid (EPSG:32638, 30 m).

Two products, matching Table 3/4 and Sec. 7.1 of the concept note:

* :func:`export_vitality_composites` -- yearly growing-season (Jul-Aug) median kNDVI
  (Camps-Valls et al. 2021) from HLS v2.0 (2013-) merged with Landsat Collection 2
  Surface Reflectance (2000-, where HLS is unavailable), with a per-pixel valid-
  observation-count band, plus LandTrendr (Kennedy et al. 2010) temporal segmentation
  of the resulting annual series. This is the input :func:`antar.hazard.observation.
  dieback_event` is written against; that function still owns the threshold rule.
* :func:`export_disturbance_ancillary` -- Hansen Global Forest Change loss-year and
  MODIS MCD64A1 burned-area, combined into a per-year "no disturbance" boolean layer
  -- exactly the ``no_disturbance`` input :func:`dieback_event` expects, so that
  dieback (a climate response) is never confused with harvest or fire.

Six more, covering the rest of Table 4's realistically-obtainable rows (everything not
gated behind a national archive, field survey, or specialist trait database -- those are
documented gaps, not silently substituted):

* :func:`export_terrain` -- Copernicus GLO-30 DEM, slope, aspect (TWI/TPI/cold-air-pooling
  index are left to local post-processing: they need neighbourhood operations GEE's
  per-pixel terrain functions do not provide directly).
* :func:`export_soils` -- SoilGrids 2.0 clay/sand/silt fraction and organic carbon, for the
  Saxton & Rawls (2006) pedotransfer functions TOPOHYDRO's soil bucket needs.
* :func:`export_era5land_forcing` -- wind speed, shortwave/longwave radiation, dewpoint
  temperature and surface pressure from ERA5-Land, exactly the variables Table 4 sources
  from ERA5-Land rather than CHELSA.
* :func:`export_vegetation_state` -- MODIS LAI/FPAR, land-surface phenology (start/end of
  season), and ESA WorldCover tree-cover fraction (complementing Hansen GFC).
* :func:`export_snow` -- MODIS daily snow-cover fraction.
* :func:`export_land_tenure` -- WDPA protected-area boundaries, rasterised to the master
  grid -- feeds the eligibility mask directly (``configs/study_area.yaml``'s
  ``eligibility_mask_exclusions``).

Not pulled, and not fakeable: national forest inventory plots, provenance/genetic trial
locations, insect/pathogen outbreak records, treeline field-survey transects (all
institution- or field-survey-only), road network/accessibility (no public Earth Engine
asset found for this region under any plausible id -- checked, not assumed), and trait
databases (XFT, TRY -- specialist data services, not Earth Engine assets).

All export functions submit asynchronous Earth Engine batch tasks to Google Drive
and return the started :class:`ee.batch.Task` object(s) -- they do not block until
the exports finish (that can take hours for a study-area-wide, multi-decade pull).
Poll ``task.status()`` to check progress.
"""
from __future__ import annotations

import ee


# ---------------------------------------------------------------------------
# kNDVI and per-sensor cloud/shadow/snow masking
# ---------------------------------------------------------------------------

def _kndvi(nir, red):
    """kNDVI = tanh(NDVI^2) (Camps-Valls et al. 2021).

    The paper's kernel form is tanh(((NIR-Red)/(2*sigma))^2); with the
    per-pixel adaptive sigma = 0.5*(NIR+Red) it recommends, (NIR-Red)/(2*sigma)
    is exactly NDVI, so this *is* the kernel formula, not an approximation.
    """
    ndvi = nir.subtract(red).divide(nir.add(red))
    return ndvi.pow(2).tanh().rename("kndvi")


def _hls_kndvi(image, nir_band: str, red_band: str):
    """One HLS v2.0 image -> masked kNDVI. Fmask bits (NASA HLS v2.0 User Guide):
    1=cloud, 2=adjacent-to-cloud/shadow, 3=cloud shadow, 4=snow/ice.
    """
    fmask = image.select("Fmask").toInt()
    bad = (
        fmask.bitwiseAnd(1 << 1).neq(0)
        .Or(fmask.bitwiseAnd(1 << 2).neq(0))
        .Or(fmask.bitwiseAnd(1 << 3).neq(0))
        .Or(fmask.bitwiseAnd(1 << 4).neq(0))
    )
    nir = image.select(nir_band).toFloat().multiply(0.0001)
    red = image.select(red_band).toFloat().multiply(0.0001)
    return _kndvi(nir, red).updateMask(bad.Not())


def _landsat_c2sr_kndvi(image, nir_band: str, red_band: str):
    """One Landsat Collection 2 Level-2 SR image -> masked kNDVI. QA_PIXEL bits
    (USGS Landsat Collection 2 Level-2 QA_PIXEL): 1=dilated cloud, 2=cirrus,
    3=cloud, 4=cloud shadow, 5=snow. SR scaling: reflectance = DN*0.0000275 - 0.2
    (the additive offset does not cancel in the NDVI ratio, unlike HLS's pure
    multiplicative scale, so it must be applied before computing kNDVI).
    """
    qa = image.select("QA_PIXEL").toInt()
    bad = (
        qa.bitwiseAnd(1 << 1).neq(0)
        .Or(qa.bitwiseAnd(1 << 2).neq(0))
        .Or(qa.bitwiseAnd(1 << 3).neq(0))
        .Or(qa.bitwiseAnd(1 << 4).neq(0))
        .Or(qa.bitwiseAnd(1 << 5).neq(0))
    )
    nir = image.select(nir_band).toFloat().multiply(0.0000275).add(-0.2)
    red = image.select(red_band).toFloat().multiply(0.0000275).add(-0.2)
    return _kndvi(nir, red).updateMask(bad.Not())


def _merged_kndvi_collection(aoi, date_start: str, date_end: str):
    """kNDVI ImageCollection over [date_start, date_end), merging every sensor whose
    archive actually covers that window (HLS v2.0, Landsat 5/7/8/9 Collection 2 SR).
    Each output image carries a single 'kndvi' band and 'system:time_start'.
    """
    def tag(img, kndvi_img):
        return kndvi_img.copyProperties(img, ["system:time_start"])

    hlsl = ee.ImageCollection("NASA/HLS/HLSL30/v002").filterBounds(aoi).filterDate(date_start, date_end)
    hlss = ee.ImageCollection("NASA/HLS/HLSS30/v002").filterBounds(aoi).filterDate(date_start, date_end)
    l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2").filterBounds(aoi).filterDate(date_start, date_end)
    l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2").filterBounds(aoi).filterDate(date_start, date_end)
    l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2").filterBounds(aoi).filterDate(date_start, date_end)
    l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2").filterBounds(aoi).filterDate(date_start, date_end)

    parts = [
        hlsl.map(lambda img: tag(img, _hls_kndvi(img, "B5", "B4"))),
        hlss.map(lambda img: tag(img, _hls_kndvi(img, "B8", "B4"))),
        l5.map(lambda img: tag(img, _landsat_c2sr_kndvi(img, "SR_B4", "SR_B3"))),
        l7.map(lambda img: tag(img, _landsat_c2sr_kndvi(img, "SR_B4", "SR_B3"))),
        l8.map(lambda img: tag(img, _landsat_c2sr_kndvi(img, "SR_B5", "SR_B4"))),
        l9.map(lambda img: tag(img, _landsat_c2sr_kndvi(img, "SR_B5", "SR_B4"))),
    ]
    out = parts[0]
    for p in parts[1:]:
        out = out.merge(p)
    return out


# ---------------------------------------------------------------------------
# Yearly growing-season composites + LandTrendr
# ---------------------------------------------------------------------------

_LANDTRENDR_DEFAULTS = dict(
    maxSegments=6,
    spikeThreshold=0.9,
    vertexCountOvershoot=3,
    preventOneYearRecovery=False,
    recoveryThreshold=0.25,
    pvalThreshold=0.05,
    bestModelProportion=0.75,
    minObservationsNeeded=6,
)  # Kennedy et al. (2010); the standard published LT-GEE defaults, not tuned here.


def export_vitality_composites(
    bbox_wgs84,
    year_start: int,
    year_end: int,
    crs: str,
    scale_m: float,
    drive_folder: str = "antar_gee_exports",
    growing_season_months: tuple[int, int] = (7, 8),
    landtrendr_params: dict | None = None,
):
    """Yearly growing-season median kNDVI + valid-observation count, and LandTrendr
    segmentation of the resulting series, both exported to Google Drive.

    Returns ``(composite_task, landtrendr_task)`` -- both already started
    (``ee.batch.Task.start()`` called); poll ``task.status()`` for progress.
    """
    aoi = ee.Geometry.Rectangle(list(bbox_wgs84))
    m0, m1 = growing_season_months
    years = list(range(year_start, year_end + 1))

    composite_bands = {}
    lt_input_images = []
    for y in years:
        coll = _merged_kndvi_collection(aoi, f"{y}-{m0:02d}-01", f"{y}-{m1:02d}-{31 if m1 in (1,3,5,7,8,10,12) else 30}")
        median = coll.select("kndvi").median().rename(f"kndvi_{y}")
        count = coll.select("kndvi").count().rename(f"valid_count_{y}")
        composite_bands[f"kndvi_{y}"] = median
        composite_bands[f"valid_count_{y}"] = count
        yearly = coll.select("kndvi").median().rename("kndvi").set("system:time_start", ee.Date.fromYMD(y, m0, 1).millis())
        lt_input_images.append(yearly)

    composite_image = ee.Image.cat(list(composite_bands.values())).toFloat()

    lt = ee.Algorithms.TemporalSegmentation.LandTrendr(
        timeSeries=ee.ImageCollection(lt_input_images), **(landtrendr_params or _LANDTRENDR_DEFAULTS)
    )
    arr = lt.select("LandTrendr")
    # LandTrendr's per-pixel output array drops years where that pixel's composite was
    # masked (no valid satellite observation that Jul-Aug), so its length along axis 1
    # varies across the raster -- arrayFlatten needs a constant length. Pad every pixel
    # up to len(years) with a sentinel (-9999, distinguishable from any real kNDVI/year/
    # vertex value) so the export has a well-defined, constant band count; the sentinel
    # marks "this pixel had fewer valid years than the record," not a real observation.
    arr = arr.arrayPad([4, len(years)], -9999)
    year_labels = [str(y) for y in years]
    lt_flat = arr.arrayFlatten([["year", "original", "fitted", "vertex"], year_labels])

    composite_task = ee.batch.Export.image.toDrive(
        image=composite_image, description="antar_vitality_composites", folder=drive_folder,
        region=aoi, crs=crs, scale=scale_m, maxPixels=1e13,
    )
    landtrendr_task = ee.batch.Export.image.toDrive(
        image=lt_flat.toFloat(), description="antar_landtrendr_segmentation", folder=drive_folder,
        region=aoi, crs=crs, scale=scale_m, maxPixels=1e13,
    )
    composite_task.start()
    landtrendr_task.start()
    return composite_task, landtrendr_task


# ---------------------------------------------------------------------------
# Disturbance ancillary: Hansen GFC loss-year + MODIS burned area
# ---------------------------------------------------------------------------

def export_disturbance_ancillary(
    bbox_wgs84,
    year_start: int,
    year_end: int,
    crs: str,
    scale_m: float,
    drive_folder: str = "antar_gee_exports",
    hansen_asset: str = "UMD/hansen/global_forest_change_2023_v1_11",
):
    """Per-year 'no_disturbance' boolean layer: NOT(Hansen GFC loss that year OR any
    MODIS MCD64A1 burn that year). This is exactly the ``no_disturbance`` input
    :func:`antar.hazard.observation.dieback_event` expects -- computed here, never
    inferred from the vitality signal itself, so that a real dieback event is never
    misread as (or masked by) a harvest/fire event and vice versa.
    """
    aoi = ee.Geometry.Rectangle(list(bbox_wgs84))
    gfc = ee.Image(hansen_asset)
    lossyear = gfc.select("lossyear")  # 0 = no loss; 1..N = loss in (2000 + value)

    bands = {}
    for y in range(year_start, year_end + 1):
        loss_y = lossyear.eq(y - 2000) if y > 2000 else ee.Image(0)
        burned = (
            ee.ImageCollection("MODIS/061/MCD64A1")
            .filterBounds(aoi)
            .filterDate(f"{y}-01-01", f"{y+1}-01-01")
            .select("BurnDate")
            .max()
            .gt(0)
            .unmask(0)
        )
        no_disturbance = loss_y.Or(burned).Not().rename(f"no_disturbance_{y}")
        bands[f"no_disturbance_{y}"] = no_disturbance

    image = ee.Image.cat(list(bands.values())).toByte()
    task = ee.batch.Export.image.toDrive(
        image=image, description="antar_disturbance_ancillary", folder=drive_folder,
        region=aoi, crs=crs, scale=scale_m, maxPixels=1e13,
    )
    task.start()
    return task


def export_terrain(bbox_wgs84, crs: str, scale_m: float, drive_folder: str = "antar_gee_exports"):
    """DEM, slope and aspect (Table 3/4). Returns the started task.

    The concept note specifies Copernicus GLO-30, but Earth Engine's copy of it
    (``COPERNICUS/DEM/GLO30``) has a real coverage gap over Armenia -- checked
    directly: only 3 of the ~14 tiles needed for the study bbox exist there,
    covering just the southwest corner. Used ``USGS/SRTMGL1_003`` (SRTM 30 m)
    instead, confirmed to have full bbox coverage (106,256 valid 1 km-sampled
    pixels, sensible 70-4978 m range) -- a substitution, not silently assumed;
    flagged here and in IMPLEMENTATION_LOG.md.

    TWI, TPI and the cold-air-pooling concavity index are not computed here:
    they are neighbourhood/flow-routing operations GEE's per-pixel
    ``ee.Terrain`` functions do not provide, and are left to local
    post-processing on the downloaded DEM.
    """
    aoi = ee.Geometry.Rectangle(list(bbox_wgs84))
    dem = ee.Image("USGS/SRTMGL1_003").select("elevation").clip(aoi)
    terrain = ee.Terrain.products(dem)
    image = terrain.select(["elevation", "slope", "aspect"]).toFloat()
    task = ee.batch.Export.image.toDrive(
        image=image, description="antar_terrain", folder=drive_folder, region=aoi, crs=crs, scale=scale_m, maxPixels=1e13,
    )
    task.start()
    return task


def export_soils(bbox_wgs84, crs: str, scale_m: float, drive_folder: str = "antar_gee_exports"):
    """SoilGrids 2.0 (Poggio et al. 2021) clay/sand/silt fraction and organic carbon,
    0-30 cm mean, for the Saxton & Rawls (2006) pedotransfer functions. Returns the
    started task.
    """
    aoi = ee.Geometry.Rectangle(list(bbox_wgs84))
    bands = {}
    for prop in ["clay", "sand", "silt", "soc"]:
        img = ee.Image(f"projects/soilgrids-isric/{prop}_mean")
        depths = ["0-5cm", "5-15cm", "15-30cm"]
        bands[prop] = img.select([f"{prop}_{d}_mean" for d in depths]).reduce(ee.Reducer.mean()).rename(f"{prop}_0_30cm_mean")
    image = ee.Image.cat(list(bands.values())).toFloat()
    task = ee.batch.Export.image.toDrive(
        image=image, description="antar_soils", folder=drive_folder, region=aoi, crs=crs, scale=scale_m, maxPixels=1e13,
    )
    task.start()
    return task


def export_era5land_forcing(bbox_wgs84, year_start: int, year_end: int, crs: str, scale_m: float,
                             drive_folder: str = "antar_gee_exports"):
    """Annual-mean wind speed, shortwave/longwave radiation, dewpoint temperature and
    surface pressure from ERA5-Land daily aggregates -- Table 4's ERA5-Land-sourced
    variables (never CHELSA, which does not carry them). Returns the started task.

    Reported as annual means (one band per year per variable) rather than the full
    daily series: TOPOHYDRO's own daily forcing comes from CHELSA/CORDEX (Sec. 5.1);
    these ERA5-Land variables feed terms (PET's aerodynamic term, the slope-radiation
    energy term, VPD's dew-point downscaling reference) that this annual-mean summary
    is enough to parameterise a first pass with -- the full daily series can be pulled
    the same way, at proportionally higher cost, if a later pass needs it.
    """
    aoi = ee.Geometry.Rectangle(list(bbox_wgs84))
    coll = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR").filterBounds(aoi)

    def _daily_wind_speed(img):
        return img.select("u_component_of_wind_10m").hypot(img.select("v_component_of_wind_10m")).rename("wind_speed")

    bands = {}
    for y in range(year_start, year_end + 1):
        yearly = coll.filterDate(f"{y}-01-01", f"{y + 1}-01-01")
        # Daily speed (hypot of that day's u,v) computed *before* averaging over the year --
        # averaging u and v separately first would let opposite-direction days cancel toward
        # zero instead of reflecting the actual mean wind speed.
        wind = yearly.map(_daily_wind_speed).mean()
        bands[f"wind_speed_{y}"] = wind.rename(f"wind_speed_{y}")
        bands[f"ssrd_{y}"] = yearly.select("surface_solar_radiation_downwards_sum").mean().rename(f"ssrd_{y}")
        bands[f"strd_{y}"] = yearly.select("surface_thermal_radiation_downwards_sum").mean().rename(f"strd_{y}")
        bands[f"dewpoint_{y}"] = yearly.select("dewpoint_temperature_2m").mean().rename(f"dewpoint_{y}")
        bands[f"surface_pressure_{y}"] = yearly.select("surface_pressure").mean().rename(f"surface_pressure_{y}")

    image = ee.Image.cat(list(bands.values())).toFloat()
    task = ee.batch.Export.image.toDrive(
        image=image, description="antar_era5land_forcing", folder=drive_folder,
        region=aoi, crs=crs, scale=scale_m, maxPixels=1e13,
    )
    task.start()
    return task


def export_vegetation_state(bbox_wgs84, year_start: int, year_end: int, crs: str, scale_m: float,
                             drive_folder: str = "antar_gee_exports"):
    """Growing-season-mean MODIS LAI/FPAR, annual land-surface phenology (start/end of
    season), and ESA WorldCover tree-cover class (a single, most-recent-epoch layer;
    WorldCover is not produced annually). Returns the started task.
    """
    aoi = ee.Geometry.Rectangle(list(bbox_wgs84))
    lai_coll = ee.ImageCollection("MODIS/061/MOD15A2H").filterBounds(aoi)
    pheno_coll = ee.ImageCollection("MODIS/061/MCD12Q2").filterBounds(aoi)

    bands = {}
    for y in range(year_start, year_end + 1):
        lai = lai_coll.filterDate(f"{y}-07-01", f"{y}-09-01").select("Lai_500m").mean()
        bands[f"lai_growing_season_{y}"] = lai.rename(f"lai_growing_season_{y}")
        pheno = pheno_coll.filterDate(f"{y}-01-01", f"{y + 1}-01-01").select(["Greenup_1", "Dormancy_1"]).first()
        bands[f"greenup_doy_{y}"] = pheno.select("Greenup_1").rename(f"greenup_doy_{y}")
        bands[f"dormancy_doy_{y}"] = pheno.select("Dormancy_1").rename(f"dormancy_doy_{y}")

    worldcover = ee.ImageCollection("ESA/WorldCover/v200").first().select("Map").rename("worldcover_class")
    image = ee.Image.cat(list(bands.values()) + [worldcover]).toFloat()
    task = ee.batch.Export.image.toDrive(
        image=image, description="antar_vegetation_state", folder=drive_folder,
        region=aoi, crs=crs, scale=scale_m, maxPixels=1e13,
    )
    task.start()
    return task


def export_snow(bbox_wgs84, year_start: int, year_end: int, crs: str, scale_m: float,
                 drive_folder: str = "antar_gee_exports"):
    """Annual snow-cover fraction (share of days with NDSI-flagged snow) and duration
    (day count) from MODIS MOD10A1. Returns the started task.
    """
    aoi = ee.Geometry.Rectangle(list(bbox_wgs84))
    coll = ee.ImageCollection("MODIS/061/MOD10A1").filterBounds(aoi)

    bands = {}
    for y in range(year_start, year_end + 1):
        yearly = coll.filterDate(f"{y}-01-01", f"{y + 1}-01-01").select("NDSI_Snow_Cover")
        is_snow = yearly.map(lambda img: img.gte(40).rename("snow"))  # NDSI >= 40 is the standard MOD10A1 snow flag
        bands[f"snow_days_{y}"] = is_snow.sum().rename(f"snow_days_{y}")
        bands[f"snow_frac_{y}"] = is_snow.mean().rename(f"snow_frac_{y}")

    image = ee.Image.cat(list(bands.values())).toFloat()
    task = ee.batch.Export.image.toDrive(
        image=image, description="antar_snow", folder=drive_folder, region=aoi, crs=crs, scale=scale_m, maxPixels=1e13,
    )
    task.start()
    return task


def export_land_tenure(bbox_wgs84, crs: str, scale_m: float, drive_folder: str = "antar_gee_exports"):
    """WDPA protected-area boundaries (Table 4), rasterised to a boolean "protected" band
    on the master grid -- feeds ``configs/study_area.yaml``'s eligibility-mask exclusions
    directly. Returns the started task.
    """
    aoi = ee.Geometry.Rectangle(list(bbox_wgs84))
    wdpa = ee.FeatureCollection("WCMC/WDPA/current/polygons").filterBounds(aoi)
    protected = ee.Image(0).paint(wdpa, 1).rename("protected_area").toByte().clip(aoi)
    task = ee.batch.Export.image.toDrive(
        image=protected, description="antar_land_tenure", folder=drive_folder,
        region=aoi, crs=crs, scale=scale_m, maxPixels=1e13,
    )
    task.start()
    return task


def export_structure(
    bbox_wgs84,
    crs: str,
    scale_m: float,
    drive_folder: str = "antar_gee_exports",
    gedi_year_start: int = 2019,
    gedi_year_end: int = 2023,
    canopy_height_asset: str = "users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1",
):
    """Canopy-structure training data for MERISTEM's attainable-height model, exported
    to Google Drive.

    GEDI L2A RH98 (quality-filtered: quality_flag=1, degrade_flag=0, sensitivity>0.9)
    is the direct lidar measurement of canopy-top height the quantile regression
    should be trained against; its footprints are sparse (along orbital tracks), so
    the Lang et al. (2023) global canopy height model (10 m, wall-to-wall, itself
    partly trained on GEDI) is exported alongside it to fill the gaps. WSL S2-VHM,
    named in this module's original stub docstring, has no confirmed public Earth
    Engine asset (checked directly, not found under any plausible id) and is not
    included -- a real gap, not silently substituted.

    Returns the started ``ee.batch.Task``.
    """
    aoi = ee.Geometry.Rectangle(list(bbox_wgs84))

    gedi = (
        ee.ImageCollection("LARSE/GEDI/GEDI02_A_002_MONTHLY")
        .filterBounds(aoi)
        .filterDate(f"{gedi_year_start}-01-01", f"{gedi_year_end + 1}-01-01")
    )

    def quality_masked_rh98(img):
        good = (
            img.select("quality_flag").eq(1)
            .And(img.select("degrade_flag").eq(0))
            .And(img.select("sensitivity").gt(0.9))
        )
        return img.select("rh98").updateMask(good)

    gedi_rh98 = gedi.map(quality_masked_rh98)
    gedi_median = gedi_rh98.median().rename("gedi_rh98_median")
    gedi_count = gedi_rh98.count().rename("gedi_valid_count")

    chm = ee.Image(canopy_height_asset).rename("eth_chm_2020").clip(aoi)

    image = ee.Image.cat([gedi_median, gedi_count, chm]).toFloat()
    task = ee.batch.Export.image.toDrive(
        image=image, description="antar_structure", folder=drive_folder,
        region=aoi, crs=crs, scale=scale_m, maxPixels=1e13,
    )
    task.start()
    return task
