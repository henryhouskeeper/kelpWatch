#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 11:12:30 2022

@author: floatingforest_mbp

HFH notes: Use conda environment kw:
(base) $ conda create -n kw -c conda-forge python=3.9
(base) $ conda activate kw
(kw) $ conda config --set channel_priority strict
(kw) $ conda install -c conda-forge netcdf4
(kw) $ conda install -c conda-forge gdal
(kw) $ conda install -c conda-forge shapely
(kw) $ pip install rio-cogeo
(kw) $ conda install spyder-kernels=1.10.0
(kw) $ pip install rasterio

## Did not run: (kw) $ conda install -c conda-forge rasterio
## bc it caused multiple package incompatibilities

"""

from warnings import filterwarnings
filterwarnings("ignore")

from IPython import get_ipython;   
get_ipython().magic('reset -sf')

import os
import netCDF4 as nc
import numpy as np
#from pyproj import CRS
from shapely.geometry import MultiPoint

#from osgeo import gdal

import rasterio as rio
from rasterio.io import MemoryFile
from rasterio.transform import from_bounds
#from rasterio.features import rasterize
#from rasterio.features import shapes

from rio_cogeo.cogeo import cog_translate
from rio_cogeo.profiles import cog_profiles


##############################################################################
#              Define parameters for script and GeoTiff output               #
##############################################################################
VIZ_TOGGLE = False
nanval = -32768
res= 30
nBands = 8
bandNames = ['Area','Area_SE','Biomass','Biomass_SE',
               'Passes','Passes4_5','Passes7','Passes8']
bandUnits = ['m^2','m^2','wet Kg/900m^2 pixel','wet Kg/900m^2 pixel',
                  'number of images used','number of images used',
                  'number of images used','number of images used']

##############################################################################
# Set working dir to location of this script version and define output dir:  #
##############################################################################
path = os.path.abspath(__file__).rsplit('/',1)[0]
print(f'\nWorking in:\n\n\t{path}')

in_dir = f'{path}/input'

#out_dir = f'{path}/output'
out_dir = '/volumes/Planktos/kelpWatch'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    
##############################################################################
#                            Import utm data:                                #
##############################################################################
fname_utm_coords = f'{in_dir}/utm_coords.csv'

utm_coords = np.genfromtxt(fname_utm_coords, delimiter=',',skip_header=1)

#utm_coords = utm_coords[0:500,:]

#x = utm_coords[:,0]
#y = utm_coords[:,1]
utm = utm_coords[:,2].astype(int)

print('\nUTM information read from:')
print('\n\t',fname_utm_coords)

##############################################################################
#                          Import netcdf data:                               #
##############################################################################

ncname = f'{in_dir}/kelpCanopyFromLandsat_2021_v2.nc'

ds = nc.Dataset(ncname)

print('\nNetCDF information read (netcdf4) from:')
print('\n\t',ncname)

print('\nNetCDF variables:')
print('\n\t',ds.variables.keys())


year = ds['year']
quarter = ds['quarter']

##########################################################
#           Import data with dim: (geo)                  #
##########################################################

# index by utm
us = np.unique(utm)
for ku in range(len(us)):
  print('\n-------------------------------')
  print(f'|\tProcessing UTM zone: {us[ku]}\t|')
  print('-------------------------------')
  
  iu = np.where(utm == us[ku])[0]
  
  x = utm_coords[iu,0].astype(int)
  y = utm_coords[iu,1].astype(int)
  lat = ds['lat'][iu]
  lon = ds['lon'][iu]
  
  # Define CRS:
  if (us[ku] == 10) & (np.nanmin(lat) >= 0):
     crsOut = rio.crs.CRS.from_epsg(32610)
  elif (us[ku] == 11) & (np.nanmin(lat) >= 0):
    crsOut = rio.crs.CRS.from_epsg(32611)   
  else:
    raise Exception('CRS for this location not yet coded.')
  '''crsOut = rio.crs.CRS.from_dict(
      proj='utm', zone=us[ku], datum='WGS84')'''
  
  if (np.nanmin(lat) >= 0):
    hemi = 'N'
  else:
    hemi = 'S'
    
  xmin = np.amin(x)
  xmax = np.amax(x)
  ymin = np.amin(y)
  ymax = np.amax(y)
  
  width = len(np.arange(xmin,xmax,res))
  height = len(np.arange(ymin,ymax,res))
  
  print(f'\n\tUTM: {us[ku]}')
  print(f'\tx bounds: {xmin} -- {xmax}')
  print(f'\ty bounds: {ymin} -- {ymax}')
  print(f'\traster dims: {width} x {height}')
  print(f'\t(N points in vector: {len(x)})')
  print(f'\t(N points in raster: {height*width})')
  
  xyPoints = MultiPoint(list(tuple(zip(x,y))))
  
  transform = (rio.transform.Affine.translation(xmin - res / 2, ymin - res / 2)
    * rio.transform.Affine.scale(res, res))
  
  ##########################################################
  #         Import data with dim: (tim,geo)                #
  ##########################################################
  
  # index by year: 28
  years = np.unique(year)
  for ky in range(len(years)):
    iy = np.where(year == years[ky])[0]
    print(f'\n\tRasterizing year [{years[ky]}]: quarter ',end='')
    
    # index by quarter: 0
    quarters = np.unique(quarter)
    for kq in range(len(quarters)):
      # merge year and quarter indices:
      it = np.intersect1d(iy,np.where(quarter == quarters[kq])[0])
      
      print(f'[{quarters[kq]}]',end='')
      
      fout = (f'{out_dir}/kelpArea_'
        f'{us[ku]:02d}{hemi:s}'
        f'_{years[ky]:04d}'
        f'_{quarters[kq]:02d}.tif')
      
      # extract data for 1 utm, 1 quarter, 1 year:
      area = ds['area'][it,iu].astype(np.int16)[0] #int16
      area_se = ds['area_se'][it,iu].astype(np.int16)[0] #int16
      biomass = ds['biomass'][it,iu].astype(np.int16)[0] #int16
      biomass_se = ds['biomass_se'][it,iu].astype(np.int16)[0] #int16
      passes = ds['passes'][it,iu].astype(np.int16)[0] #int8
      passes4_5 = ds['passes4_5'][it,iu].astype(np.int16)[0] #int8
      passes7 = ds['passes7'][it,iu].astype(np.int16)[0] #int8
      passes8 = ds['passes8'][it,iu].astype(np.int16)[0] #int8
      
      ##########################################################
      #                 Rasterize point data                   #
      ##########################################################
      
      src_transform = from_bounds(
        west=xmin,east=xmax,
        south=ymin,north=ymax,
        width=width, height=height)
      
      src_profile = dict(
        driver="GTiff",
        dtype="int16",
        count=nBands,
        height=height,
        width=width,
        crs=crsOut,
        transform=src_transform)
      
      areaRaster = rio.features.rasterize(
        ((geom,value) for geom, value in tuple(zip(xyPoints, area))),
        out_shape=tuple((height,width)),
        fill=nanval,
        all_touched = True,
        transform = src_transform).astype(np.int16)
      
      areaSeRaster = rio.features.rasterize(
        ((geom,value) for geom, value in tuple(zip(xyPoints, area_se))),
        out_shape=tuple((height,width)),
        fill=nanval,
        all_touched = True,
        transform = src_transform).astype(np.int16)
      
      biomassRaster = rio.features.rasterize(
        ((geom,value) for geom, value in tuple(zip(xyPoints, biomass))),
        out_shape=tuple((height,width)),
        fill=nanval,
        all_touched = True,
        transform = src_transform).astype(np.int16)
      
      biomassSeRaster = rio.features.rasterize(
        ((geom,value) for geom, value in tuple(zip(xyPoints, biomass_se))),
        out_shape=tuple((height,width)),
        fill=nanval,
        all_touched = True,
        transform = src_transform).astype(np.int16)
      
      passesRaster = rio.features.rasterize(
        ((geom,value) for geom, value in tuple(zip(xyPoints, passes))),
        out_shape=tuple((height,width)),
        fill=nanval,
        all_touched = True,
        transform = src_transform).astype(np.int16)
      
      passes4_5Raster = rio.features.rasterize(
        ((geom,value) for geom, value in tuple(zip(xyPoints, passes4_5))),
        out_shape=tuple((height,width)),
        fill=nanval,
        all_touched = True,
        transform = src_transform).astype(np.int16)
      
      passes7Raster = rio.features.rasterize(
        ((geom,value) for geom, value in tuple(zip(xyPoints, passes7))),
        out_shape=tuple((height,width)),
        fill=nanval,
        all_touched = True,
        transform = src_transform).astype(np.int16)
      
      passes8Raster = rio.features.rasterize(
        ((geom,value) for geom, value in tuple(zip(xyPoints, passes8))),
        out_shape=tuple((height,width)),
        fill=nanval,
        all_touched = True,
        transform = src_transform).astype(np.int16)
      
      dataStack = np.stack((areaRaster,areaSeRaster,biomassRaster,biomassSeRaster,
        passesRaster,passes4_5Raster,passes7Raster,passes8Raster))
      
      del areaRaster, areaSeRaster, biomassRaster, biomassSeRaster
      del passes4_5Raster, passes7Raster, passes8Raster
      
      ##########################################################
      #              Output Raster Files as COGs               #
      ##########################################################
      
      with MemoryFile() as memfile:
        with memfile.open(**src_profile) as mem:
          mem.write(dataStack)
          for kb in range(0,nBands):
            mem.set_band_description(kb+1,bandNames[kb])
            mem.set_band_unit(kb+1,bandUnits[kb])
          
          #dst_profile = cog_profiles.get("deflate")
          dst_profile = cog_profiles.get("lzw")
          cog_translate(
            mem,
            fout,
            dst_profile,
            forward_band_tags=True,
            quiet=True,
            nodata=nanval)
  
      del mem, dataStack

if VIZ_TOGGLE:
  os.system("rio cogeo info output/kelpArea_zone_11N_year_2012_quarter_01.tif")
  # https://github.com/developmentseed/rio-viz
  os.system("rio viz output/kelpArea_zone_10N_year_2012_quarter_01.tif")

def _translate(src_path, dst_path, profile="webp", profile_options={}, **options):
    """Convert image to COG."""
    # Format creation option (see gdalwarp `-co` option)
    output_profile = cog_profiles.get(profile)
    output_profile.update(dict(BIGTIFF="IF_SAFER"))
    output_profile.update(profile_options)

    # Dataset Open option (see gdalwarp `-oo` option)
    config = dict(
        GDAL_NUM_THREADS="ALL_CPUS",
        GDAL_TIFF_INTERNAL_MASK=True,
        GDAL_TIFF_OVR_BLOCKSIZE="128",
    )

    cog_translate(
        src_path,
        dst_path,
        output_profile,
        config=config,
        in_memory=False,
        quiet=True,
        **options,
    )
    return True