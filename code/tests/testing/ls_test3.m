function ls_test3()
%function to test stuff... can be deleted afterwards


%get bigger study area (srtm1)
[SRTM,srtm_info] = climada_srtm1_get([8.111086 8.341411 46.815369 46.946240],'','','',1);

Z = SRTM.h;
lon = SRTM.x;
lat = SRTM.y;

DEM = GRIDobj(lon,lat,Z);

GRIDobj2geotiff(DEM);

disp('hier')

end