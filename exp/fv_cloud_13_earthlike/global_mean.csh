#!/bin/csh -f

#set filename = earth_e.nc #6000w_mean.nc
set filename = day9490h00/day9490h00.atmos_echeck.nc
ncap2 -h -O -s "weights=cos(grid_yt*3.1415/180)" $filename in.nc
ncwa -h -O -w weights -a grid_yt,grid_xt in.nc e_global_mean.nc
#ncdump -v OLR global_mean.nc
#foreach var (OLR, raddown, radup, flux_t)
foreach var (OLR, OSR)
   ncdump -v $var e_global_mean.nc | grep "$var ="
end

#set filename = day3200h00/day3200h00.atmos_static.nc

#ncap -h -O -s "weights=cos(grid_yt*3.1415/180)" $filename in.nc
#ncwa -h -O -w weights -a grid_yt,grid_xt in.nc global_mean.nc
#ncdump -v solar,area global_mean.nc
