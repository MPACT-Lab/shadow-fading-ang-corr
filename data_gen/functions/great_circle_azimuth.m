function [azim_wrt_BS] = great_circle_azimuth(lat_BS,lon_BS,lat_UAV,lon_UAV)
% https://dtcenter.org/sites/default/files/community-code/met/docs/write-ups/gc_simple.pdf
%GREAT_CIRCLE_AZIMUTH Summary of this function goes here
%   we will calculate azimuth wrt BS
% first we will calculate the bearing wrt the BS
% but bearing is calculated eastward from north, north is 0 degree, east is
% 90 degree, west is 270 degree
% for calculating del L, Lat2 - Lat1, where east Longitude is considered
% positive
S = cos( pi/180*lat_UAV ) .* sin( pi/180*(lon_UAV - lon_BS) ) ;
C = cos( pi/180*lat_BS ) .* sin( pi/180*lat_UAV ) - sin( pi/180*lat_BS ) .* cos( pi/180*lat_UAV ) .* cos( pi/180*(lon_UAV - lon_BS) );
bearing = atan2(S , C) ;
azim_wrt_BS = - bearing + pi/2; % - bearing makes counter clockwise and +pi/2 shifts zero to +x
azim_wrt_BS = azim_wrt_BS - 2*pi.*(azim_wrt_BS > pi); % within -pi and pi
end

