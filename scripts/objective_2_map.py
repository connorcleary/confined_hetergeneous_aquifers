import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt

figsize = [4, 5]
text_size = 10.
fig = plt.figure(figsize=figsize)
margin_x = [0.01, 0.01]
margin_y = [0.01, 0.01]
ax = plt.axes([margin_x[0], margin_y[0], 1-sum(margin_x), 1-sum(margin_y)], projection=ccrs.PlateCarree())
x_extent = [166, 179]              # longitude plot limits
y_center = -41                   # center of lattitude axis (limits computed to be consistent with lon extent) 
y_delta = figsize[1]/figsize[0]*(x_extent[1]-x_extent[0])
y_extent = [y_center-y_delta/2, y_center+y_delta/2]
ax.set_extent(x_extent+y_extent, ccrs.PlateCarree())

# add lines
ax.coastlines(resolution='50m')

# add coloured land and ocean
land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='k',facecolor=cfeature.COLORS['land'])
ax.add_feature(land_50m, linewidth=.01)
ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m', edgecolor='face',facecolor=cfeature.COLORS['water'])
ax.add_feature(ocean_50m)

hb = {'lon': 176.7, 'lat':-39.4, 'txt':'Heretaunga Plains', 'side':"left"}
hv = {'lon': 174.9, 'lat':-41.2, 'txt':'Hutt Valley', 'side':"right"}
mv = {'lon': 173.0, 'lat':-41.2, 'txt':"Moutere Valley", 'side':"left"}
wps = {'lon': 173.0, 'lat':-41.4, 'txt':"Waimea Plains", 'side':"left"}
wp = {'lon': 174.0, 'lat':-41.5, 'txt':'Wairau Plain', 'side':"left"}
wap = {'lon': 172.5, 'lat':-43.3, 'txt':'Waimakariri-Ashley Plains', 'side':"left"}
cp = {'lon': 172.6, 'lat':-43.6, 'txt':'Central Plains', 'side':"left"}
arp = {'lon': 171.6, 'lat':-44, 'txt':'Ashburton-Rangitata Plains', 'side':"left"}



for EGS in [hb, hv, mv, wps, wp, wap, cp, arp]:
    # we'll use 'scatter' to add the points
    ax.scatter(EGS['lon'], EGS['lat'], marker = 'o', c='r', zorder = 5, s = 9, lw=.1, edgecolors='k')
    # and 'text' to annotate them
    # if EGS["side"] == 'left':
    #     ax.text(EGS['lon']-0.5, EGS['lat'], EGS['txt'], fontsize = text_size, va='center', ha='right')
    # else:
    #     ax.text(EGS['lon']+0.5, EGS['lat'], EGS['txt'], fontsize = text_size, va='center', ha='left')

plt.show()