# Directory locations
HDFDIR = '/home/erwin/MOD35L2/hdfs/'
SAVEDIR = '/home/erwin/MOD35L2/plots/'

# Load all necessary modules / packages
from pylab import *
from pyhdf.SD import SD
from scipy.interpolate import RectBivariateSpline

# Function that extracts desired matrices from MOD35L2 HDF SDS
def readMOD35L2(fname, geoloc_only=False):
    hdf_file = SD(HDFDIR + fname)
    if not geoloc_only:
        cloud_mask = hdf_file.select('Cloud_Mask').get()
    lon = hdf_file.select('Longitude').get()
    lat = hdf_file.select('Latitude').get()
    hdf_file.end()
    
    if not geoloc_only:
        cld_msk = uint8(cloud_mask[0])
        cloud = cld_msk & 6 # 0, 2, 4, 6
        land = cld_msk & 192 # 0, 64, 128, 192
    
        cloud[cloud==0] = 1 # 0 -> Cloud
        cloud[cloud!=1] = 0 # 2, 4, 6 -> No cloud

        coast = land
        coast[coast==64] = 1 # 64 -> Coast
        coast[coast!=1] = 0 # 0, 128, 192 -> Not coast

        land[land!=0] = 1 # 64, 128, 192 -> Land, 0 -> Water
        
        return lon, lat, cloud, land, coast

    return lon, lat

# Function that auto-inverts geographically upside-down MOD35L2 HDF
def autoinvertMOD35L2(lon, lat, cloud, land, coast):
    if lon[0][0] > lon[0][-1]: # left-right inverted
        lon = fliplr(lon)
        lat = fliplr(lat)
        cloud = fliplr(cloud)
        land = fliplr(land)
        coast = fliplr(coast)
    if lat[0][0] < lat[-1][0]: # up-down inverted
        lon = flipud(lon)
        lat = flipud(lat)
        cloud = flipud(cloud)
        land = flipud(land)
        coast = flipud(coast)
    
    return lon, lat, cloud, land, coast

# Function that previews cloud & geographic image of a MOD35L2 HDF
def previewMOD35L2(fname):
    lon, lat, cld, lnd, cst = readMOD35L2(fname)
    lon, lat, cld, lnd, cst = autoinvertMOD35L2(lon, lat, cld, lnd, cst)
    cld[cld==1] = 100 # Cloud
    cld[cld==0] = 50 # No cloud
    cld[cst==1] = 0 # Overlay coast
    
    tl_lon = round(lon[0][0], 2) # top left
    tl_lat = round(lat[0][0], 2)
    tr_lon = round(lon[0][-1], 2) # top right
    tr_lat = round(lat[0][-1], 2)
    bl_lon = round(lon[-1][0], 2) # bot left
    bl_lat = round(lat[-1][0], 2)
    br_lon = round(lon[-1][-1], 2) # bot right
    br_lat = round(lat[-1][-1], 2)
    
    tl = '(' + str(tl_lon) + ',' + str(tl_lat) + ')'
    tr = '(' + str(tr_lon) + ',' + str(tr_lat) + ')'
    bl = '(' + str(bl_lon) + ',' + str(bl_lat) + ')'
    br = '(' + str(br_lon) + ',' + str(br_lat) + ')'

    fig = figure()
    img = imshow(cld, cmap=cm.hot)
    xticks([])
    yticks([])
    x, y = shape(cld)
    m, n = y, x
    text(0, -5, tl, ha='left', va='bottom')
    text(m, -5, tr, ha='right', va='bottom')
    text(0, n+5, bl, ha='left', va='top')
    text(m, n+5, br, ha='right', va='top')
    savefig(SAVEDIR + fname + '_preview.png')
    close(fig)

# Function that finds corners of MOD35L2 group's common area
def MOD35L2group_corners(hdf_fnames):
    lons = []
    lats = []
    for hdf_fname in hdf_fnames:
        print hdf_fname
        lon, lat = readMOD35L2(hdf_fname, geoloc_only=True)
        dummy = zeros((2,2))
        lon, lat, dummy, dummy, dummy = autoinvertMOD35L2(lon, lat, \
                                                          dummy, dummy, dummy)
        if shape(lon) != (406, 270):
            lon = reshape_matrix(lon, 406, 270)
            lat = reshape_matrix(lat, 406, 270)
        lons.append(lon)
        lats.append(lat)
    lons = array(lons)
    lats = array(lats)
    gr_tl_lon = lons[:,0,0].max()
    gr_tl_lat = lats[:,0,0].min()
    gr_br_lon = lons[:,-1,-1].min()
    gr_br_lat = lats[:,-1,-1].max()

    print gr_tl_lon, gr_tl_lat, gr_br_lon, gr_br_lat

    if gr_tl_lon > gr_br_lon or gr_tl_lat < gr_br_lat:
        return NaN, NaN, NaN, NaN # No common area

    return gr_tl_lon, gr_tl_lat, gr_br_lon, gr_br_lat

# Function that intersects area bounded by tl & br to Phil Area
def intersect_phil(tl_lon, tl_lat, br_lon, br_lat):
    intrsct_tl_lon = max(tl_lon, 115.)
    intrsct_tl_lat = min(tl_lat, 22.5)
    intrsct_br_lon = min(br_lon, 130.)
    intrsct_br_lat = max(br_lat, 2.5)

    if intrsct_tl_lon > intrsct_br_lon or intrsct_tl_lat < intrsct_br_lat:
        return NaN, NaN, NaN, NaN # Null intersection

    return intrsct_tl_lon, intrsct_tl_lat, intrsct_br_lon, intrsct_br_lat

# Function that
def processMOD35L2group(hdf_fnames, l, t, r, b):
    lons = []
    lats = []
    clds = []
    lnds = []
    csts = []

    print 'Extracting matrices...'
    i = 1
    for hdf_fname in hdf_fnames:
        print i, hdf_fname
        i += 1
        lon, lat, cld, lnd, cst = readMOD35L2(hdf_fname)
        lon, lat, cld, lnd, cst = autoinvertMOD35L2(lon, lat, cld, lnd,cst)
        lons.append(lon)
        lats.append(lat)
        clds.append(cld)
        lnds.append(lnd)
        csts.append(cst)
    lons = array(lons)
    lats = array(lats)
    clds = array(clds)
    lnds = array(lnds)
    csts = array(csts)
    
    land_cfs = []
    water_cfs = []
    total_cfs = []
    x = []
    y = []
    N = len(clds)
    print 'Processing matrices...'
    for i in arange(N):
        print i + 1, hdf_fname
        m, n = shape(clds[i])
        lons[i] = reshape_matrix(lons[i], m, n)
        lats[i] = reshape_matrix(lats[i], m, n)
        I, J, K, L = index_of_closest_corners(lons[i], lats[i], l, t, r, b)
        
        lons[i] = trim_matrix(lons[i], I, J, K, L)
        lats[i] = trim_matrix(lats[i], I, J, K, L)
        clds[i] = trim_matrix(clds[i], I, J, K, L)
        lnds[i] = trim_matrix(lnds[i], I, J, K, L)
        csts[i] = trim_matrix(csts[i], I, J, K, L)
        m1, n1 = shape(clds[i])
        x.append(m1)
        y.append(n1)

        land_loc = find(lnds[i] == 1)
        water_loc = find(lnds[i] == 0)
        land_cfs.append(sum(clds[i].flat[land_loc]) / size(land_loc))
        water_cfs.append(sum(clds[i].flat[water_loc]) / size(water_loc))
        total_cfs.append(sum(clds[i]) / size(clds[i]))
    land_cfs = array(land_cfs) # Cloud fraction over land
    water_cfs = array(water_cfs) # Cloud fraction over water
    total_cfs = array(total_cfs) # Total cloud fraction on cloud map/image
    x = array(x)
    y = array(y)
    m2 = max(x)
    n2 = max(y)
   
    print 'Reshaping matrices...'
    for i in arange(N):
        print i + 1, hdf_fnames[i]
        clds[i] = reshape_matrix(clds[i], m2, n2)
        csts[i] = reshape_matrix(csts[i], m2, n2)

    print 'Plotting...'
    tl_lon = round(mean(lons[:][0][0]), 2)
    tl_lat = round(mean(lats[:][0][0]), 2)
    tr_lon = round(mean(lons[:][0][-1]), 2)
    tr_lat = round(mean(lats[:][0][-1]), 2)
    bl_lon = round(mean(lons[:][-1][0]), 2)
    bl_lat = round(mean(lats[:][-1][0]), 2)
    br_lon = round(mean(lons[:][-1][-1]), 2)
    br_lat = round(mean(lats[:][-1][-1]), 2)

    tl = '(' + str(tl_lon) + ',' + str(tl_lat) + ')'
    tr = '(' + str(tr_lon) + ',' + str(tr_lat) + ')'
    bl = '(' + str(bl_lon) + ',' + str(bl_lat) + ')'
    br = '(' + str(br_lon) + ',' + str(br_lat) + ')'

    fig = figure()
    cld_map = clds.mean(0) * 100
    cst_map = csts.mean(0)
    cld_map[cst_map>0] = 0 # Overlay coastlines to cloud map
    img = imshow(cld_map, cmap=cm.jet)
    xticks([])
    yticks([])
    colorbar(mappable=img)
    x, y = shape(cld_map)
    m, n = y, x
    text(0, -5, tl, ha='left', va='bottom')
    text(m, -5, tr, ha='right', va='bottom')
    text(0, n+5, bl, ha='left', va='top')
    text(m, n+5, br, ha='right', va='top')
    savefig(SAVEDIR + hdf_fnames[0] + '_cloudmap.png')
    close(fig)

    print 'Writing Cloud Fraction values...'
    land_cf = round(land_cfs.mean() * 100, 2)
    water_cf = round(water_cfs.mean() * 100, 2)
    total_cf = round(total_cfs.mean() * 100, 2)
    output_file = file(SAVEDIR + 'cf_values.txt', 'a')
    output_file.write(hdf_fnames[0] + '-' \
                      + str(land_cf) + '-' \
                      + str(water_cf) + '-' \
                      + str(total_cf) + '\n')
    output_file.close()

    print 'Done.'

# Function that reshapes a matrix using linear rect bivariate spline interp
def reshape_matrix(input_matrix, x, y):
    m, n = shape(input_matrix)
    i = arange(m)
    j = arange(n)
    ii = linspace(i.min(), i.max(), x)
    jj = linspace(j.min(), j.max(), y)
    matrix_spline = RectBivariateSpline(i, j, input_matrix, kx=1, ky=1)
    
    return matrix_spline(ii, jj)

# Function that finds the indices of points nearest to specified corners
def index_of_closest_corners(lon, lat, l, t, r, b):
    # define distance^2 as abs val of the difference bet two complex numbers
    # then, find the minimum distance and get the index of that closest corner
    I, J = unravel_index(abs((lon + lat*1j) - (l + t*1j)).argmin(), shape(lon))
    K, L = unravel_index(abs((lon + lat*1j) - (r + b*1j)).argmin(), shape(lon))
    return I, J, K, L

# Function that trims a matrix
def trim_matrix(input_matrix, I, J, K, L):
    m, n = shape(input_matrix)
    input_matrix = delete(input_matrix, arange(K + 1, m), 0)
    input_matrix = delete(input_matrix, arange(0, I), 0)
    input_matrix = delete(input_matrix, arange(L + 1, n), 1)
    input_matrix = delete(input_matrix, arange(0, J), 1)

    return input_matrix

