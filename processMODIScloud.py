# Directory locations
HDF_DIR = '/media/Expansion Drive/erwin/ladsweb/hdfs/'
LIST_DIR = '/home/erwin/gits/crap/groups/pruned/'
SAVE_DIR = '/home/erwin/gits/crap/save/pruned/'
#HDF_DIR = '/home/erwin/gits/crap/rats/'
#LIST_DIR = '/home/erwin/gits/crap/rats/lists/'
#SAVE_DIR = '/home/erwin/gits/crap/rats/pics/'

# Load all necessary modules/packages
from pylab import *
from pyhdf.SD import SD

# Main program
def process_HDF_group(HDF_list, intersect_phil=True):
    print HDF_list + ':\n\tLOOKING FOR COMMON AREA'
    top, bottom, left, right = check_common_area(HDF_list)
    
    if isnan(top and bottom) or isnan(left and right):
        print 'ERROR: No common area.\n'
        return
    
    print '\tCommon area (top, bottom, left, right):' + '\n\t\t' \
            + str(round(top,2)) + ', ' \
            + str(round(bottom,2)) + ', ' \
            +str(round(left,2)) + ', ' \
            +str(round(right,2))
    
    if intersect_phil:
        print '\tINTERSECTING TO PHILIPPINE AREA'
        top, bottom, left, right = intersect_phil_area(top, bottom, left, right)
        if isnan(top and bottom) or isnan(left and right):
            print 'ERROR: No common area.\n'
            return
        else:
            print '\tIntersection (top, bottom, left, right): ' + '\n\t\t' \
                    + str(round(top,2)) + ', ' \
                    + str(round(bottom,2)) + ', ' \
                    +str(round(left,2)) + ', ' \
                    +str(round(right,2))
    
    print '\n\tCOMPUTING CLOUD FRACTION'
    get_cloud_fraction(HDF_list, top, bottom, left, right)
    print 'Done.\n'

# Use to extract mulitple Scientific Data Sets from mulitple HDF files
def read_HDF_SDS(HDF_list, *data_sets):
    input_file = file(LIST_DIR + HDF_list, 'r')
    HDF_file_names = input_file.readlines()
    input_file.close()
    
    HDF_SDS = []
    sds = ''
    for data_set in data_sets:
        sds += data_set + ' '
        HDF_SDS.append([])
    
    print '\tReading HDF SDS:\n\t\t' + sds
    for HDF_file_name in HDF_file_names:
        HDF_file = SD(HDF_DIR + HDF_file_name.rstrip())
        SDS_index = 0
        for data_set in data_sets:
            HDF_data_set = HDF_file.select(data_set)
            HDF_SDS[SDS_index].append(HDF_data_set.get())
            SDS_index += 1
        HDF_data_set.endaccess()
        HDF_file.end()
    
    return HDF_SDS

# Use to generate an image preview of the HDF files
def preview_HDF_files(HDF_list):
    lat, lon, cloud_mask = read_HDF_SDS(HDF_list, \
                                                                'Latitude',\
                                                                'Longitude',\
                                                                'Cloud_Mask')
    input_file = file(LIST_DIR + HDF_list, 'r')
    file_names = input_file.readlines()
    input_file.close()
    
    print 'Generating image previews...'
    lat, lon, cloud_mask = correct_inversion(lat, lon, cloud_mask)
    cloudiness = bitmask_cloud_mask(cloud_mask, 6)
    boundary_mask = bitmask_cloud_mask(cloud_mask, 192)
    del cloud_mask
    N = len(lat)
    for r in range(N):
        print file_names[r][0:22] + '.png'
        tl_lat = round(lat[r][0][0],2)
        tl_lon = round(lon[r][0][0],2)
        m, n = shape(lat[r])
        br_lat = round(lat[r][m-1][n-1],2)
        bl_lat = round(lat[r][m-1][0],2)
        tr_lat = round(lat[r][0][n-1],2)
        m, n = shape(lon[r])
        br_lon = round(lon[r][m-1][n-1],2)
        bl_lon = round(lon[r][m-1][0],2)
        tr_lon = round(lon[r][0][n-1],2)
        boundary_mask[r][boundary_mask[r]!=64] = 0
        boundary_mask[r][boundary_mask[r]==64] = 1  # boundary
        cloudiness[r][cloudiness[r]!=0] = 50  # land or desert
        cloudiness[r][cloudiness[r]==0] = 100  # cloud
        cloudiness[r][boundary_mask[r]==1] = 0  # overlay boundary
        save_HDF_image(cloudiness[r], file_names[r], \
                                    tl_lat, tl_lon, tr_lat, tr_lon, \
                                    bl_lat, bl_lon, br_lat, br_lon, \
                                    cm.hot)
    print 'Done.\n'

def save_HDF_image(image_map, file_name, \
                                        tl_lat, tl_lon, tr_lat, tr_lon, \
                                        bl_lat, bl_lon, br_lat, br_lon, \
                                        color_map):
    fig = figure()
    img = imshow(image_map, cmap=color_map)
    xticks([])
    yticks([])
    x, y = shape(image_map)
    m, n = y, x
    text(0, -5, '(' + str(tl_lon) + ',' + str(tl_lat) + ')', \
            ha='left', va='bottom')
    text(m, -5, '(' + str(tr_lon) + ',' + str(tr_lat) + ')', \
            ha='right', va='bottom')
    text(0, n+5, '(' + str(bl_lon) + ',' + str(bl_lat) + ')', \
            ha='left', va='top')
    text(m, n+5, '(' + str(br_lon) + ',' + str(br_lat) + ')', \
            ha='right', va='top')
    colorbar(mappable=img)
    savefig(SAVE_DIR + file_name[0:22] + '.png')
    close(fig)

# Use to check for and get coordinates of the common area of a group of HDF files
def check_common_area(HDF_list):
    lat, lon = read_HDF_SDS(HDF_list, 'Latitude', 'Longitude')
    lat, lon, dummy = correct_inversion(lat, lon)
    del dummy
    
    N = len(lat)
    top_lat = zeros((N,1))
    bottom_lat = zeros((N,1))
    left_lon = zeros((N,1))
    right_lon = zeros((N,1))
    for r in range(N):
        m, n = shape(lat[r])
        top_lat[r] = lat[r][0][0]
        bottom_lat[r] = lat[r][m - 1][n - 1]
        left_lon[r] = lon[r][0][0]
        right_lon[r] = lon[r][m - 1][n - 1]
    
    common_area_top = min(top_lat)
    common_area_bottom = max(bottom_lat)
    common_area_left = max(left_lon)
    common_area_right = min(right_lon)
    
    if common_area_top < common_area_bottom:
        common_area_top = common_area_bottom = NaN
    if common_area_left > common_area_right:
        common_area_left = common_area_right = NaN
    
    return common_area_top, \
                common_area_bottom, \
                common_area_left, \
                common_area_right

# Use to intersect area bounded by coordinates to Philippine area
def intersect_phil_area(top, bottom, left, right):
    intersection_top = min(top, 22.5)
    intersection_bottom = max(bottom, 2.5)
    intersection_left = max(left, 115.)
    intersection_right = min(right, 130.)
    
    if intersection_top < intersection_bottom:
        intersection_top = intersection_bottom = NaN
    if intersection_left > intersection_right:
        intersection_left = intersection_right = NaN
    
    return intersection_top, \
                intersection_bottom, \
                intersection_left, \
                intersection_right

# Use to group HDF files into groups with roughly the same area
def group_HDF_files(file_list):
    input_file = file(file_list, 'r')
    HDF_file_names = input_file.readlines()
    input_file.close()
    
    HDF_groups = []
    
    for HDF_file_name in HDF_file_names:
        HDF_file_satellite = HDF_file_name[0:3]
        HDF_file_year = HDF_file_name[10:14]
        HDF_file_day = HDF_file_name[14:17]
        HDF_file_time = HDF_file_name[18:22]
        
        new_HDF_group = True
        for HDF_group in HDF_groups:
            HDF_group_satellite = HDF_group[0][0:3]
            HDF_group_year = HDF_group[0][10:14]
            HDF_group_day = HDF_group[0][14:17]
            HDF_group_time = HDF_group[0][18:22]
            
            if HDF_file_satellite == HDF_group_satellite \
            and HDF_file_year == HDF_group_year \
            and (abs(int(HDF_file_day) - int(HDF_group_day)) % 16) == 0 \
            and HDF_file_time == HDF_group_time:
                HDF_group.append(HDF_file_name)
                new_HDF_group = False
                break
        
        if new_HDF_group:
            HDF_groups.append([HDF_file_name])
    
    for HDF_group in HDF_groups:
        HDF_group_satellite = HDF_group[0][0:3]
        HDF_group_year = HDF_group[0][10:14]
        HDF_group_day = HDF_group[0][14:17]
        HDF_group_time = HDF_group[0][18:22]
        HDF_group_name = HDF_group_satellite + '-'\
                        + HDF_group_year + '-'\
                        + HDF_group_day + '-'\
                        + HDF_group_time
        output_file = file(LIST_DIR + HDF_group_name, 'w')
        output_file.writelines(HDF_group)
        output_file.close()

# Use to generate cloud fraction images and values
def get_cloud_fraction(HDF_list, top, bottom, left, right):
    lat, lon, cloud_mask = read_HDF_SDS(HDF_list, \
                                                                'Latitude', \
                                                                'Longitude', \
                                                                'Cloud_Mask')
    lat, lon, cloud_mask = correct_inversion(lat, lon, cloud_mask)
    
    input_file = file(LIST_DIR + HDF_list, 'r')
    HDF_file_names = input_file.readlines()
    input_file.close()
    N = len(lat)
    print '\t' + str(N) + ' HDF files in group.'
    
    print '\tBitmasking cloud_mask matrices...'
    cloudiness = bitmask_cloud_mask(cloud_mask, 6)
    land_water = bitmask_cloud_mask(cloud_mask, 192)
    boundary_mask = uint8(cloud_mask[0][0]) & 192
    boundary_mask[boundary_mask==192] = 0  # land
    boundary_mask[boundary_mask==128] = 0  # desert
    boundary_mask[boundary_mask==64] = 1  # coastal
    boundary_mask[boundary_mask==0] = 0  # water
    del cloud_mask
    
    for r in range(N):
        cloudiness[r][cloudiness[r]==6] = 1  # confident clear
        cloudiness[r][cloudiness[r]==4] = 1  # probably clear
        cloudiness[r][cloudiness[r]==2] = 1  # probably cloudy
        cloudiness[r][cloudiness[r]==0] = 100.0 / N  # confident cloudy
        cloudiness[r][cloudiness[r]==1] = 0
        land_water[r][land_water[r]==192] = 1  # land
        land_water[r][land_water[r]==128] = 1  # desert
        land_water[r][land_water[r]==64] = 1  # coastal
        land_water[r][land_water[r]==0] = 0  # water
    
    print '\tIntersecting and stretching matrices...'
    tl_lat = zeros((N,1))  # top left
    tl_lon = zeros((N,1))
    tr_lat = zeros((N,1))  # top right
    tr_lon = zeros((N,1))
    bl_lat = zeros((N,1))  # bottom left
    bl_lon = zeros((N,1))
    br_lat = zeros((N,1))  # bottom right
    br_lon = zeros((N,1))
    x = zeros((N,1), int)  # number of rows
    y = zeros((N,1), int)  # number of columns
    land_cloud_fraction = zeros((N,1))
    water_cloud_fraction = zeros((N,1))
    total_cloud_fraction = zeros((N,1))
    
    for r in range(N):
        print '\t\t' + HDF_file_names[r].rstrip()
        m, n = shape(cloudiness[r])
        lat[r] = expand_matrix(lat[r], m, n)
        lon[r] = expand_matrix(lon[r], m, n)
        I, J, K, L = index_of_corners(lat[r], lon[r], top, left, bottom, right)
        if r == 0:
            boundary_mask = trim_matrix(boundary_mask, I, J, K, L)
        cloudiness[r] = trim_matrix(cloudiness[r], I, J, K, L)
        land_water[r] = trim_matrix(land_water[r], I, J, K, L)
        x[r], y[r] = shape(cloudiness[r])
        tl_lat[r] = lat[r][I][J]
        tl_lon[r] = lon[r][I][J]
        tr_lat[r] = lat[r][I][L]
        tr_lon[r] = lon[r][I][L]
        bl_lat[r] = lat[r][K][J]
        bl_lon[r] = lon[r][K][J]
        br_lat[r] = lat[r][K][L]
        br_lon[r] = lon[r][K][L]
        
        cloud_temp = cloudiness[r] * N
        land_loc = find(land_water[r] > 0)
        water_loc = find(land_water[r] == 0)
        land_cloud_fraction[r] = sum(cloud_temp.flat[land_loc]) / size(land_loc)
        water_cloud_fraction[r] = sum(cloud_temp.flat[water_loc]) / size(water_loc)
        total_cloud_fraction[r] = sum(cloud_temp) / (x[r] * y[r])
    del lat, lon
    
    print '\tCompiling...'
    land_cloud_fraction = round(sum(land_cloud_fraction) / N, 2)
    water_cloud_fraction = round(sum(water_cloud_fraction) / N, 2)
    total_cloud_fraction = round(sum(total_cloud_fraction) / N, 2)
    out_file = file(SAVE_DIR + 'cloud_fraction_values.txt', 'a')
    out_file.write(HDF_list + '-' \
                        + str(land_cloud_fraction) + '-' \
                        + str(water_cloud_fraction) + '-' \
                        + str(total_cloud_fraction) + '\n')
    out_file.close()
    
    print '\tPlotting...'
    tl_lat = round(average(tl_lat), 2)
    tl_lon = round(average(tl_lon), 2)
    tr_lat = round(average(tr_lat), 2)
    tr_lon = round(average(tr_lon), 2)
    bl_lat = round(average(bl_lat), 2)
    bl_lon = round(average(bl_lon), 2)
    br_lat = round(average(br_lat), 2)
    br_lon = round(average(br_lon), 2)
    x = max(x)
    y = max(y)
    boundary_mask = expand_matrix(boundary_mask, x, y)
    for r in range(N):
        cloudiness[r] = expand_matrix(cloudiness[r], x, y)
        if r != 0:
            cloudiness[0] += cloudiness[r]
    cloudiness[0][boundary_mask==1] = 0  # overlay boundaries
    save_HDF_image(cloudiness[0], HDF_list, \
                                tl_lat, tl_lon, tr_lat, tr_lon, \
                                bl_lat, bl_lon, br_lat, br_lon, \
                                cm.jet)

# Use to correct inverted orientations
def correct_inversion(lat, lon, cloud_mask=[]):
    N = len(lat)
    n = len(cloud_mask)
    for r in range(N):
        if lat[r][0][0] < lat[r][1][0]:
            lat[r] = flipud(lat[r])
            lon[r] = flipud(lon[r])
            if n == N:
                cloud_mask[r][0] = flipud(cloud_mask[r][0])
        if lon[r][0][0] > lon[r][0][1]:
            lat[r] = fliplr(lat[r])
            lon[r] = fliplr(lon[r])
            if n == N:
                cloud_mask[r][0] = fliplr(cloud_mask[r][0])
    
    return lat, lon, cloud_mask

# Use to bitmask cloud_mask matrices
def bitmask_cloud_mask(mask_matrices, bit_mask):
    N = len(mask_matrices)
    output_matrices = []
    for r in range(N):
        output_matrices.append(float_(uint8(mask_matrices[r][0]) & bit_mask))
    return output_matrices

# Use to expand m by n matrix to x by y using linear interpolation
def expand_matrix(map_matrix, x, y):
    m, n = shape(map_matrix)
    map_hor_exp = zeros((m, y))
    map_interp = zeros((x, y))
    for i in range(m):
        map_hor_exp[i] = interp(linspace(0,n-1,y), linspace(0,n-1,n), map_matrix[i])
    for j in range(y):
        map_interp[:,j] = interp(linspace(0,m-1,x), linspace(0,m-1,m), map_hor_exp[:,j])
    
    return map_interp

# Use to find index of corners
def index_of_corners(lat, lon, top, left, bottom, right):
# print 'Looking for optimal top left corner...'
    corner_radius = 0.015625
    alg_ceiling = corner_radius
    alg_floor = 0.0
    found_max = False
    intersect_count = 2
    while intersect_count != 1:
        intersect_lat = intersect1d(find((top - corner_radius) < lat), \
                                                    find(lat < top))
        intersect_lon = intersect1d(find(left < lon), \
                                                    find(lon < (left + corner_radius)))
        intersection = intersect1d(intersect_lat, \
                                                    intersect_lon)
        intersect_count = size(intersection)
        if intersect_count > 1:
            corner_radius = corner_radius - (alg_ceiling - alg_floor) / 2.0
            alg_ceiling = corner_radius
            found_max = True
        elif intersect_count == 0:
            if found_max:
                corner_radius = corner_radius + (alg_ceiling - alg_floor) / 2.0
            else:
                corner_radius *= 2.0
            alg_floor = alg_ceiling
            alg_ceiling = corner_radius
    I, J = unravel_index(intersection[0], shape(lat))
# print 'Done.'
    
# print 'Looking for optimal bottom right corner...'
    corner_radius = 0.015625
    alg_ceiling = corner_radius
    alg_floor = 0.0
    found_max = False
    intersect_count = 2
    while intersect_count != 1:
        intersect_lat = intersect1d(find(bottom < lat), \
                                                    find(lat < (bottom + corner_radius)))
        intersect_lon = intersect1d(find((right - corner_radius) < lon), \
                                                    find(lon < right))
        intersection = intersect1d(intersect_lat, \
                                                    intersect_lon)
        intersect_count = size(intersection)
        if intersect_count > 1:
            corner_radius = corner_radius - (alg_ceiling - alg_floor) / 2.0
            alg_ceiling = corner_radius
            found_max = True
        elif intersect_count == 0:
            if found_max:
                corner_radius = corner_radius + (alg_ceiling - alg_floor) / 2.0
            else:
                corner_radius *= 2.0
            alg_floor = alg_ceiling
            alg_ceiling = corner_radius
    K, L = unravel_index(intersection[0], shape(lat))
# print 'Done.'
    
    return I, J, K, L

def trim_matrix(map_matrix, I, J, K, L):
    m,n = shape(map_matrix)
    map_matrix = delete(map_matrix,range(K+1,m),0)
    map_matrix = delete(map_matrix,range(0,I),0)
    map_matrix = delete(map_matrix,range(L+1,n),1)
    map_matrix = delete(map_matrix,range(0,J),1)
    
    return map_matrix
