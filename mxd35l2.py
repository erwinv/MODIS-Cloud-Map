#!/usr/bin/python2 
import numpy as NP
import matplotlib.pyplot as PLT
import pyhdf.SD as SD
import datetime as DT
import copy as CP
from scipy.interpolate import RectBivariateSpline

class MXD35L2:
    """Base class for MOD35_L2 (Terra) or MYD35_L2 (Aqua) products"""
    def __init__(self):
        """Initialize all data to default values"""
        self.name = "noname"
        self.datetimestamp = DT.datetime(2000, 1, 1, 0, 0) # 1 Jan 2000 00:00

        # Initialize all scientific datasets to all zeros
        # with appropriate shape and datatype
        self.lon = NP.zeros((406, 270), dtype=NP.float32)
        self.lat = NP.zeros((406, 270), dtype=NP.float32)
        self.shape = (2030, 1354)
        self.cloud = NP.zeros(self.shape, dtype=NP.uint8)
        self.water = NP.zeros(self.shape, dtype=NP.uint8)
        self.coast = NP.zeros(self.shape, dtype=NP.uint8)

        # Initialize bounding box to top-left and bot-right lon&lat coordinates
        self.top = self.lat[0]
        self.bot = self.lat[-1]
        self.lef = self.lon[0]
        self.ryt = self.lon[-1]

    def interpLonLat(self):
        """Interpolate lon&lat to match the size of cloud/water/coast"""
        print("> Interpolating Lon & Lat...")
        m2, n2 = self.shape
        m1, n1 = self.lon.shape
        i = range(m1)
        j = range(n1)
        ii = NP.linspace(NP.min(i), NP.max(i), m2)
        jj = NP.linspace(NP.min(j), NP.max(j), n2)
        lon = RectBivariateSpline(i, j, self.lon, kx=1, ky=1)(ii, jj)
        lat = RectBivariateSpline(i, j, self.lat, kx=1, ky=1)(ii, jj)
        
        self.lon = lon
        self.lat = lat

    def cutLonLat(self):
        """Cut lon & lat according to specified bounding box"""
        print("> Cutting Lon & Lat...")
        top, bot, lef, ryt = self.top, self.bot, self.lef, self.ryt
        lon = self.lon
        lat = self.lat

        distance_to_toplef = NP.abs((lon + lat*1j) - (lef + top*1j))
        distance_to_botryt = NP.abs((lon + lat*1j) - (ryt + bot*1j))
        T, L = NP.unravel_index(distance_to_toplef.argmin(), self.shape)
        B, R = NP.unravel_index(distance_to_botryt.argmin(), self.shape)
       
        self.lon = cut(self.lon, T, L, B, R)
        self.lat = cut(self.lat, T, L, B, R)

    def cutCloudWaterCoast(self):
        """Cut cloud, water, and coast according to specified bounding box"""
        print("> Cutting Cloud, Water, & Coast...")
        top, bot, lef, ryt = self.top, self.bot, self.lef, self.ryt
        lon = self.lon
        lat = self.lat

        distance_to_toplef = NP.abs((lon + lat*1j) - (lef + top*1j))
        distance_to_botryt = NP.abs((lon + lat*1j) - (ryt + bot*1j))
        T, L = NP.unravel_index(distance_to_toplef.argmin(), self.shape)
        B, R = NP.unravel_index(distance_to_botryt.argmin(), self.shape)
       
        self.cloud = cut(self.cloud, T, L, B, R)
        self.water = cut(self.water, T, L, B, R)
        self.coast = cut(self.coast, T, L, B, R)

    def interpCloudWaterCoast(self):
        """Interpolate cloud, water, and coast to new shape"""
        print("> Interpolating Cloud, Water, & Coast...")
        m2, n2 = self.shape
        m1, n1 = self.cloud.shape
        i = range(m1)
        j = range(n1)
        ii = NP.linspace(NP.min(i), NP.max(i), m2)
        jj = NP.linspace(NP.min(j), NP.max(j), n2)
        self.cloud = RectBivariateSpline(i, j, self.cloud, kx=1, ky=1)(ii, jj)
        self.water = RectBivariateSpline(i, j, self.water, kx=1, ky=1)(ii, jj)
        self.coast = RectBivariateSpline(i, j, self.coast, kx=1, ky=1)(ii, jj)

    def imshowCloudCoast(self, ext='', clrmap=PLT.cm.hot):
        """Show Cloud Fraction image while also showing coast lines"""
        fname = self.name + ext + '.png'
        print("> Plotting: " + fname)
        cloud = CP.deepcopy(self.cloud)

        if clrmap == PLT.cm.jet:
            cloud = cloud * 100
        else:
            cloud = cloud * 50
            cloud = cloud + 50
        cloud[self.coast >= 0.5] = 0

        fig = PLT.figure()
        PLT.imshow(cloud, cmap=clrmap)
        PLT.xticks([0, cloud.shape[1]], [self.lef, self.ryt])
        PLT.yticks([0, cloud.shape[0]], [self.top, self.bot])
        PLT.title(self.datetimestamp.strftime('%G %B %d %R')) # linux
        #PLT.title(self.datetimestamp.strftime('%Y %B %d %I:%M %p')) # windows
        if(clrmap == PLT.cm.jet):
            PLT.colorbar()

        PLT.savefig(fname, dpi=300)
        PLT.close(fig)
    
    def computeCloudFrac(self, ext=''):
        """Compute Cloud Fraction over Water, Land, and Total"""
        fname = self.name + '-CF.txt'
        print("> Computing Cloud Fraction: " + fname)
        totalcloud = self.cloud
        self.totalcloudfrac = totalcloud.sum() / totalcloud.size

        watercloud = CP.deepcopy(self.cloud)
        watercloud[self.water >= 0.5] = 0
        self.watercloudfrac = watercloud.sum() / watercloud.size

        landcloud = CP.deepcopy(self.cloud)
        landcloud[self.water <= 0.5] = 0
        self.landcloudfrac = landcloud.sum() / landcloud.size

        out_text = ext + ' '\
                + str(self.totalcloudfrac) + ' '\
                + str(self.landcloudfrac) + ' '\
                + str(self.watercloudfrac) + '\n'

        out_file = open(fname, 'a')
        out_file.write(out_text)
        out_file.close()

class MXD35L2File(MXD35L2):
    """Derived class that represents a single HDF file's data"""
    def __init__(self, fname):
        """Initialize the HDF file"""
        MXD35L2.__init__(self)
        print("> Reading: " + fname)
        name = fname

        start = fname.find('D35_L2.A') + len('D35_L2.A')
        end = start + len('YYYYDDD.HHMM')
        datetimestamp = DT.datetime.strptime(fname[start:end], '%Y%j.%H%M')
        
        hdf_file = SD.SD(fname)
        lon = hdf_file.select('Longitude').get()
        lat = hdf_file.select('Latitude').get()
        cloud_mask = NP.uint8(hdf_file.select('Cloud_Mask').get()[0])
        hdf_file.end()
        
        cloud = cloud_mask & 6 # get bits 1 and 2
        cloud[cloud == 0] = 1 # 00 = confident cloudy
        cloud[cloud != 1] = 0
        
        water = cloud_mask & 192 # get bits 6 and 7
        water[water == 0] = 1 # 00 = water
        water[water != 1] = 0
        
        coast = cloud_mask & 192 # get bits 6 and 7
        coast[coast == 64] = 1 # 01 = coastal
        coast[coast != 1] = 0
        
        lon, lat, cloud, water, coast = autoFlip(lon, 'horizontal', lon, lat,
                cloud, water, coast)
        lon, lat, cloud, water, coast = autoFlip(lat, 'vertical', lon, lat,
                cloud, water, coast)
        
        self.name = name
        self.datetimestamp = datetimestamp
        self.lon = lon
        self.lat = lat
        self.cloud = cloud
        self.water = water
        self.coast = coast
        self.top = self.lat[0]
        self.bot = self.lat[-1]
        self.lef = self.lon[0]
        self.ryt = self.lon[-1]
        
class MXD35L2Group(MXD35L2File):
    """Derived class that represents a group of HDF files"""
    def __init__(self, *fnames):
        """Initialize the HDF group using the first HDF file"""
        N = len(fnames)

        #Get the lon&lat coordinates of the common area of the whole group
        tops = []
        bots = []
        lefs = []
        ryts = []

        for fname in fnames:
            hdf_file = SD.SD(fname)
            lon = hdf_file.select('Longitude').get()
            lat = hdf_file.select('Latitude').get()
            hdf_file.end()

            lon, lat = autoFlip(lon, 'horizontal', lon, lat)
            lon, lat = autoFlip(lat, 'vertical', lon, lat)

            tops.append(lat[0][0])
            bots.append(lat[-1][-1])
            lefs.append(lon[0][0])
            ryts.append(lon[-1][-1])

        #Intersect with Philippine area, 2.5-22.5 deg lat, 115-130 deg lon
        top = min(NP.min([tops]), 22.5)
        bot = max(NP.max([bots]), 2.5)
        lef = max(NP.max([lefs]), 115)
        ryt = min(NP.min([ryts]), 130)

        print("Common area (top, bot)(left, right): "
                + '(' + str(top) + ',' + str(bot) + ')'
                + '(' + str(lef) + ',' + str(ryt) + ')')

        if top < bot or lef > ryt:
            self.proceed = False
            print("No common area found")
        else:
            self.proceed = True
            print("Processing: 1/" + str(N))
            MXD35L2File.__init__(self, fnames[0])

            self.top = top
            self.bot = bot
            self.lef = lef
            self.ryt = ryt

            self.interpLonLat()
            self.imshowCloudCoast('-1RAW')
            self.computeCloudFrac('1RAW')
            self.cutCloudWaterCoast()
            self.imshowCloudCoast('-2CUT')
            self.computeCloudFrac('2CUT')
            self.interpCloudWaterCoast()
            self.imshowCloudCoast('-3INTERP')
            self.computeCloudFrac('3INTERP')

        self.N = N
        self.members = fnames

    def processGroup(self):
        """Process the rest of the HDF files in the group"""
        if self.proceed:
            idx = 2
            for fname in self.members[1:]:
                print("Processing: " + str(idx) + "/" + str(self.N))
                newmember = MXD35L2File(fname)
                newmember.top = self.top
                newmember.bot = self.bot
                newmember.lef = self.lef
                newmember.ryt = self.ryt
                newmember.interpLonLat()
                newmember.imshowCloudCoast('-1RAW')
                newmember.computeCloudFrac('1RAW')
                newmember.cutCloudWaterCoast()
                newmember.imshowCloudCoast('-2CUT')
                newmember.computeCloudFrac('2CUT')
                newmember.interpCloudWaterCoast()
                newmember.imshowCloudCoast('-3INTERP')
                newmember.computeCloudFrac('3INTERP')
                
                self.cloud = self.cloud + newmember.cloud
                self.water = self.water + newmember.water
                self.coast = self.coast + newmember.coast
                
                idx = idx + 1

            self.cloud = self.cloud / self.N
            self.water = self.water / self.N
            self.coast = self.coast / self.N
            
            return True
        else:
            return False

def autoFlip(reference, orientation, *datasets):
    """Automatically flips datasets to upright orientation"""
    if orientation == 'horizontal':
        if reference[0][0] > reference[0][-1]:
            datasets = [NP.fliplr(dataset) for dataset in datasets]
    elif orientation == 'vertical':
        if reference[0][0] < reference[-1][0]:
            datasets = [NP.flipud(dataset) for dataset in datasets]
    else:
        raise("Invalid orientation: use either 'horizontal' or 'vertical'")
    return datasets

def cut(dataset, top, lef, bot, ryt):
    m, n = dataset.shape
    dataset = NP.delete(dataset, range(bot + 1, m), 0)
    dataset = NP.delete(dataset, range(0, top), 0)
    dataset = NP.delete(dataset, range(ryt + 1, n), 1)
    dataset = NP.delete(dataset, range(0, lef), 1)
    return dataset

if __name__ == '__main__':
    import sys
    if sys.argv[1] == '--plot-files':
        for arg in sys.argv[2:]:
            MyMXD35L2File = MXD35L2File(arg)
            MyMXD35L2File.imshowCloudCoast()
    elif sys.argv[1] == '--process-group':
        MyMXD35L2Group = MXD35L2Group(*sys.argv[2:])
        if MyMXD35L2Group.processGroup():
            print("Done.")
            MyMXD35L2Group.imshowCloudCoast('-4GROUP', PLT.cm.jet)
            MyMXD35L2Group.computeCloudFrac('4GROUP')
    elif sys.argv[1] == '--input-file':
        in_file = open(sys.argv[2], 'rt')
        fnames = in_file.read().split()
        in_file.close()
        MyMXD35L2Group = MXD35L2Group(*fnames)
        MyMXD35L2Group.processGroup()
        print("Done.")
        MyMXD35L2Group.imshowCloudCoast('-4GROUP', PLT.cm.jet)
        MyMXD35L2Group.computeCloudFrac('4GROUP')
    elif sys.argv[1] == '--help':
        print("Valid options: --plot-files, --process-group, --input-file")
    else:
        print("Invalid option: " + sys.argv[1])
        print("Pass --help option to get all valid options")
