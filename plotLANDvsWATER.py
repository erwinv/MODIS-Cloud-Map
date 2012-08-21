SAVE_DIR = '/home/erwin/plots/'

from pylab import *

def plot_cf_land_water(cloud_fraction_values):
    N = len(cloud_fraction_values)
    
    land_cfrac = zeros((N,1))
    water_cfrac = zeros((N,1))
    
    for i in range(N):
        land_cfrac[i] = cloud_fraction_values[i][4]
        water_cfrac[i] = cloud_fraction_values[i][5]
    
    fig = figure()
    cfrac_plot = plot(land_cfrac, water_cfrac, 'o',
		      range(101), range(101), '-')
    xlabel('Cloud Fraction over Land (%)')
    ylabel('Cloud Fraction over Water (%)')
    grid(True)
    
    title('CF_Land vs CF_Water')
    
    plot_name = 'landvswater.png'
    savefig(SAVE_DIR + plot_name)
    close(fig)

