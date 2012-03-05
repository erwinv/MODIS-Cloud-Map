SAVE_DIR = '/home/erwin/plots/'

from pylab import *

def plot_cloud_fraction(cloud_fraction_group):
    N = len(cloud_fraction_group)
    
    land_cfrac = zeros((N,1))
    water_cfrac = zeros((N,1))
    total_cfrac = zeros((N,1))
    cfrac_years = []
    
    for i in range(N):
        land_cfrac[i] = cloud_fraction_group[i][4]
        water_cfrac[i] = cloud_fraction_group[i][5]
        total_cfrac[i] = cloud_fraction_group[i][6]
        cfrac_years.append(cloud_fraction_group[i][1])
    
    fig = figure()
    cfrac_x_axis = range(N)
    cfrac_plot = plot(cfrac_x_axis, land_cfrac, 'ro-', \
                      cfrac_x_axis, water_cfrac, 'bo-', \
                      cfrac_x_axis, total_cfrac, 'go-')
    xticks(cfrac_x_axis, cfrac_years, rotation=20)
    xlabel('Year')
    ylabel('Cloud Fraction (%)')
    grid(True)
    
    acq_satellite = cloud_fraction_group[0][0]
    acq_day = cloud_fraction_group[0][2]
    acq_time = cloud_fraction_group[0][3]
    title(acq_satellite + ' ' + \
          'Group ' + acq_day + \
          ' ' + acq_time + ' Hours')
    legend(('Over Land', 'Over Water', 'Total'), \
            'lower right', shadow=True, fancybox=True)
    
    plot_name = acq_satellite + '-' + \
                acq_day + '-' + \
                acq_time + '.png'
    savefig(SAVE_DIR + plot_name)
    close(fig)

