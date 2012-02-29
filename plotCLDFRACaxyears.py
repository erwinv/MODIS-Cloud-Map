SAVE_DIR = '/home/erwin/gits/crap/save/plots/final/'

from pylab import *

def plot_cfracs_across_years(cloud_fraction_groups, cfrac='total'):
    num_of_groups = 0

    for every_group in cloud_fraction_groups:
        if len(every_group) == 9:
            num_of_groups += 1

    total_cfracs = zeros((num_of_groups, 9))
    group_names = []
    idx = 0

    for every_group in cloud_fraction_groups:
        if len(every_group) == 9:
            for i in range(9):
                if cfrac == 'land':
                    j = 4
                elif cfrac == 'water':
                    j = 5
                elif cfrac == 'total':
                    j = 6
                else:
                    j = 6
                total_cfracs[idx][i] = every_group[i][j]
            idx += 1
            group_satellite = every_group[0][0]
            group_day = every_group[0][2]
            group_time = every_group[0][3]
            group_names.append(group_satellite + '-' + \
                                                group_day + '-' + \
                                                group_time)
    
    fig = figure()
    x = range(9)
    for i in range(num_of_groups):
        if i < 7:
            line_pattern = 'x-'
        elif i < 14:
            line_pattern = '+-.'
        elif i < 21:
            line_pattern = '*-'
        else:
            line_pattern = 'o--'
        plot(x, total_cfracs[i], line_pattern)
        
    grid(True)
    xticks(x, range(2002, 2011), rotation=20)
    xlabel('Year')
    ylabel('Cloud Fraction (%)')
    
    if cfrac == 'land':
        plot_title = 'Land Cloud Fractions'
    elif cfrac == 'water':
        plot_title = 'Water Cloud Fractions'
    elif cfrac == 'total':
        plot_title = 'Total Cloud Fractions'
    else:
        plot_title = 'Total Cloud Fractions'
    
    title(plot_title)
    #legend((group_names), 'lower right',
    #            shadow=True, fancybox=True)
    
    if cfrac == 'land':
        save_name = 'land_cfracs.png'
    elif cfrac == 'water':
        save_name = 'water_cfracs.png'
    elif cfrac == 'total':
        save_name = 'total_cfracs.png'
    else:
        save_name = 'total_cfracs.png'
    
    savefig(SAVE_DIR + save_name)
    close(fig)