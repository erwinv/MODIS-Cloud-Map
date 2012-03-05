from pylab import *

def get_cloud_fraction_groups(my_list):
    list_of_groups = []

    for every_file in my_list:
        acq_satellite = every_file[0]
        acq_year = every_file[1]
        acq_day = every_file[2]
        acq_time= every_file[3]
        
        new_group = True
        
        for every_group in list_of_groups:
            group_satellite = every_group[0][0]
            group_year = every_group[0][1]
            group_day = every_group[0][2]
            group_time= every_group[0][3]
            
            if acq_satellite == group_satellite and \
                acq_day == group_day and \
                acq_time == group_time:
                    every_group.append(every_file)
                    new_group = False
                    break;
        
        if new_group:
            list_of_groups.append([every_file])

    return list_of_groups

