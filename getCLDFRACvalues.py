from pylab import *

def get_cloud_fraction_values(file_name):
    input_file = file(file_name, 'r')
    text_file_readlines = input_file.readlines()
    input_file.close()
    del text_file_readlines[0]
    
    cloud_fraction_values = []

    for every_line in text_file_readlines:
        # Identifiers
        acq_satellite = every_line[0:3]
        acq_year = every_line[4:8]
        acq_day = every_line[9:12]
        acq_time = every_line[13:17]
        
        # Cloud fraction values
        cloud_fraction = [[], [], []]
        idx = 0
        init = 17
        for i in range(18, 17 + len(every_line[17:])):
            if every_line[i] == '-' or every_line[i] == '\n':
                cloud_fraction[idx] = float(every_line[init+1:i])
                idx += 1
                init = i
        
        cloud_fraction_values.append([acq_satellite, \
                                                        acq_year, \
                                                        acq_day, \
                                                        acq_time, \
                                                        cloud_fraction[0], \
                                                        cloud_fraction[1], \
                                                        cloud_fraction[2]])
    
    return cloud_fraction_values