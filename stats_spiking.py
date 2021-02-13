

import common_functions                 as cf
import numpy                            as np





modulation = 1

cell_type = 'dspm'
mod_type = 'ACh'


if not modulation:
    
    delta = list(np.arange(0,100+1,20))
    
    # load data
    spike_data = {}
    for i, r in enumerate(delta):
        data = cf.load_data('Data/{}_HFI[1]+{}_validation.json'.format(cell_type, delta[i]))
        spike_data[r] = data['all']
        
    clus_info = data['meta']['clustered']
    clus_labels = clus_info['label']
    
    # collects iteration-wise spike data into single vector (control data)
    spiked = {}
    for delt in delta:
        spiked[delt] = {}
        for clus in clus_labels:
            spiked[delt][clus] = []
            for i in range(len(spike_data[delt][clus]['spiked'])):
                for j in range(len(spike_data[delt][clus]['spiked'][i])):
                    spiked[delt][clus].append(spike_data[delt][clus]['spiked'][i][j])
    
    
    # normality checking =====
    
    # gets stats
    spike_data['stats'] = {}
    
    # tests on data =====
    
    spike_data['stats']['diff'] = {'spiked':{}}
    
    for delt in delta:        
        # tests whether number of cells spiking significantly different
        #tab = cf.binary_to_table(spike_data[delt][clus_labels[0]]['spiked'],spike_data[delt][clus_labels[1]]['spiked'])
        tab = cf.binary_to_table(spiked[delt][clus_labels[0]],spiked[delt][clus_labels[1]])
        spike_data['stats']['diff']['spiked'][delt] = cf.McNemar(tab)
        


else:
    
    
    delta = list(np.arange(0,100+1,20))
    
    # load data
    spike_data = {}
    spike_data_ctrl = {}
    for delt in delta:
        data = cf.load_data('Data/{}_HFI[1]+{}_validation.json'.format(cell_type, delt))
        spike_data_ctrl[delt] = data['all']
        
        data = cf.load_data('Data/{}_HFI[1]+{}_{}-modulation.json'.format(cell_type, delt, mod_type))
        spike_data[delt] = data['all']
    
    clus_info = data['meta']['clustered']
    clus_labels = clus_info['label']
    
    mod_info = data['meta'][mod_type + ' info']
    mod_labels = mod_info['label']
    
    # collects iteration-wise spike data into single vector (modulation data)
    spiked = {}
    for delt in delta:
        spiked[delt] = {}
        for clus in clus_labels:
            spiked[delt][clus] = {}
            spiked[delt][clus][clus] = []
            for i in range(len(spike_data[delt][clus][clus]['spiked'])):
                for j in range(len(spike_data[delt][clus][clus]['spiked'][i])):
                    spiked[delt][clus][clus].append(spike_data[delt][clus][clus]['spiked'][i][j])
            for mod in mod_labels:
                spiked[delt][clus][mod] = []
                for i in range(len(spike_data[delt][clus][mod]['spiked'])):
                    for j in range(len(spike_data[delt][clus][mod]['spiked'][i])):
                        spiked[delt][clus][mod].append(spike_data[delt][clus][mod]['spiked'][i][j])
    
    # collects iteration-wise spike data into single vector (control data)
    spiked_ctrl = {}
    for delt in delta:
        spiked_ctrl[delt] = {}
        for clus in clus_labels:
            spiked_ctrl[delt][clus] = []
            for i in range(len(spike_data_ctrl[delt][clus]['spiked'])):
                for j in range(len(spike_data_ctrl[delt][clus]['spiked'][i])):
                    spiked_ctrl[delt][clus].append(spike_data_ctrl[delt][clus]['spiked'][i][j])
                        
         
        
    # tests on data =====
    
    spike_data['stats'] = {}       
    
    # tests whether distal and proximal stimulation data is significantly different to each other
    spike_data['stats']['diff'] = {}
    for d, delt in enumerate(delta):
        if d == 0:
            spike_data['stats']['diff']['on-site'] = {}
        
        # whether spike occured
        #tab = cf.binary_to_table(spike_data[delt][clus_labels[0]][clus_labels[0]]['spiked'], spike_data[delt][clus_labels[1]][clus_labels[1]]['spiked'])
        tab = cf.binary_to_table(spiked[delt][clus_labels[0]][clus_labels[0]], spiked[delt][clus_labels[1]][clus_labels[1]])
        spike_data['stats']['diff']['on-site'][delt] = cf.McNemar(tab)
        
        for j, mod_lab in enumerate(mod_labels):
            if d == 0:
                spike_data['stats']['diff'][mod_lab] = {}
            
            # whether spike occured
            #tab = cf.binary_to_table(spike_data[delt][clus_labels[0]][ACh_lab]['spiked'], spike_data[delt][clus_labels[1]][ACh_lab]['spiked'])
            tab = cf.binary_to_table(spiked[delt][clus_labels[0]][mod_lab], spiked[delt][clus_labels[1]][mod_lab])
            spike_data['stats']['diff'][mod_lab][delt] = cf.McNemar(tab)
                
                
        
    
    # tests whether modulation data is significantly different to control
    spike_data['stats']['diff_ctrl'] = {}
    for d, delt in enumerate(delta):
        
        for i, clus_lab in enumerate(clus_labels):
            if d == 0:
                spike_data['stats']['diff_ctrl'][clus_lab]= {}
                spike_data['stats']['diff_ctrl'][clus_lab][clus_lab] = {}
        
            # whether spike occured
            #tab = cf.binary_to_table(spike_data_ctrl[d][clus_lab]['spiked'], spike_data[d][clus_lab][clus_lab]['spiked'])
            tab = cf.binary_to_table(spiked_ctrl[delt][clus_lab], spiked[delt][clus_lab][clus_lab])
            spike_data['stats']['diff_ctrl'][clus_lab][clus_lab][delt] = cf.McNemar(tab)
        
            for j, mod_lab in enumerate(mod_labels):
                if d == 0:
                    spike_data['stats']['diff_ctrl'][clus_lab][mod_lab] = {}
                
                # whether spike occured
                #tab = cf.binary_to_table(spike_data_ctrl[delt][clus_lab]['spiked'], spike_data[delt][clus_lab][ACh_lab]['spiked'])
                tab = cf.binary_to_table(spiked_ctrl[delt][clus_lab], spiked[delt][clus_lab][mod_lab])
                spike_data['stats']['diff_ctrl'][clus_lab][mod_lab][delt] = cf.McNemar(tab)
        
    
    
