

import common_functions                 as cf
import scipy.stats                      as stats
import matplotlib.pyplot                as plt
import numpy                            as np





modulation = 0



if not modulation:
    
    delta = list(np.arange(0,100+1,10))
    
    # load data
    spike_data = {}
    for i, r in enumerate(delta):
        data = cf.load_data('Data/ispn_HFI[1]+{}_validation.json'.format(delta[i]))
        spike_data[r] = data['all']
        
    clus_info = data['meta']['clustered']
    clus_labels = clus_info['label']
    cell_type = data['meta']['cell type']
    


    # normality checking =====
    
    # gets stats
    spike_data['stats'] = {}
    '''
    for r in delta:
        spike_data['stats'][r] = {}
        
        for clus_lab in clus_labels:
            
            spike_data[r][clus_lab]['first_spike'] = [x for x in spike_data[r][clus_lab]['first_spike'] if x]
            spike_data[r][clus_lab]['spike_n'] = [x for x in spike_data[r][clus_lab]['spike_n'] if x]
            
            spike_data['stats'][r][clus_lab] = {}
            # first spike time data
            spike_data['stats'][r][clus_lab]['first_spike'] = cf.norm_dist(spike_data[r][clus_lab]['first_spike'])
            # num spikes data
            spike_data['stats'][r][clus_lab]['spike_n'] = cf.norm_dist(spike_data[r][clus_lab]['spike_n'])
            
        
    # plot histograms
    for i, clus_lab in enumerate(clus_labels):
        
        fig, axs = plt.subplots(2,len(delta))
        fig.suptitle(cell_type + ', ' + clus_lab)
        
        for d, delt in enumerate(delta):
            
            # first spike time data
            axs[0,d].hist(spike_data[delt][clus_lab]['first_spike'])
            axs[0,d].set_title('time,{},{}'.format(delt, str(spike_data['stats'][delt][clus_lab]['first_spike'])))
            
            # num spikes data
            axs[1,d].hist(spike_data[delt][clus_lab]['spike_n'])
            axs[1,d].set_title('num,{},{}'.format(delt, str(spike_data['stats'][delt][clus_lab]['spike_n'])))
        
        plt.tight_layout()
    '''
    
    # tests on data =====
    
    spike_data['stats']['diff'] = {'spiked':{}, 'first_spike':{}, 'spike_n':{}}
    
    for delt in delta:
        
        # tests whether number of cells spiking significantly different
        tab = cf.binary_to_table(spike_data[delt][clus_labels[0]]['spiked'],spike_data[delt][clus_labels[1]]['spiked'])
        spike_data['stats']['diff']['spiked'][delt] = cf.McNemar(tab)
        '''
        # tests whether time to first spike significantly different
        spike_data['stats']['diff']['first_spike'][delt] = test.test(spike_data[delt][clus_labels[0]]['first_spike'], \
                  spike_data[delt][clus_labels[1]]['first spike'])
            
        # tests whether number of spikes significantly different
        spike_data['stats']['diff']['spike_n'][delt] = test.test(spike_data[delt][clus_labels[0]]['spike_n'], \
                  spike_data[delt][clus_labels[1]]['spike_n'])
        '''
        

else:
    
    delta = list(np.arange(0,100+1,20))
    
    # load data
    spike_data = {}
    spike_data_ctrl = {}
    for i, r in enumerate(delta):
        data = cf.load_data('Data/dspn_HFI[1]+{}_validation.json'.format(delta[i]))
        spike_data_ctrl[r] = data['all']
        
        data = cf.load_data('Data/dspn_HFI[1]+{}_modulation.json'.format(delta[i]))
        spike_data[r] = data['all']
    
    clus_info = data['meta']['clustered']
    clus_labels = clus_info['label']
    cell_type = data['meta']['cell type']
    
    ACh_info = data['meta']['ACh info']
    ACh_labels = ACh_info['label']


    # normality checking =====
    
    # gets stats
    spike_data['stats'] = {}
    for r in delta:
        spike_data['stats'][r] = {}
        
        for clus_lab in clus_labels:
            
            spike_data['stats'][r][clus_lab] = {}
            spike_data['stats'][r][clus_lab][clus_lab] = {}
            # first spike time data
            spike_data['stats'][r][clus_lab][clus_lab]['first_spike'] = cf.norm_dist( \
                      spike_data[r][clus_lab][clus_lab]['first_spike'])
            # num spikes data
            spike_data['stats'][r][clus_lab][clus_lab]['spike_n'] = cf.norm_dist( \
                      spike_data[r][clus_lab][clus_lab]['spike_n'])
            
            for ACh_lab in ACh_labels:
        
                spike_data['stats'][r][clus_lab][ACh_lab] = {}
                # first spike time data
                spike_data['stats'][r][clus_lab][ACh_lab]['first_spike'] = cf.norm_dist( \
                          spike_data[r][clus_lab][ACh_lab]['first_spike'])
                # num spikes data
                spike_data['stats'][r][clus_lab][ACh_lab]['spike_n'] = cf.norm_dist( \
                          spike_data[r][clus_lab][ACh_lab]['spike_n'])
    
    
    
    # plot histograms
    for d, delt in enumerate(delta):
        
        for i, clus_lab in enumerate(clus_labels):
        
            fig, axs = plt.subplots(2,len(delta))
            fig.suptitle(cell_type + ', {}, {}'.format(clus_lab,ACh_lab))
            
            # first spike time data
            axs[0,0].hist(spike_data[r][clus_lab][clus_lab]['first_spike'])
            axs[0,0].set_title('time, +{}, {}'.format(delt, str(data['stats'][d][clus_lab][clus_lab]['first_spike'])))
            
            # num spikes data
            axs[1,0].hist(spike_data[r][clus_lab][clus_lab]['spike_n'])
            axs[1,0].set_title('num, +{}, {}'.format(delt, str(data['stats'][d][clus_lab][clus_lab]['spike_n'])))
        
            for j, ACh_lab in enumerate(ACh_labels):
                
                # first spike time data
                axs[0,j+1].hist(spike_data[r][clus_lab][ACh_lab]['first_spike'])
                axs[0,j+1].set_title('time, +{}, {}'.format(delt, str(data['stats'][d][clus_lab][ACh_lab]['first_spike'])))
                
                # num spikes data
                axs[1,j+1].hist(spike_data[r][clus_lab][ACh_lab]['spike_n'])
                axs[1,j+1].set_title('num, +{}, {}'.format(delt, str(data['stats'][d][clus_lab][ACh_lab]['spike_n'])))
            
            plt.tight_layout()
            
        
    # tests on data =====
    
    # tests whether distal and proximal stimulation data is significantly different to each other
    for d, delt in enumerate(delta):
        
        spike_data['stats'][d]['diff'] = {}
        spike_data['stats'][d]['diff']['on site'] = {}
        
        # whether spike occured
        tab = cf.binary_to_table(spike_data[d][clus_labels[0]][clus_labels[0]]['spiked'], \
                                 spike_data[d][clus_labels[1]][clus_labels[1]]['spiked'])
        spike_data['stats'][d]['diff']['on site']['spiked'] = cf.McNemar(tab)
        # time of first spike
        spike_data['stats'][d]['diff']['on site']['first_spike'] = test.test( \
                  spike_data[d][clus_labels[0]][clus_labels[0]]['first_spike'], \
                  spike_data[d][clus_labels[1]][clus_labels[1]]['first_spike'])
        # number of spikes
        spike_data['stats'][d]['diff']['on site']['spike_n'] = test.test( \
                  spike_data[d][clus_labels[0]][clus_labels[0]]['spike_n'], \
                  spike_data[d][clus_labels[1]][clus_labels[1]]['spike_n'])
        
        for j, ACh_lab in enumerate(ACh_labels):
            
            spike_data['stats'][d]['diff'][ACh_lab] = {}
            
            # whether spike occured
            tab = cf.binary_to_table(spike_data[d][clus_labels[0]][ACh_lab]['spiked'], \
                                     spike_data[d][clus_labels[1]][ACh_lab]['spiked'])
            spike_data['stats'][d]['diff'][ACh_lab]['spiked'] = cf.McNemar(tab)
            # time of first spike
            spike_data['stats'][d]['diff'][ACh_lab]['first_spike'] = test.test( \
                      spike_data[d][clus_labels[0]][ACh_lab]['first_spike'], \
                      spike_data[d][clus_labels[1]][ACh_lab]['first_spike'])
            # number of spikes
            spike_data['stats'][d]['diff'][ACh_lab]['spike_n'] = test.test( \
                      spike_data[d][clus_labels[0]][ACh_lab]['spike_n'], \
                      spike_data[d][clus_labels[1]][ACh_lab]['spike_n'])
                
                
        
    
    # tests whether modulation data is significantly different to control
    for d, delt in enumerate(delta):
        
        spike_data['stats'][d]['diff_ctrl'] = {}
        
        for i, clus_lab in enumerate(clus_labels):
            
            spike_data['stats'][d]['diff_ctrl'][clus_lab]= {}
            spike_data['stats'][d]['diff_ctrl'][clus_lab][clus_lab] = {}
        
            # whether spike occured
            tab = cf.binary_to_table(spike_data_ctrl[d][clus_lab]['spiked'], \
                                     spike_data[d][clus_lab][clus_lab]['spiked'])
            spike_data['stats'][d]['diff_ctrl'][clus_lab][clus_lab]['spiked'] = cf.McNemar(tab)
            # time of first spike
            spike_data['stats'][d]['diff'][clus_lab][clus_lab]['first_spike'] = test.test( \
                      spike_data_ctrl[d][clus_lab]['first_spike'], spike_data[d][clus_lab][clus_lab]['first_spike'])
            # number of spikes
            spike_data['stats'][d]['diff'][clus_lab][clus_lab]['spike_n'] = test.test( \
                      spike_data_ctrl[d][clus_lab]['spike_n'], spike_data[d][clus_lab][clus_lab]['spike_n'])
        
            for j, ACh_lab in enumerate(ACh_labels):
                
                spike_data['stats'][d]['diff_ctrl'][clus_lab][ACh_lab] = {}
                
                # whether spike occured
                tab = cf.binary_to_table(spike_data_ctrl[d][clus_lab]['spiked'], \
                                         spike_data[d][clus_lab][ACh_lab]['spiked'])
                spike_data['stats'][d]['diff_ctrl'][clus_lab][ACh_lab]['spiked'] = cf.McNemar(tab)
                # time of first spike
                spike_data['stats'][d]['diff'][clus_lab][ACh_lab]['first_spike'] = test.test( \
                          spike_data_ctrl[d][clus_lab]['first_spike'], spike_data[d][clus_lab][ACh_lab]['first_spike'])
                # number of spikes
                spike_data['stats'][d]['diff'][clus_lab][ACh_lab]['spike_n'] = test.test( \
                          spike_data_ctrl[d][clus_lab]['spike_n'], spike_data[d][clus_lab][ACh_lab]['spike_n'])
        
    
    
