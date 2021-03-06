

import common_functions                 as cf
import scipy.stats                      as stats
import matplotlib.pyplot                as plt





modulation = 1

cell_type = 'dspn'
mod_type = 'ACh'

if not modulation:
    
    # load data
    data = cf.load_data('Data/{}_HFI[0]+0_validation.json'.format(cell_type))
    clus_info = data['meta']['clustered']
    clus_labels = clus_info['label']
    
    # normality checking =====
    
    # gets stats
    data['stats'] = {}
    for clus_lab in clus_labels:
        
        data['stats'][clus_lab] = {'dur':{}, 'amp':{}}
        # duration data
        data['stats'][clus_lab]['dur'] = cf.norm_dist(data['all'][clus_lab]['dur'])
        # amplitude data
        data['stats'][clus_lab]['amp'] = cf.norm_dist(data['all'][clus_lab]['amp'])
        
    
    # plot histograms
    fig, axs = plt.subplots(2,2)
    fig.suptitle(cell_type)
    for i, clus_lab in enumerate(clus_labels):
        
        # duration data
        axs[i,0].hist(data['all'][clus_lab]['dur'])
        axs[i,0].set_title('dur, ' + clus_lab + ', ' + str(data['stats'][clus_lab]['dur']))
        
        # amplitude data
        axs[i,1].hist(data['all'][clus_lab]['amp'])
        axs[i,1].set_title('amp, ' + clus_lab + ', ' + str(data['stats'][clus_lab]['amp']))
        
        plt.tight_layout()
    
    
    
    # tests on data =====
        # for dspn, dur and amp data was normally distributed; for ispn, a bit iffy but still within reason for normal distribution
    
    
    data['stats']['diff'] = {'dur':{}, 'amp':{}}
    
    # duration
    data['stats']['diff']['dur']['result'] = stats.anderson_ksamp([data['all'][clus_labels[1]]['dur'],
        data['all'][clus_labels[0]]['dur']])[-1]
    data['stats']['diff']['dur']['label'] = 'Anderson-Darling; proximal:distal duration'
    # ttest_rel() doesn't have one-sided implemented, so to see if duration of one greater than other, 
        # t-statistic must be > 0 and p-value/2 < .05
        
    # amplitude
    data['stats']['diff']['amp']['result'] = stats.anderson_ksamp([data['all'][clus_labels[1]]['amp'],
        data['all'][clus_labels[0]]['amp']])[-1]
    data['stats']['diff']['amp']['label'] = 'Anderson-Darling; proximal:distal amplitude'
    
    
else:
    
    
    # load data
    data = cf.load_data('Data/{}_HFI[0]+0_{}-modulation.json'.format(cell_type, mod_type))
    ctrl_data = cf.load_data('Data/{}_HFI[0]+0_validation.json'.format(cell_type))
    clus_info = data['meta']['clustered']
    clus_labels = clus_info['label']
    
    mod_info = data['meta'][mod_type + ' info']
    mod_labels = mod_info['label']
    
    cell_type = data['meta']['cell type']
    
    # normality checking =====
    
    # gets stats
    data['stats'] = {}
    for clus_lab in clus_labels:
        
        data['stats'][clus_lab] = {}
        data['stats'][clus_lab][clus_lab] = {'dur':{}, 'amp':{}}
        # duration data
        data['stats'][clus_lab][clus_lab]['dur'] = cf.norm_dist(data['all'][clus_lab][clus_lab]['dur'])
        # amplitude data
        data['stats'][clus_lab][clus_lab]['amp'] = cf.norm_dist(data['all'][clus_lab][clus_lab]['amp'])
        
        for mod_lab in mod_labels:
            
            data['stats'][clus_lab][mod_lab] = {'dur':{}, 'amp':{}}
            # duration data
            data['stats'][clus_lab][mod_lab]['dur'] = cf.norm_dist(data['all'][clus_lab][mod_lab]['dur'])
            # amplitude data
            data['stats'][clus_lab][mod_lab]['amp'] = cf.norm_dist(data['all'][clus_lab][mod_lab]['amp'])
    
    
    # plot histograms
    for i, clus_lab in enumerate(clus_labels):
        
        fig, axs = plt.subplots(2,len(mod_labels)+1)
        fig.suptitle(cell_type + ', ' + clus_lab)
        
        # duration data
        axs[0,0].hist(data['all'][clus_lab][clus_lab]['dur'])
        axs[0,0].set_title('dur, on-site, ' + str(data['stats'][clus_lab][clus_lab]['dur']))
        
        # amplitude data
        axs[1,0].hist(data['all'][clus_lab][clus_lab]['amp'])
        axs[1,0].set_title('amp, on-site, ' + str(data['stats'][clus_lab][clus_lab]['amp']))
        
        for j, mod_lab in enumerate(mod_labels):
        
            # duration data
            axs[0,j+1].hist(data['all'][clus_lab][mod_lab]['dur'])
            axs[0,j+1].set_title('dur, ' + mod_lab + ', ' + str(data['stats'][clus_lab][mod_lab]['dur']))
            
            # amplitude data
            axs[1,j+1].hist(data['all'][clus_lab][mod_lab]['amp'])
            axs[1,j+1].set_title('amp, ' + mod_lab + ', ' + str(data['stats'][clus_lab][mod_lab]['amp']))
        
        plt.tight_layout()  
        
    
    
    # tests on data =====
        # for dspn and ispn, dur and amp data was pretty normally distributed
    
    data['stats']['diff'] = {}
    
    
    # tests whether duration and amplitude values for each modulation target are signficantly different to control values
    
    for i, clus_lab in enumerate(clus_labels):
        
        data['stats']['diff'][clus_lab] = {}
        data['stats']['diff'][clus_lab][clus_lab] = {'dur':{}, 'amp':{}}
        # duration data
        data['stats']['diff'][clus_lab][clus_lab]['dur']['result'] = stats.anderson_ksamp([ctrl_data['all'][clus_lab]['dur'],
            data['all'][clus_lab][clus_lab]['dur']])[-1]
        data['stats']['diff'][clus_lab][clus_lab]['dur']['label'] = 'Anderson-Darling; control:modulated duration'
        # amplitude data
        data['stats']['diff'][clus_lab][clus_lab]['amp']['result'] = stats.anderson_ksamp([ctrl_data['all'][clus_lab]['amp'],
            data['all'][clus_lab][clus_lab]['amp']])[-1]
        data['stats']['diff'][clus_lab][clus_lab]['amp']['label'] = 'Anderson-Darling; control:modulated amplitude'
        
        for j, mod_lab in enumerate(mod_labels):
            
            data['stats']['diff'][clus_lab][mod_lab] = {'dur':{}, 'amp':{}}
            # duration data
            data['stats']['diff'][clus_lab][mod_lab]['dur']['result'] = stats.anderson_ksamp([ctrl_data['all'][clus_lab]['dur'],
                data['all'][clus_lab][mod_lab]['dur']])[-1]
            data['stats']['diff'][clus_lab][mod_lab]['dur']['label'] = 'Anderson-Darling; control:modulated duration'
            # amplitude data
            data['stats']['diff'][clus_lab][mod_lab]['amp']['result'] = stats.anderson_ksamp([ctrl_data['all'][clus_lab]['amp'],
                data['all'][clus_lab][mod_lab]['amp']])[-1]
            data['stats']['diff'][clus_lab][mod_lab]['amp']['label'] = 'Anderson-Darling; control:modulated amplitude'
    

    # tests whether duration and amplitude values for each modulation target are signficantly different to each other for
        # proximal and distal clustered stimulation
        
    mod_targets = ['on-site']
    mod_targets.extend(mod_labels)
    
    variables = ['dur', 'amp']
    
    for i, tar in enumerate(mod_targets):
        
        data['stats']['diff'][tar] = {'dur':{}, 'amp':{}}
        
        if tar == 'on-site':
            for var in variables:
                data['stats']['diff'][tar][var]['result'] = stats.anderson_ksamp([data['all'][clus_labels[0]][clus_labels[0]][var],
                     data['all'][clus_labels[1]][clus_labels[1]][var]])[-1]
                data['stats']['diff'][tar][var]['label'] = 'Anderson-Darling; proximal:distal duration'
            
        else:
            for var in variables:
                data['stats']['diff'][tar][var]['result'] = stats.anderson_ksamp([data['all'][clus_labels[0]][tar][var],
                     data['all'][clus_labels[1]][tar][var]])[-1]
                data['stats']['diff'][tar][var]['label'] = 'Anderson-Darling; proximal:distal duration'
        

    
    