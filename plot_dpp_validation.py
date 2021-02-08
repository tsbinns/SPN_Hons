'''
Plots voltage traces as well as potential duration and amplitude information.
    - Currently implemented for plotting basic plateau potentials
'''


import numpy                as np
import matplotlib.pyplot    as plt
import scipy.stats          as stats
import common_functions     as cf
     

HFI = 0

if not HFI:
    
    # load data
    data = cf.load_data('Data/dspn_HFI[0]+0_validation.json')
    
    
    # plotting =================
    
    clus_info = data['meta']['clustered']
    
    # simulation data
    stim_n = clus_info['params']['stim_n']
    stim_t = clus_info['params']['stim_t']
    stop_t = clus_info['params']['stop_t']
    pre_t = clus_info['params']['pre_t']
    isi = clus_info['params']['isi']
    cell_type = data['meta']['cell type']    
    targets = clus_info['target']
    target_labels = clus_info['label']
    model_iterator = data['meta']['iterations']
    n_rounds = data['meta']['n rounds']
    
    labels = ['proximal dendrite', 'distal dendrite']
    
    # plot voltage traces =====
    
    plt.figure()
    axs = plt.subplot(111)
    for i, lab in enumerate(target_labels):
        
        plt.plot(data['meta']['tm'],data['avg'][lab]['vm'], label=labels[i])
        
        # calculate sem
        sem = stats.sem(data['all'][lab]['vm'],axis=0)
        sem_plus = data['avg'][lab]['vm'] + sem
        sem_minus = data['avg'][lab]['vm'] - sem
        # plot std dev shading
        axs.fill_between(data['meta']['tm'], sem_plus, sem_minus, alpha=.1)
        
    plt.legend()
    plt.show()
    
    # ignores data at start of simulation before voltage reaches baseline
    plt.xlim(stim_t+pre_t, stop_t)
    plt.xticks(ticks=np.arange(stim_t+pre_t, stop_t+1, step=50), \
               labels=np.arange(pre_t, stop_t+pre_t+1, step=50))
    
    plt.xlabel('time (ms)')
    plt.ylabel('membrane potential (mV)')
    plt.title(cell_type)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    
    # underscore area of stimulation
    axs.plot([stim_t,stim_t+stim_n*isi],[plt.ylim()[0],plt.ylim()[0]], \
             linewidth=5,color='red',solid_capstyle='butt')
    # shade area of stimulation
    #shade_x = [stim_t,stim_t+stim_n*isi,stim_t+stim_n*isi,stim_t]
    #shade_y = [plt.ylim()[0],plt.ylim()[0],plt.ylim()[1],plt.ylim()[1]]
    #plt.fill(shade_x,shade_y,color='darkgrey',alpha=.2)
    plt.tight_layout(True)
    
    # plot duration and amplitude data =====
    if len(model_iterator) > 1 or n_rounds > 1:
        
        fig, axs = plt.subplots(1,2)
        fig.suptitle(cell_type)
        plt.tight_layout(True)
        
        # plots duration data
        axs[0].boxplot([data['all'][target_labels[0]]['dur'],data['all'][target_labels[1]]['dur']],widths=.6)
        axs[0].set_xticklabels(labels)
        axs[0].set_ylabel('duration (ms)')
        axs[0].spines['right'].set_visible(False)
        axs[0].spines['top'].set_visible(False)
        plt.tight_layout(True)
        
        # plots amplitude data
        axs[1].boxplot([data['all'][target_labels[0]]['amp'],data['all'][target_labels[1]]['amp']],widths=.6)
        axs[1].set_xticklabels(labels)
        axs[1].set_ylabel('amplitude (mV)')
        axs[1].spines['right'].set_visible(False)
        axs[1].spines['top'].set_visible(False)
        plt.tight_layout(True)


else:
    
    colors = (plt.rcParams['axes.prop_cycle']).by_key()['color']
    delta = [0, 20, 40, 60, 80, 100]
    delta_labels = []
    spiking = {'avg':{}, 'sem':{}}
    
    for r in range(len(delta)):
        
        # load data
        data = cf.load_data('Data/dspn_HFI[1]+{}_validation.json'.format(delta[r]))
        '''
        new_spiked = {'proximal dend': [], 'distal dend': []}
        keys = list(data.keys())
        tm_stop = data['meta']['tm'].index(min(data['meta']['tm'], key=lambda x:abs(x-(HFI_info['stim_t']+75))))
        for i, key in enumerate(keys):
            if key != 'all' and key != 'meta':
                for i, lab in enumerate(target_labels):
                    if max(data[key]['0'][lab]['vm'][:tm_stop]) > 0:
                        new_spiked[lab].append(1)
                    else:
                        new_spiked[lab].append(0)
        
        for i, lab in enumerate(target_labels):
            data['all'][lab]['spiked'] = new_spiked[lab]
        '''
        delta_labels.append('+{}'.format(delta[r]))
        
        # plotting =================
        
        clus_info = data['meta']['clustered']
        HFI_info = data['meta']['HFI']
        
        # simulation data
        stim_n = clus_info['params']['stim_n']
        stim_t = clus_info['params']['stim_t']
        stop_t = clus_info['params']['stop_t']
        pre_t = clus_info['params']['pre_t']
        isi = clus_info['params']['isi']
        cell_type = data['meta']['cell type']
        targets = clus_info['target']
        target_labels = clus_info['label']
        model_iterator = data['meta']['iterations']
        n_rounds = data['meta']['n rounds']
        
        labels = ['proximal dendrite', 'distal dendrite']
        
        # plot voltage traces =====
        
        fig, axs = plt.subplots(2,1)
        fig.suptitle(cell_type + ' ({})'.format(delta_labels[r]))
        for i, lab in enumerate(target_labels):
            
            # avg firing probability at each time point
            if r == 0:
                spiking['avg'][lab] = []
            spiking['avg'][lab].append(np.mean(data['all'][lab]['spiked']))
            
            # std dev of firing probability at each time point
            if r == 0:
                spiking['sem'][lab] = []
            spiking['sem'][lab].append(stats.sem(data['all'][lab]['spiked']))
                
            for j in range(len(data['meta']['iterations'])*data['meta']['n rounds']):
                if data['all'][lab]['spiked'][j] == 1:
                    col = colors[i]
                else:
                    col = 'grey'
                axs[i].plot(data['meta']['tm'],data['all'][lab]['vm'][j], color=col)
         
        for i, lab in enumerate(target_labels):
            
            axs[i].set_xlabel('time (ms)')
            axs[i].set_ylabel('membrane potential (mV)')
            axs[i].set_title(labels[i])
            axs[i].spines['right'].set_visible(False)
            axs[i].spines['top'].set_visible(False)
            
            # underscore area of clustered stimulation
            axs[i].plot([stim_t,stim_t+stim_n*isi],[plt.ylim()[0],plt.ylim()[0]], \
                     linewidth=5,color='black',solid_capstyle='butt')
            
            # underscore area of HFI
            axs[i].plot([HFI_info['stim_t'],HFI_info['stop_t']],[plt.ylim()[0],plt.ylim()[0]], \
                     linewidth=5,color='grey',solid_capstyle='butt')
            
            # ignores data at start of simulation before voltage reaches baseline
            axs[i].set_xlim(stim_t+pre_t, stop_t+delta[r]+(stim_n*isi))
            axs[i].set_xticks(np.arange(stim_t+pre_t, stop_t+delta[r]+(stim_n*isi)+1, step=50))
            axs[i].set_xticklabels(np.arange(pre_t,stop_t-stim_t+(stim_n*isi)+1+delta[r],step=50))
            
            plt.tight_layout(True)
            
    
    # plot duration and amplitude data =====
    if len(model_iterator) > 1 or n_rounds > 1:
        
        plt.figure()
        axs = plt.subplot(111)
        
        # plots spike probability data
        for i, lab in enumerate(target_labels):
            
            axs.errorbar(delta, spiking['avg'][lab], yerr=spiking['sem'][lab], color=colors[i], label=labels[i])
            
        axs.set_xticklabels([0]+delta_labels)
        axs.set_xlabel('delta (ms)')
        axs.set_ylabel('spike probability')
        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)
        axs.set_ylim(0,1)
        
        plt.tight_layout(True)
        
        
        