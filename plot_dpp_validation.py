'''
Plots voltage traces as well as potential duration and amplitude information.
    - Currently implemented for plotting basic plateau potentials
'''


import numpy                as np
import matplotlib.pyplot    as plt
import common_functions     as cf
     


HFI = 1

if not HFI:
    
    # load data
    data = cf.load_data('Data/ispn_HFI[0]+0_validation.json')
    
    
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
    
    # plot voltage traces =====
    
    plt.figure()
    plt.tight_layout(True)
    axs = plt.subplot(111)
    for i, lab in enumerate(target_labels):
            
        plt.plot(data['meta']['tm'],data['avg'][lab]['vm'], label=lab)
        
        # calculate std dev shading
        std = np.std(data['all'][lab]['vm'],axis=0)
        std_plus = data['avg'][lab]['vm'] + std
        std_minus = data['avg'][lab]['vm'] - std
        # plot std dev shading
        axs.fill_between(data['meta']['tm'], std_plus, std_minus, alpha=.1)
        
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
             linewidth=5,color='black',solid_capstyle='butt')
    # shade area of stimulation
    #shade_x = [stim_t,stim_t+stim_n*isi,stim_t+stim_n*isi,stim_t]
    #shade_y = [plt.ylim()[0],plt.ylim()[0],plt.ylim()[1],plt.ylim()[1]]
    #plt.fill(shade_x,shade_y,color='darkgrey',alpha=.2)
    
    
    # plot duration and amplitude data =====
    if len(model_iterator) > 1 or n_rounds > 1:
        
        fig, axs = plt.subplots(1,2)
        fig.suptitle(cell_type)
        plt.tight_layout(True)
        
        # plots duration data
        axs[0].boxplot([data['all'][target_labels[0]]['dur'],data['all'][target_labels[1]]['dur']],widths=.6)
        axs[0].set_xticklabels(target_labels)
        axs[0].set_ylabel('duration (ms)')
        axs[0].spines['right'].set_visible(False)
        axs[0].spines['top'].set_visible(False)
        plt.tight_layout(True)
        
        # plots amplitude data
        axs[1].boxplot([data['all'][target_labels[0]]['amp'],data['all'][target_labels[1]]['amp']],widths=.6)
        axs[1].set_xticklabels(target_labels)
        axs[1].set_ylabel('amplitude (mV)')
        axs[1].spines['right'].set_visible(False)
        axs[1].spines['top'].set_visible(False)
        plt.tight_layout(True)


else:
    
    # load data
    data = cf.load_data('Data/dspn_HFI[1]+0_validation.json')
    
    
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
    
    
    # plot voltage traces =====
    
    colors = (plt.rcParams['axes.prop_cycle']).by_key()['color']
    
    plt.figure()
    plt.tight_layout(True)
    axs = plt.subplot(111)
    for i, lab in enumerate(target_labels):
            
        for j in range(len(data['meta']['iterations'])*data['meta']['n rounds']):
            
            if j == 1:
                plt.plot(data['meta']['tm'],data['all'][lab]['vm'][j], color=colors[i], label=lab)
            else:
                plt.plot(data['meta']['tm'],data['all'][lab]['vm'][j], color=colors[i])
        
    plt.legend()
    plt.show()
    
    plt.xlabel('time (ms)')
    plt.ylabel('membrane potential (mV)')
    plt.title(cell_type)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    
    # underscore area of clustered stimulation
    axs.plot([stim_t,stim_t+stim_n*isi],[plt.ylim()[0],plt.ylim()[0]], \
             linewidth=5,color='black',solid_capstyle='butt')
    
    # underscore area of HFI
    axs.plot([HFI_info['stim_t'],HFI_info['stop_t']],[plt.ylim()[0],plt.ylim()[0]], \
             linewidth=5,color='grey',solid_capstyle='butt')
    
    # ignores data at start of simulation before voltage reaches baseline
    plt.xlim(stim_t+pre_t, stop_t)
    plt.xticks(ticks=np.arange(stim_t+pre_t, stop_t+1, step=50), \
               labels=np.arange(pre_t, stop_t+pre_t+1, step=50))
    
    # plot duration and amplitude data =====
    if len(model_iterator) > 1 or n_rounds > 1:
        
        plt.figure()
        plt.tight_layout(True)
        axs = plt.subplot(111)
        
        # plots spike probability data
        axs.boxplot([data['all'][target_labels[0]]['spiked'],data['all'][target_labels[1]]['spiked']],widths=.6)
        axs.set_xticklabels(target_labels)
        axs.set_ylabel('spike probability')
        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)
        axs.set_ylim(0,1)
        plt.tight_layout(True)
        
        
        