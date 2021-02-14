'''
Plots voltage traces as well as potential duration and amplitude information
    for the cholinergic modulation of the plateau potential
'''


import numpy                as np
import matplotlib.pyplot    as plt
import common_functions     as cf
import scipy.stats          as stats
     



HFI = 0
cell_type = 'dspn'
mod_type = 'DA'
mod_tar = 'all'

if HFI == 0:
    
    if mod_tar == 'indiv':
        
        # load data
        data = cf.load_data('Data/{}_HFI[0]+0_{}-modulation-{}.json'.format(cell_type, mod_type, mod_tar))#cf.load_data('C:/Users/tomth/OneDrive/Documents/Work/Courses/Level 4/Honours/Data/dspn_n16.json')
        ctrl_data = cf.load_data('Data/{}_HFI[0]+0_validation.json'.format(cell_type))
        
        
        # ===== organise data =====
        
        # clustered stimulation data
        stim_n = data['meta']['clustered']['params']['stim_n']
        pre_t = data['meta']['clustered']['params']['pre_t']
        isi = data['meta']['clustered']['params']['isi']
        clus_stim_t = data['meta']['clustered']['params']['stim_t']
        clus_stop_t = data['meta']['clustered']['params']['stop_t']
        clus_targets = data['meta']['clustered']['target']
        clus_labels = data['meta']['clustered']['label']
        
        # modulation data
        mod_stim_t = data['meta'][mod_type + ' info']['params']['stim_t']
        mod_stop_t = data['meta'][mod_type + ' info']['params']['stop_t']
        mod_targets = data['meta'][mod_type + ' info']['target']
        mod_labels = data['meta'][mod_type + ' info']['label']
        
        # simulation data
        cell_type = data['meta']['cell type']
        model_iterator = data['meta']['iterations']
        
        
        # gets control data for the cell iterations being modelled
        ctrl_avg = {}
        dur_data = {}
        dur_data['control'] = []
        amp_data = {}
        amp_data['control'] = []
        for i, clus_lab in enumerate(clus_labels):
            ctrl_avg[clus_lab] = {}
            # voltage values
            #ctrl_vm = [ctrl_data['avg'][clus_lab]['vm'][x] for x in model_iterator]
            #ctrl_avg[clus_lab]['avg_vm'] = np.ndarray.tolist(np.mean(ctrl_vm, axis=0))
            ctrl_avg[clus_lab]['avg_vm'] = ctrl_data['avg'][clus_lab]['vm']
            # duration
            dur_data['control'].append([ctrl_data['all'][clus_lab]['dur'][x] for x in model_iterator])
            # amplitude
            amp_data['control'].append([ctrl_data['all'][clus_lab]['amp'][x] for x in model_iterator])
        
        
        # ===== plot voltage traces =====
                
        fig_v, axs = plt.subplots(1,2)
        fig_v.suptitle(cell_type)
        
        # colours of lines
        cols = ['black']
        cols.extend((plt.rcParams['axes.prop_cycle']).by_key()['color'])
        
        for i, clus_lab in enumerate(clus_labels):
            
            axs[i].set_title(clus_lab)
            col_i = 0
            
            # plots control voltage traces (clustered input only, no modulation)
            lab = 'control'
            axs[i].plot(data['meta']['tm'], ctrl_avg[clus_lab]['avg_vm'], linewidth=5, color=cols[col_i], label=lab)
            col_i += 1
            
            # plots voltage values for modulation applied to same site as clustered input
            axs[i].plot(data['meta']['tm'], data['avg'][clus_lab][clus_lab]['vm'], label='on-site', color=cols[col_i])
            # plots sem shading
            sem = stats.sem(data['avg'][clus_lab][clus_lab]['vm'], axis=0)
            sem_plus = data['avg'][clus_lab][clus_lab]['vm'] + sem
            sem_minus = data['avg'][clus_lab][clus_lab]['vm'] - sem
            axs[i].fill_between(data['meta']['tm'], sem_plus, sem_minus, alpha=.1, color=cols[col_i])
            col_i += 1
            
            for j, mod_lab in enumerate(mod_labels):
                
                # plots voltage values for modulation applied to non-clustered input sites
                lab = mod_lab
                axs[i].plot(data['meta']['tm'], data['avg'][clus_lab][mod_lab]['vm'], label=lab, color=cols[col_i])
                # plots sem shading
                sem = stats.sem(data['avg'][clus_lab][mod_lab]['vm'], axis=0)
                sem_plus = data['avg'][clus_lab][mod_lab]['vm'] + sem
                sem_minus = data['avg'][clus_lab][mod_lab]['vm'] - sem
                axs[i].fill_between(data['meta']['tm'], sem_plus, sem_minus, alpha=.1, color=cols[col_i])
                col_i += 1
            
            axs[i].legend()
            
            # ignores data at start of simulation before voltage reaches baseline
            axs[i].set_xlim(clus_stim_t+pre_t, clus_stop_t)
            axs[i].set_xticks(np.arange(clus_stim_t+pre_t, clus_stop_t+1, step=50))
            axs[i].set_xticklabels(np.arange(pre_t, clus_stop_t+pre_t+1, step=50))
            
            axs[i].set_xlabel('time (ms)')
            axs[i].set_ylabel('membrane potential (mV)')
            axs[i].spines['right'].set_visible(False)
            axs[i].spines['top'].set_visible(False)
            
            # underscore time of clustered stimulation
            axs[i].plot([clus_stim_t,clus_stim_t+stim_n*isi],[axs[i].get_ylim()[0],axs[i].get_ylim()[0]], \
                linewidth=5,color='red',solid_capstyle='butt')
            
            # underscore time of cholinergic modulation
            axs[i].plot([mod_stim_t,mod_stop_t],[axs[i].get_ylim()[0],axs[i].get_ylim()[0]], \
                linewidth=3,color='cyan',linestyle='--',solid_capstyle='butt')
        
        
        # sets axes of subplots to equal one another
        axs_ylims = []
        axs_ylims = np.vstack([axs[i].get_ylim() for i in range(len(axs))])
        new_ylims = np.append(min(axs_ylims[:,0]),max(axs_ylims[:,1]))
        for i in range(len(axs)):
            axs[i].set_ylim(new_ylims)
            
        plt.show()
        
        
        
        # ===== plot duration and amplitude data =====
        
        if len(model_iterator) > 1:
               
            fig, axs = plt.subplots(1,2)
            fig.suptitle(cell_type)
            
            # groups duration data for plotting
            dur_data['on-site'] = []
            for i, mod_lab in enumerate(mod_labels):
                dur_data[mod_lab] = []
            for i, clus_lab in enumerate(clus_labels): # clustered input site duration data
                dur_data['on-site'].append(data['all'][clus_lab][clus_lab]['dur'])
                for i, mod_lab in enumerate(mod_labels): # off-site and soma duration data
                    dur_data[mod_lab].append(data['all'][clus_lab][mod_lab]['dur'])
            
            # plots duration data
            bp_labels = {'groups':['proximal dend','distal dend'], 'axes':[]}
            bp_labels['axes'] = {'y':'duration (ms)'}
            cf.grouped_boxplot(dur_data, axs[0], bp_labels, cols[:len(dur_data)])
            axs[0].spines['right'].set_visible(False)
            axs[0].spines['top'].set_visible(False)
            
            # groups amplitude data for plotting
            amp_data['on-site'] = []
            for i, mod_lab in enumerate(mod_labels):
                amp_data[mod_lab] = []
            for i, clus_lab in enumerate(clus_labels): # clustered input site duration data
                amp_data['on-site'].append(data['all'][clus_lab][clus_lab]['amp'])
                for i, mod_lab in enumerate(mod_labels): # off-site and soma duration data
                    amp_data[mod_lab].append(data['all'][clus_lab][mod_lab]['amp'])
            
            # plots amplitude data
            bp_labels = {'groups':['proximal dend','distal dend'], 'axes':[]}
            bp_labels['axes'] = {'y':'amplitude (mV)'}
            cf.grouped_boxplot(amp_data, axs[1], bp_labels, cols[:len(amp_data)])
            axs[1].spines['right'].set_visible(False)
            axs[1].spines['top'].set_visible(False)
        
            plt.tight_layout()
            
    else:
        # load data
        data = cf.load_data('Data/{}_HFI[0]+0_{}-modulation-{}.json'.format(cell_type, mod_type, mod_tar))#cf.load_data('C:/Users/tomth/OneDrive/Documents/Work/Courses/Level 4/Honours/Data/dspn_n16.json')
        ctrl_data = cf.load_data('Data/{}_HFI[0]+0_validation.json'.format(cell_type))
        
        
        # ===== organise data =====
        
        # clustered stimulation data
        stim_n = data['meta']['clustered']['params']['stim_n']
        pre_t = data['meta']['clustered']['params']['pre_t']
        isi = data['meta']['clustered']['params']['isi']
        clus_stim_t = data['meta']['clustered']['params']['stim_t']
        clus_stop_t = data['meta']['clustered']['params']['stop_t']
        clus_targets = data['meta']['clustered']['target']
        clus_labels = data['meta']['clustered']['label']
        
        # modulation data
        mod_stim_t = data['meta'][mod_type + ' info']['params']['stim_t']
        mod_stop_t = data['meta'][mod_type + ' info']['params']['stop_t']
        mod_targets = data['meta'][mod_type + ' info']['target']
        mod_labels = data['meta'][mod_type + ' info']['label']
        
        # simulation data
        cell_type = data['meta']['cell type']
        model_iterator = data['meta']['iterations']
        
        
        # gets control data for the cell iterations being modelled
        ctrl_avg = {}
        dur_data = {}
        dur_data['control'] = []
        amp_data = {}
        amp_data['control'] = []
        for i, clus_lab in enumerate(clus_labels):
            ctrl_avg[clus_lab] = {}
            # voltage values
            #ctrl_vm = [ctrl_data['avg'][clus_lab]['vm'][x] for x in model_iterator]
            #ctrl_avg[clus_lab]['avg_vm'] = np.ndarray.tolist(np.mean(ctrl_vm, axis=0))
            ctrl_avg[clus_lab]['avg_vm'] = ctrl_data['avg'][clus_lab]['vm']
            # duration
            dur_data['control'].append([ctrl_data['all'][clus_lab]['dur'][x] for x in model_iterator])
            # amplitude
            amp_data['control'].append([ctrl_data['all'][clus_lab]['amp'][x] for x in model_iterator])
        
        
        # ===== plot voltage traces =====
                
        fig_v, axs = plt.subplots(1,2)
        fig_v.suptitle(cell_type)
        
        # colours of lines
        cols = ['black']
        cols.extend((plt.rcParams['axes.prop_cycle']).by_key()['color'])
        
        for i, clus_lab in enumerate(clus_labels):
            
            axs[i].set_title(clus_lab)
            col_i = 0
            
            # plots control voltage traces (clustered input only, no modulation)
            lab = 'control'
            axs[i].plot(data['meta']['tm'], ctrl_avg[clus_lab]['avg_vm'], linewidth=5, color=cols[col_i], label=lab)
            col_i += 1
            
            for j, mod_lab in enumerate(mod_labels):
                
                # plots voltage values for modulation applied to non-clustered input sites
                lab = mod_lab
                axs[i].plot(data['meta']['tm'], data['avg'][clus_lab][mod_lab]['vm'], label=lab, color=cols[col_i])
                # plots sem shading
                sem = stats.sem(data['avg'][clus_lab][mod_lab]['vm'], axis=0)
                sem_plus = data['avg'][clus_lab][mod_lab]['vm'] + sem
                sem_minus = data['avg'][clus_lab][mod_lab]['vm'] - sem
                axs[i].fill_between(data['meta']['tm'], sem_plus, sem_minus, alpha=.1, color=cols[col_i])
                col_i += 1
            
            axs[i].legend()
            
            # ignores data at start of simulation before voltage reaches baseline
            axs[i].set_xlim(clus_stim_t+pre_t, clus_stop_t)
            axs[i].set_xticks(np.arange(clus_stim_t+pre_t, clus_stop_t+1, step=50))
            axs[i].set_xticklabels(np.arange(pre_t, clus_stop_t+pre_t+1, step=50))
            
            axs[i].set_xlabel('time (ms)')
            axs[i].set_ylabel('membrane potential (mV)')
            axs[i].spines['right'].set_visible(False)
            axs[i].spines['top'].set_visible(False)
            
            # underscore time of clustered stimulation
            axs[i].plot([clus_stim_t,clus_stim_t+stim_n*isi],[axs[i].get_ylim()[0],axs[i].get_ylim()[0]], \
                linewidth=5,color='red',solid_capstyle='butt')
            
            # underscore time of cholinergic modulation
            axs[i].plot([mod_stim_t,mod_stop_t],[axs[i].get_ylim()[0],axs[i].get_ylim()[0]], \
                linewidth=3,color='cyan',linestyle='--',solid_capstyle='butt')
        
        
        # sets axes of subplots to equal one another
        axs_ylims = []
        axs_ylims = np.vstack([axs[i].get_ylim() for i in range(len(axs))])
        new_ylims = np.append(min(axs_ylims[:,0]),max(axs_ylims[:,1]))
        for i in range(len(axs)):
            axs[i].set_ylim(new_ylims)
            
        plt.show()
        
        
        
        # ===== plot duration and amplitude data =====
        
        if len(model_iterator) > 1:
               
            fig, axs = plt.subplots(1,2)
            fig.suptitle(cell_type)
            
            # groups duration data for plotting
            for i, mod_lab in enumerate(mod_labels):
                dur_data[mod_lab] = []
            for i, clus_lab in enumerate(clus_labels):
                for i, mod_lab in enumerate(mod_labels): # off-site and soma duration data
                    dur_data[mod_lab].append(data['all'][clus_lab][mod_lab]['dur'])
            
            # plots duration data
            bp_labels = {'groups':['proximal dend','distal dend'], 'axes':[]}
            bp_labels['axes'] = {'y':'duration (ms)'}
            cf.grouped_boxplot(dur_data, axs[0], bp_labels, cols[:len(dur_data)])
            axs[0].spines['right'].set_visible(False)
            axs[0].spines['top'].set_visible(False)
            
            # groups amplitude data for plotting
            for i, mod_lab in enumerate(mod_labels):
                amp_data[mod_lab] = []
            for i, clus_lab in enumerate(clus_labels):
                for i, mod_lab in enumerate(mod_labels): # off-site and soma duration data
                    amp_data[mod_lab].append(data['all'][clus_lab][mod_lab]['amp'])
            
            # plots amplitude data
            bp_labels = {'groups':['proximal dend','distal dend'], 'axes':[]}
            bp_labels['axes'] = {'y':'amplitude (mV)'}
            cf.grouped_boxplot(amp_data, axs[1], bp_labels, cols[:len(amp_data)])
            axs[1].spines['right'].set_visible(False)
            axs[1].spines['top'].set_visible(False)
        
            plt.tight_layout()


else: # ===================

    
    colors = (plt.rcParams['axes.prop_cycle']).by_key()['color']
    delta = list(np.arange(0,100+1,20))
    delta_labels = []
    spiking = {}
    ctrl_spiking = {}
    
    labels = ['proximal dendrite', 'distal dendrite']
    ctrl_labels = ['proximal (control)', 'distal (control)']
    
    for d in range(len(delta)):
        
        # load data
        data = cf.load_data('Data/{}_HFI[1]+{}_{}-modulation.json'.format(cell_type, delta[d], mod_type))
        delta_labels.append('+{}'.format(delta[d]))
        
        ctrl_data = cf.load_data('Data/{}_HFI[1]+{}_validation.json'.format(cell_type, delta[d]))
        
        
        # plotting =================
        
        clus_info = data['meta']['clustered']
        HFI_info = data['meta']['HFI']
        mod_info = data['meta'][mod_type + ' info']
        
        # simulation data
        stim_n = clus_info['params']['stim_n']
        stim_t = clus_info['params']['stim_t']
        stop_t = clus_info['params']['stop_t']
        pre_t = clus_info['params']['pre_t']
        isi = clus_info['params']['isi']
        cell_type = data['meta']['cell type']
        targets = clus_info['target']
        clus_labels = clus_info['label']
        model_iterator = data['meta']['iterations']
        n_rounds = data['meta']['n rounds']
        
        mod_labels = mod_info['label']
        
        mod_targets = ['on-site']
        mod_targets.extend(mod_labels)
        
        labels = ['proximal dendrite', 'distal dendrite']
        
        
        # collects control data =====
        for i, clus_lab in enumerate(clus_labels):
        
            if d == 0:
                ctrl_spiking[clus_lab] = {'spiked':{'avg':[], 'sem':[]}}
            
            # avg spiking data at each time point and sem
            ctrl_spiking[clus_lab]['spiked']['avg'].append(np.mean(ctrl_data['all'][clus_lab]['spiked']))
            ctrl_spiking[clus_lab]['spiked']['sem'].append(stats.sem(ctrl_data['all'][clus_lab]['spiked']))
            
        
        # plot voltage traces =====
        
        for j, mod_lab in enumerate(mod_targets):

            use_clus = 0
            if mod_lab == 'on-site':
                use_clus = 1
            '''
            fig, axs = plt.subplots(2,1)
            fig.suptitle(cell_type + ', {} ({})'.format(ACh_lab,delta_labels[d]))
            '''
            for i, clus_lab in enumerate(clus_labels):

                if use_clus:
                    mod_lab = clus_lab
                
                if d == 0 and j == 0:
                    spiking[clus_lab] = {clus_lab:{}}
                if d == 0:
                    spiking[clus_lab][mod_lab] = {'spiked':{'avg':[], 'sem':[]}, 'first_spike':{'avg':[], 'sem':[]}, \
                           'spike_n':{'avg':[], 'sem':[]}}
                
                # avg spiking data at each time point and sem
                spiking[clus_lab][mod_lab]['spiked']['avg'].append(np.mean(data['all'][clus_lab][mod_lab]['spiked_avg']))
                spiking[clus_lab][mod_lab]['spiked']['sem'].append(stats.sem(data['all'][clus_lab][mod_lab]['spiked_avg']))
                
                '''
                for k in range(len(model_iterator)*(n_rounds)):
                    if data['all'][clus_lab][ACh_lab]['spiked'][k] == 1:
                        col = colors[i]
                    else:
                        col = 'grey'
                    axs[i].plot(data['meta']['tm'], data['all'][clus_lab][ACh_lab]['vm'][k], color=col)
            
            ylim = plt.ylim()
            ylim = [ylim[0]] * 2
                
            for i, lab in enumerate(clus_labels):
                
                axs[i].set_xlabel('time (ms)')
                axs[i].set_ylabel('membrane potential (mV)')
                axs[i].set_title(labels[i])
                axs[i].spines['right'].set_visible(False)
                axs[i].spines['top'].set_visible(False)
                
                # underscore area of clustered stimulation
                axs[i].plot([stim_t,stim_t+stim_n*isi], ylim, linewidth=5,color='red',solid_capstyle='butt')
                
                # underscore area of HFI
                axs[i].plot([HFI_info['stim_t'],HFI_info['stop_t']], ylim, linewidth=5,color='grey',solid_capstyle='butt')
                
                # ignores data at start of simulation before voltage reaches baseline
                axs[i].set_xlim(stim_t+pre_t, stop_t+delta[d]+(stim_n*isi))
                axs[i].set_xticks(np.arange(stim_t+pre_t, stop_t+delta[d]+(stim_n*isi)+1, step=50))
                axs[i].set_xticklabels(np.arange(pre_t,stop_t-stim_t+(stim_n*isi)+1+delta[d],step=50))
                
                plt.tight_layout()
                '''
    
    # plot spiking data =====
    #if len(model_iterator) > 1 or n_rounds > 1:
        
    # baseline firing rate for HFI without clustered inputs
    baseline = cf.load_data('Data/{}_HFI[1]+0_baseline.json'.format(cell_type))
    baseline_spiking = baseline['all']['proximal dend']['spiked']
    baseline_spiking.extend(baseline['all']['distal dend']['spiked'])
    baseline_spiking = np.mean(baseline_spiking)
    
    '''
    for j, ACh_lab in enumerate(ACh_targets):
    
        use_clus = 0
        if ACh_lab == 'on-site':
            use_clus = 1
        
        
        # plots spike probability
        plt.figure()
        axs = plt.subplot(111)
        axs.set_title(cell_type + ', {}'.format(ACh_targets[j]))
        
        # modulation data
        for i, clus_lab in enumerate(clus_labels):
            if use_clus:
                ACh_lab = clus_lab
            axs.errorbar(delta, spiking[clus_lab][ACh_lab]['spiked']['avg'], \
                         yerr=spiking[clus_lab][ACh_lab]['spiked']['sem'], color=colors[i], 
                         label=labels[i], capsize=5)
        
        # control data
        for i, clus_lab in enumerate(clus_labels):
            axs.plot(delta, ctrl_spiking[clus_lab]['spiked']['avg'], color=colors[i], alpha=.5, label=ctrl_labels[i])
        
        # baseline data
        axs.plot(axs.get_xlim(), [baseline_spiking]*2, linestyle='--', color='grey')
        
        axs.set_xticklabels([0]+delta_labels)
        axs.set_xlabel('delta (ms)')
        axs.set_ylabel('spike probability')
        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)
        axs.set_ylim(0,1)
        axs.legend()
        
        plt.tight_layout()
    '''
    
    fig, axs = plt.subplots(2,2)
    fig.suptitle(cell_type)
    
    for j, mod_lab in enumerate(mod_targets):
    
        use_clus = 0
        if mod_lab == 'on-site':
            use_clus = 1
        
        
        # plots spike probability
        pos1 = j//2
        pos2 = j%2
        axs[pos1,pos2].set_title(mod_lab)
        
        # modulation data
        for i, clus_lab in enumerate(clus_labels):
            if use_clus:
                mod_lab = clus_lab
            axs[pos1,pos2].errorbar(delta, spiking[clus_lab][mod_lab]['spiked']['avg'], \
                         yerr=spiking[clus_lab][mod_lab]['spiked']['sem'], color=colors[i], 
                         label=labels[i], capsize=5)
        
        # control data
        for i, clus_lab in enumerate(clus_labels):
            '''
            axs.errorbar(delta, ctrl_spiking[clus_lab]['spiked']['avg'], \
                         yerr=ctrl_spiking[clus_lab]['spiked']['sem'], color=colors[i], alpha=.5,
                         label=ctrl_labels[i], capsize=5)
            '''
            axs[pos1,pos2].plot(delta, ctrl_spiking[clus_lab]['spiked']['avg'], color=colors[i], alpha=.5, label=ctrl_labels[i])
        
        # baseline data
        axs[pos1,pos2].plot(axs[pos1,pos2].get_xlim(), [baseline_spiking]*2, linestyle='--', color='grey')
        
        axs[pos1,pos2].set_xticklabels([0]+delta_labels)
        axs[pos1,pos2].set_xlabel('delta (ms)')
        axs[pos1,pos2].set_ylabel('spike probability')
        axs[pos1,pos2].spines['right'].set_visible(False)
        axs[pos1,pos2].spines['top'].set_visible(False)
        axs[pos1,pos2].set_ylim(0,1)
        axs[pos1,pos2].legend()
        
    plt.tight_layout()
