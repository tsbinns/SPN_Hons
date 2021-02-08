'''
Plots voltage traces as well as potential duration and amplitude information
    for the cholinergic modulation of the plateau potential
'''


import numpy                as np
import matplotlib.pyplot    as plt
import common_functions     as cf
import scipy.stats          as stats
     


# load data
data = cf.load_data('Data/ispn_HFI[0]+0_modulation.json')#cf.load_data('C:/Users/tomth/OneDrive/Documents/Work/Courses/Level 4/Honours/Data/dspn_n16.json')
ctrl_data = cf.load_data('Data/ispn_HFI[0]+0_validation.json')


# ===== organise data =====

# clustered stimulation data
stim_n = data['meta']['clustered']['params']['stim_n']
pre_t = data['meta']['clustered']['params']['pre_t']
isi = data['meta']['clustered']['params']['isi']
clus_stim_t = data['meta']['clustered']['params']['stim_t']
clus_stop_t = data['meta']['clustered']['params']['stop_t']
clus_targets = data['meta']['clustered']['target']
clus_labels = data['meta']['clustered']['label']

# cholinergic stimulation data
ACh_stim_t = data['meta']['ACh info']['params']['stim_t']
ACh_stop_t = data['meta']['ACh info']['params']['stop_t']
ACh_targets = data['meta']['ACh info']['target']
ACh_labels = data['meta']['ACh info']['label']

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
    
    for j, ACh_lab in enumerate(ACh_labels):
        
        # plots voltage values for modulation applied to non-clustered input sites
        lab = ACh_lab
        axs[i].plot(data['meta']['tm'], data['avg'][clus_lab][ACh_lab]['vm'], label=lab, color=cols[col_i])
        # plots sem shading
        sem = stats.sem(data['avg'][clus_lab][ACh_lab]['vm'], axis=0)
        sem_plus = data['avg'][clus_lab][ACh_lab]['vm'] + sem
        sem_minus = data['avg'][clus_lab][ACh_lab]['vm'] - sem
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
    axs[i].plot([ACh_stim_t,ACh_stop_t],[axs[i].get_ylim()[0],axs[i].get_ylim()[0]], \
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
    for i, ACh_lab in enumerate(ACh_labels):
        dur_data[ACh_lab] = []
    for i, clus_lab in enumerate(clus_labels): # clustered input site duration data
        dur_data['on-site'].append(data['all'][clus_lab][clus_lab]['dur'])
        for i, ACh_lab in enumerate(ACh_labels): # off-site and soma duration data
            dur_data[ACh_lab].append(data['all'][clus_lab][ACh_lab]['dur'])
    
    # plots duration data
    bp_labels = {'groups':['proximal dend','distal dend'], 'axes':[]}
    bp_labels['axes'] = {'y':'duration (ms)'}
    cf.grouped_boxplot(dur_data, axs[0], bp_labels, cols[:len(dur_data)])
    axs[0].spines['right'].set_visible(False)
    axs[0].spines['top'].set_visible(False)
    
    # groups amplitude data for plotting
    amp_data['on-site'] = []
    for i, ACh_lab in enumerate(ACh_labels):
        amp_data[ACh_lab] = []
    for i, clus_lab in enumerate(clus_labels): # clustered input site duration data
        amp_data['on-site'].append(data['all'][clus_lab][clus_lab]['amp'])
        for i, ACh_lab in enumerate(ACh_labels): # off-site and soma duration data
            amp_data[ACh_lab].append(data['all'][clus_lab][ACh_lab]['amp'])
    
    # plots amplitude data
    bp_labels = {'groups':['proximal dend','distal dend'], 'axes':[]}
    bp_labels['axes'] = {'y':'amplitude (mV)'}
    cf.grouped_boxplot(amp_data, axs[1], bp_labels, cols[:len(amp_data)])
    axs[1].spines['right'].set_visible(False)
    axs[1].spines['top'].set_visible(False)

    plt.tight_layout()

