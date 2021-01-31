'''
Plots voltage traces as well as potential duration and amplitude information
    for the cholinergic modulation of the plateau potential
'''


import numpy                as np
import matplotlib.pyplot    as plt
import common_functions     as cf
     


# load data
data_avg = cf.load_data('C:/Users/tomth/OneDrive/Documents/Work/Courses/Level 4/Honours/Data/dspn_n16.json')


# plotting =================


# simulation data
stim_n = data_avg['meta']['stim_n']
stim_t = data_avg['meta']['stim_t']
stop_t = data_avg['meta']['stop_t']
pre_t = data_avg['meta']['pre_t']
isi = data_avg['meta']['isi']
cell_type = data_avg['meta']['cell_type']
targets = data_avg['meta']['targets']
target_labels = data_avg['meta']['labels']
model_iterator = data_avg['meta']['specs']

keys = list(data_avg.keys())


# plot voltage traces =====
if len(data_avg[keys[0]]['amp']) == 1:
    title = cell_type + ' ({}), n = '.format(model_iterator[0]) + \
            str(stim_n) + ', rheo = ' + str(int(data_avg[keys[0]]['avg_rheo']))
else:
    title = cell_type + ' (avg), n = ' + str(stim_n) + ', rheo = ' + \
            str(int(data_avg[keys[0]]['avg_rheo']))
plt.figure()
axs = plt.subplot(111)
for i in range(len(targets)):
    lab = data_avg['meta']['labels'][i] + \
              ' ({} $\mu$m)'.format(data_avg['meta']['dist'][i])
        
    plt.plot(data_avg['meta']['tm'],data_avg[keys[i]]['avg_vm'], label=lab)
    
    # calculate std dev shading
    std = np.std(data_avg[keys[i]]['vm'],axis=0)
    std_plus = data_avg[keys[i]]['avg_vm'] + std
    std_minus = data_avg[keys[i]]['avg_vm'] - std
    # plot std dev shading
    axs.fill_between(data_avg['meta']['tm'],std_plus,std_minus,alpha=.1)
    
plt.legend()
plt.show()

# ignores data at start of simulation before voltage reaches baseline
plt.xlim(stim_t+pre_t, stop_t)
plt.xticks(ticks=np.arange(stim_t+pre_t, stop_t+1, step=50), \
           labels=np.arange(pre_t, stop_t+pre_t+1, step=50))

plt.xlabel('time (ms)')
plt.ylabel('membrane potential (mV)')
plt.title(title)
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
if len(model_iterator) > 1:
    
    fig, axs = plt.subplots(1,2)
    fig.suptitle(title)
    
    # plots duration data
    axs[0].boxplot([data_avg['0']['dur'],data_avg['1']['dur']],widths=.6)
    axs[0].set_aspect(.08)
    axs[0].set_xticklabels(target_labels)
    axs[0].set_ylabel('duration (ms)')
    axs[0].spines['right'].set_visible(False)
    axs[0].spines['top'].set_visible(False)
    
    # plots amplitude data
    axs[1].boxplot([data_avg['0']['amp'],data_avg['1']['amp']],widths=.6)
    axs[1].set_aspect(.6)
    axs[1].set_xticklabels(target_labels)
    axs[1].set_ylabel('amplitude (mV)')
    axs[1].spines['right'].set_visible(False)
    axs[1].spines['top'].set_visible(False)

