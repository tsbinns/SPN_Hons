
# Minimal example on how to use the model libraries used in:
#   Lindroos & Hellgren Kotaleski 2020


from   neuron           import h
import numpy                as np
import MSN_builder          as build
import pickle
import matplotlib.pyplot    as plt
import sys
import common_functions     as cf

# Load model mechanisms
#if 'neuron' not in sys.modules and 'neuron' not in dir():
import neuron               as nrn
#nrn.load_mechanisms('mechanisms/single')


h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')

# specs
specs = {'dspn': {
                    'N': 71,
                    'lib': 'Libraries/D1_71bestFit_updRheob.pkl',
                    'par': 'params_dMSN.json',
                    'morph': 'Morphologies/WT-dMSN_P270-20_1.02_SGA1-m24.swc'},
         'ispn': {
                    'N': 34,
                    'lib': 'Libraries/D2_34bestFit_updRheob.pkl',
                    'par': 'params_iMSN.json',
                    'morph': 'Morphologies/WT-iMSN_P270-09_1.01_SGA2-m1.swc'}
        }
        




# chose cell type ('ispn' or 'dspn') and model id(s) to simulate
#---------------------------------------------------------------
cell_type         = 'ispn'    # 'dspn'/'ispn'
model_iterator    = [1]  # range(specs[cell_type]['N']) gives all models
# for dspn, 10 has lowest rheo, 54 has highest, 22 has mean, 41 has median; 22 is also average for experimental value
# for ispn, 8 has mean and median; 1 is average for experimental value

# stimulation details
if cell_type == 'dspn':
    target = ['dend[49]','dend[51]']
elif cell_type == 'ispn':
    target = ['dend[12]','dend[17]']
else:
    raise ValueError('The specified cell_type is unsupported')
target_labels = ['proximal','distal']


# open library (channel distributions etc)   
with open(specs[cell_type]['lib'], 'rb') as f:
    model_sets = pickle.load(f, encoding="latin1")

# simulate model(s)
OUT = {}
stim_n = 16
stim_t = 100
stop_t = stim_t+250 # when to stop simulating
isi = 1
pre_t = -50 # how much to time at the start to cut when plotting

for cell_n, cell_index in enumerate(model_iterator): 
    
    print('Simulating cell specification {} of {}'.format(cell_n+1,len(model_iterator)))
    
    OUT[cell_index] = {}
    tm = {}
    vm = {}
    for i, t in enumerate(target):
        # initiate cell
        cell = build.MSN(  params=specs[cell_type]['par'],
                           morphology=specs[cell_type]['morph'],
                           variables=model_sets[cell_index]['variables']   )
        
        rheobase        =   model_sets[cell_index]['rheobase']
        
        # record vectors
        tm[i]  = h.Vector()
        tm[i].record(h._ref_t)
        vm[i]  = h.Vector()
        vm[i].record(cell.soma(0.5)._ref_v)
        # add clustered inputs
        syn,stim,ncon,d2soma = cf.set_clustered_stim(   cell,t,n=stim_n, \
                                                        act_time=stim_t, \
                                                        ISI=isi          )
        # run simulation
        h.finitialize(-80)
        while h.t < stop_t:
            h.fadvance()
        tm[i] = tm[i].to_python()
        vm[i] = vm[i].to_python()
        # calculate dpp duration and amplitude
        base_t = next(i for i, x in enumerate(tm[i]) if x >= stim_t) - 1
        dur = cf.HMDur(tm[i],vm[i],base_t)
        amp = max(vm[i]) - vm[i][base_t]
        # store data
        OUT[cell_index][i] = {'tm':tm[i], 'vm':vm[i], 'dist':d2soma, \
                                  'rheo':rheobase, 'dur':dur, 'amp':amp}


print('Simulations completed, now calculating...')

# averaging ================

OUT_avg = {}

for i in range(len(target)):
    OUT_avg[i] = {'vm':[], 'avg_vm':[], 'rheo':[], 'avg_rheo':[], \
                      'dur':[], 'avg_dur':[], 'amp':[], 'avg_amp':[]}
    for cell_index in model_iterator:
        OUT_avg[i]['vm'].append(OUT[cell_index][i]['vm'])
        OUT_avg[i]['rheo'].append(OUT[cell_index][i]['rheo'])
        OUT_avg[i]['dur'].append(OUT[cell_index][i]['dur'])
        OUT_avg[i]['amp'].append(OUT[cell_index][i]['amp'])
    OUT_avg[i]['avg_vm'] = np.ndarray.tolist(np.mean(OUT_avg[i]['vm'],axis=0))
    OUT_avg[i]['avg_rheo'] = float(np.mean(OUT_avg[i]['rheo']))
    OUT_avg[i]['avg_dur'] = float(np.mean(OUT_avg[i]['dur']))
    OUT_avg[i]['avg_amp'] = float(np.mean(OUT_avg[i]['amp']))
OUT_avg['meta'] = {'cell_type':cell_type, 'tm': OUT[cell_index][0]['tm'], \
                   'dist': [OUT[cell_index][0]['dist'],OUT[cell_index][1]['dist']], \
                   'stim_n':stim_n, 'isi':isi, 'stim_t':stim_t, \
                   'stop_t':stop_t, 'pre_t':pre_t, 'labels': target_labels, \
                   'targets':target, 'specs':model_iterator}        

# plotting =================

stim_n = OUT_avg['meta']['stim_n']
stim_t = OUT_avg['meta']['stim_t']
stop_t = OUT_avg['meta']['stop_t']
pre_t = OUT_avg['meta']['pre_t']
isi = OUT_avg['meta']['isi']
cell_type = OUT_avg['meta']['cell_type']
targets = OUT_avg['meta']['targets']
target_labels = OUT_avg['meta']['labels']
model_iterator = OUT_avg['meta']['specs']

keys = list(OUT_avg.keys())
# plot Vm
if len(OUT_avg[keys[0]]['amp']) == 1:
    title = cell_type + ' ({}), n = '.format(model_iterator[0]) + \
            str(stim_n) + ', rheo = ' + str(int(OUT_avg[keys[0]]['avg_rheo']))
else:
    title = cell_type + ' (avg), n = ' + str(stim_n) + ', rheo = ' + \
            str(int(OUT_avg[keys[0]]['avg_rheo']))

plt.figure()
axs = plt.subplot(111)
for i in range(len(OUT_avg)-1):
    lab = OUT_avg['meta']['labels'][i] + \
              ' ({} $\mu$m)'.format(OUT_avg['meta']['dist'][i])
        
    plt.plot(OUT_avg['meta']['tm'],OUT_avg[keys[i]]['avg_vm'], label=lab)
    
plt.legend()
plt.show()

plt.xlim(stim_t+pre_t, stop_t)
plt.xticks(ticks=np.arange(stim_t+pre_t, stop_t+1, step=50), \
           labels=np.arange(pre_t, stop_t+pre_t+1, step=50))
plt.xlabel('time (ms)')
plt.ylabel('membrane potential (mV)')
plt.title(title)
axs.spines['right'].set_visible(False)
axs.spines['top'].set_visible(False)
# shade area of stimulation
shade_x = [stim_t,stim_t+stim_n*isi,stim_t+stim_n*isi,stim_t]
shade_y = [plt.ylim()[0],plt.ylim()[0],plt.ylim()[1],plt.ylim()[1]]
plt.fill(shade_x,shade_y,color='darkgrey',alpha=.2)

if len(model_iterator) > 1:
    # plot potential durations and amplitudes
    fig, axs = plt.subplots(1,2)
    fig.suptitle(title)
    # plots duration data
    axs[0].boxplot([OUT_avg['0']['dur'],OUT_avg['1']['dur']],widths=.6)
    axs[0].set_aspect(.08)
    axs[0].set_xticklabels(target_labels)
    axs[0].set_ylabel('duration (ms)')
    axs[0].spines['right'].set_visible(False)
    axs[0].spines['top'].set_visible(False)
    # plots amplitude data
    axs[1].boxplot([OUT_avg['0']['amp'],OUT_avg['1']['amp']],widths=.6)
    axs[1].set_aspect(.6)
    axs[1].set_xticklabels(target_labels)
    axs[1].set_ylabel('amplitude (mV)')
    axs[1].spines['right'].set_visible(False)
    axs[1].spines['top'].set_visible(False)



# save/load data ============
#cf.save_data(data,path)
#data = cf.load_data(path)
