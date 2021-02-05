'''
Provides clustered input to a SPN's dendrites (a distal dendrite and proximal dendrite, separately)
    to generate a pleateu potential (or not, in the case of the proximal dendritic stimulation),
    alongside background noise
'''

from   neuron               import h
import numpy                    as np
import pickle
import common_functions         as cf
import simulation_functions     as sf
import time



# ===== for parallelisation =====
h.nrnmpi_init()
pc = h.ParallelContext()



# ===== load model mechanisms/parameters =====
import neuron               as nrn
#nrn.load_mechanisms('mechanisms/single')

h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')

# specs
specs = {'dspn': {
                    'N': 71,
                    'lib': 'Libraries/D1_71bestFit_updRheob.pkl',
                    'par': 'Params/params_dMSN.json',
                    'morph': 'Morphologies/WT-dMSN_P270-20_1.02_SGA1-m24.swc'},
         'ispn': {
                    'N': 34,
                    'lib': 'Libraries/D2_34bestFit_updRheob.pkl',
                    'par': 'Params/params_iMSN.json',
                    'morph': 'Morphologies/WT-iMSN_P270-09_1.01_SGA2-m1.swc'}
        }
        
# chose cell type ('ispn' or 'dspn') and model id(s) to simulate...

cell_type = 'dspn'
if cell_type != 'dspn' and cell_type != 'ispn':
    raise ValueError("The requested cell type is not supported.\nOnly 'dpsn' and 'ispn' are recognised.")
    
#model_iterator = cf.iter_params(cell_type, only_ids=True)
model_iterator = [0]

iterations = model_iterator.copy()

n_rounds = 1
model_round = []
for r in range(n_rounds):
    for i in range(len(iterations)):
        model_round.append(r)
    if r > 0:
        model_iterator.extend(iterations)

if pc.id() == 0:
    print('Simulating {} cell iteration(s) of type: {}'.format(len(model_iterator),cell_type), \
          flush=True)
   
# open library (channel distributions etc)   
with open(specs[cell_type]['lib'], 'rb') as f:
    model_sets = pickle.load(f, encoding="latin1")



# ===== simulate model(s) =====
# model information to pass to simulations
model_data = {'specs':specs[cell_type], 'cell_type':cell_type, 'model_sets':model_sets}
noise = 1
HFI = 0
HFI_delay = 0
dur_and_amp = 1
spike = 0

start = time.time() # for timing simulations

pc.runworker() # start workers for parallelisation

# clear temp data folder
folder = 'temp_data'
cf.clear_folder(folder,'.json')

if pc.nhost() == 1: # use the serial form
    
    for cell_n in range(len(model_iterator)): # for each model
        cell_index = model_iterator[cell_n]
        # simulate model
        run_info = {'curr_n':cell_n, 'tot_n':len(model_iterator), 'round':model_round[cell_n]}
        data = sf.dpp_ACh_modded(model_data, cell_index, run_info, noise, HFI, HFI_delay, dur_and_amp, spike)
        # save file to folder
        name = '{}_{}-{}_modulation'.format(cell_type,data['meta']['round'], data['meta']['id'])
        cf.save_data(data,'{}/{}.json'.format(folder, name))
            
else: # use the bulleting board form
    
    for cell_n in range(len(model_iterator)): # scatter processes
        cell_index = model_iterator[cell_n]
        # simulate model
        run_info = {'curr_n':cell_n, 'tot_n':len(model_iterator), 'round':model_round[cell_n]}
        pc.submit(sf.dpp_ACh_modded, model_data, cell_index, run_info, noise, HFI, HFI_delay, dur_and_amp, spike)
        
    while pc.working(): # gather results
        data = pc.pyret()
        # save file to folder
        name = '{}_{}-{}_modulation'.format(cell_type,data['meta']['round'], data['meta']['id'])
        cf.save_data(data,'{}/{}.json'.format(folder, name))
  
  
pc.done() # end parallelisation


# for timing simulations
print('Simulations completed (took %.0f secs).\nNow performing calculations/collating data...' % (time.time()-start))

 

# ===== combine and/or average data =====

info = cf.params_for_input(cell_type, 'clustered')
clus_info = info['clustered']
info = cf.params_for_input(cell_type, 'ACh')
ACh_info = info['ACh']

# collates data loaded from files
data_all = {}

for i, iteration in enumerate(iterations):
    round_data = {}
    for r in range(n_rounds):
        # load data
        name = name = '{}_{}-{}_modulation'.format(cell_type,r,iteration)
        data = cf.load_data('{}/{}.json'.format(folder,name))
        round_data[r] = data
    # combine data
    data_all[i] = round_data
    
data = data_all
    
            
# collates data for each cell iteration
for i in range(len(iterations)):
    
    data[i]['all'] = {}
    
    for lab in clus_info['label']:
        
        data[i]['all'][lab] = {}
        data[i]['all'][lab][lab] = {'vm':[], 'dur':[], 'amp':[], 'spiked':[]}
        
        for j in range(n_rounds):
            
            data[i]['all'][lab][lab]['vm'].append(data[i][j][lab][lab]['vm'])
            if dur_and_amp:
                data[i]['all'][lab][lab]['dur'].append(data[i][j][lab][lab]['dur'])
                data[i]['all'][lab][lab]['amp'].append(data[i][j][lab][lab]['amp'])  
            if spike:
                data[i]['all'][lab][lab]['spiked'].append(data[i][j][lab][lab]['spiked'])
                
        for mod in ACh_info['label']:
            
            data[i]['all'][lab][mod] = {'vm':[], 'dur':[], 'amp':[], 'spiked':[]}
        
            for j in range(n_rounds):
                
                data[i]['all'][lab][mod]['vm'].append(data[i][j][lab][mod]['vm'])
                if dur_and_amp:
                    data[i]['all'][lab][mod]['dur'].append(data[i][j][lab][mod]['dur'])
                    data[i]['all'][lab][mod]['amp'].append(data[i][j][lab][mod]['amp'])  
                if spike:
                    data[i]['all'][lab][mod]['spiked'].append(data[i][j][lab][mod]['spiked'])
            
            
        data[i]['all']['meta'] = data[i][j]['meta']
        

# collates data across cell iterations
data['all'] = {}

for lab in clus_info['label']:
    
    data['all'][lab] = {}
    data['all'][lab][lab] = {'vm':[], 'dur':[], 'amp':[], 'spiked':[]}    

    for i in range(len(iterations)):
        
        data['all'][lab][lab]['vm'].extend(data[i]['all'][lab][lab]['vm'])
        if dur_and_amp:
            data['all'][lab][lab]['dur'].extend(data[i]['all'][lab][lab]['dur'])
            data['all'][lab][lab]['amp'].extend(data[i]['all'][lab][lab]['amp'])
        if spike:
            data['all'][lab][lab]['spiked'].extend(data[i]['all'][lab][lab]['spiked'])
            
    for mod in ACh_info['label']:
        
        data['all'][lab][mod] = {'vm':[], 'dur':[], 'amp':[], 'spiked':[]}    
    
        for i in range(len(iterations)):
            
            data['all'][lab][mod]['vm'].extend(data[i]['all'][lab][mod]['vm'])
            if dur_and_amp:
                data['all'][lab][mod]['dur'].extend(data[i]['all'][lab][mod]['dur'])
                data['all'][lab][mod]['amp'].extend(data[i]['all'][lab][mod]['amp'])
            if spike:
                data['all'][lab][mod]['spiked'].extend(data[i]['all'][lab][mod]['spiked'])
            

# collates meta data
data['meta'] = {'tm':data[i]['all']['meta']['tm'], 'cell type':cell_type, 'iterations':iterations,
                       'n rounds':n_rounds, 'clustered':clus_info, 'ACh info':ACh_info}
if noise:
    info = cf.params_for_input(cell_type, 'noise')
    data['meta']['noise'] = info['noise']
if HFI:
    info = cf.params_for_input(cell_type, 'HFI')
    info['HFI']['stim_t'] = clus_info['params']['stim_t'] + HFI_delay
    info['HFI']['stop_t'] = clus_info['params']['stop_t'] + HFI_delay
    data['meta']['HFI'] = info['HFI']
    

# grand averages data across cell iterations
if not HFI:
    data['avg'] = {}
    for lab in clus_info['label']:
        data['avg'][lab] = {}
        data['avg'][lab][lab] = {}
        data['avg'][lab][lab]['vm'] = np.ndarray.tolist(np.mean(data['all'][lab][lab]['vm'],axis=0))
        data['avg'][lab][lab]['dur'] = float(np.mean(data['all'][lab][lab]['dur']))
        data['avg'][lab][lab]['amp'] = float(np.mean(data['all'][lab][lab]['amp']))
        for mod in ACh_info['label']:
            data['avg'][lab][mod] = {}
            data['avg'][lab][mod]['vm'] = np.ndarray.tolist(np.mean(data['all'][lab][mod]['vm'],axis=0))
            data['avg'][lab][mod]['dur'] = float(np.mean(data['all'][lab][mod]['dur']))
            data['avg'][lab][mod]['amp'] = float(np.mean(data['all'][lab][mod]['amp']))
        

                    
            
        
# ===== save collated data =====
folder = 'Data/'
name = '{}_HFI[{}]+{}_modulation.json'.format(cell_type,HFI,HFI_delay)
cf.save_data(data,folder+name)
print('Saving data as {}'.format(name))


'''
h.quit()
'''

    