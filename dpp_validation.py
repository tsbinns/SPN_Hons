'''
Basic simulation of the dendritic plateau potential without any modulation.
    - Run from command line with e.g. mpiexec -n 4 python dpp_validation.py
    
Thomas Binns (author), 27/01/21
'''



from   neuron           import h
import numpy                as np
import pickle
import common_functions     as cf
import simulation_functions as sf
import time



# for parallelisation
h.nrnmpi_init()
pc = h.ParallelContext()

# Load model mechanisms
import neuron               as nrn
nrn.load_mechanisms('mechanisms/single')

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
cell_type         = 'ispn'    # 'dspn'/'ispn'
model_iterator    = list(range(5))  # range(specs[cell_type]['N']) gives all models
# for dspn, 10 has lowest rheo, 54 has highest, 22 has mean, 41 has median; 22 is also average for experimental value
# for ispn, 8 has mean and median; 1 is average for experimental value
if pc.id() == 0:
    print('Simulating {} cell(s) of type: {}'.format(len(model_iterator),cell_type))


# stimulation details
stim_info = cf.params_for_input(cell_type, 'clustered')
target = stim_info['clustered']['target']
target_labels = stim_info['clustered']['label']

# open library (channel distributions etc)   
with open(specs[cell_type]['lib'], 'rb') as f:
    model_sets = pickle.load(f, encoding="latin1")



# simulate model(s) ======================
    
OUT = {}

# model information to pass to simulations
model_data = {'specs':specs, 'cell_type':cell_type, 'model_sets':model_sets, \
              'target':target}

# stimulation information to pass to simulations
stim_data = {'stim_n':16, 'stim_t':100, 'isi':1, 'pre_t':-50}
stim_data['stop_t'] = stim_data['stim_t'] + 250

start = time.time() # for timing simulations

pc.runworker() # start workers for parallelisation


if pc.nhost() == 1: # use the serial form
    
    collate = 0
    
    for cell_n, cell_index in enumerate(model_iterator): # for each model
        # simulate model
        print('Simulating cell specification {} of {}'.format(cell_n+1,len(model_iterator)))
        OUT[cell_index] = sf.dpp_validation(model_data, stim_data, cell_index)
            
else: # use the bulleting board form
    
    collate = 1
    
    # clear temp data folder
    folder = 'temp_data'
    cf.clear_folder(folder,'.json')
    
    for cell_n, cell_index in enumerate(model_iterator): # scatter processes
        # simulate model
        sim_info = {'curr_n':cell_n, 'tot_n':len(model_iterator)}
        pc.submit(sf.dpp_validation,model_data,stim_data,cell_index,sim_info)
        #sys.stdout.flush()
        
    while pc.working(): # gather results
        OUT = pc.pyret()
        # save file to folder
        name = '{}_{}'.format(OUT[0]['cell_type'],OUT[0]['id'])
        cf.save_data(OUT,'{}/{}.json'.format(folder,name))
  
  
pc.done() # end parallelisation


# for timing simulations
end = time.time()
print('Simulations completed (took %.0f secs), now performing calculations...' % (end-start))



# combine and/or average data ===========

OUT_avg = {}

if collate:
    
    OUT_all = {}
    
    for n, i in enumerate(model_iterator):
        # load data
        name = '{}_{}'.format(cell_type,i)
        OUT = cf.load_data('{}/{}.json'.format(folder,name))
        # combine data
        OUT_all[i] = OUT
        
    OUT = OUT_all
    
    for i in range(len(target)): # collates and averages data
        
        idx = str(i)
        OUT_avg[i] = {'vm':[], 'avg_vm':[], 'rheo':[], 'avg_rheo':[], \
                          'dur':[], 'avg_dur':[], 'amp':[], 'avg_amp':[]}
        
        for cell_index in model_iterator:
            
            OUT_avg[i]['vm'].append(OUT[cell_index][idx]['vm'])
            OUT_avg[i]['rheo'].append(OUT[cell_index][idx]['rheo'])
            OUT_avg[i]['dur'].append(OUT[cell_index][idx]['dur'])
            OUT_avg[i]['amp'].append(OUT[cell_index][idx]['amp'])
            
        OUT_avg[i]['avg_vm'] = np.ndarray.tolist(np.mean(OUT_avg[i]['vm'],axis=0))
        OUT_avg[i]['avg_rheo'] = float(np.mean(OUT_avg[i]['rheo']))
        OUT_avg[i]['avg_dur'] = float(np.mean(OUT_avg[i]['dur']))
        OUT_avg[i]['avg_amp'] = float(np.mean(OUT_avg[i]['amp']))
        
    OUT_avg['meta'] = {'cell_type':cell_type, 'tm': OUT[cell_index]['0']['tm'], \
                       'dist': [OUT[cell_index]['0']['dist'],OUT[cell_index]['1']['dist']], \
                       'stim_n':stim_data['stim_n'], 'isi':stim_data['isi'], \
                       'stim_t':stim_data['stim_t'], 'stop_t':stim_data['stop_t'], \
                       'pre_t':stim_data['pre_t'], 'labels': target_labels, \
                       'targets':target, 'specs':model_iterator} 
        
else:
    
    for idx in range(len(target)): # collates and averages data
        
        OUT_avg[idx] = {'vm':[], 'avg_vm':[], 'rheo':[], 'avg_rheo':[], \
                          'dur':[], 'avg_dur':[], 'amp':[], 'avg_amp':[]}
        
        for cell_index in model_iterator:
            
            OUT_avg[idx]['vm'].append(OUT[cell_index][idx]['vm'])
            OUT_avg[idx]['rheo'].append(OUT[cell_index][idx]['rheo'])
            OUT_avg[idx]['dur'].append(OUT[cell_index][idx]['dur'])
            OUT_avg[idx]['amp'].append(OUT[cell_index][idx]['amp'])
            
        OUT_avg[idx]['avg_vm'] = np.ndarray.tolist(np.mean(OUT_avg[idx]['vm'],axis=0))
        OUT_avg[idx]['avg_rheo'] = float(np.mean(OUT_avg[idx]['rheo']))
        OUT_avg[idx]['avg_dur'] = float(np.mean(OUT_avg[idx]['dur']))
        OUT_avg[idx]['avg_amp'] = float(np.mean(OUT_avg[idx]['amp']))
        
    OUT_avg['meta'] = {'cell_type':cell_type, 'tm': OUT[cell_index][0]['tm'], \
                       'dist': [OUT[cell_index][0]['dist'],OUT[cell_index][1]['dist']], \
                       'stim_n':stim_data['stim_n'], 'isi':stim_data['isi'], \
                       'stim_t':stim_data['stim_t'], 'stop_t':stim_data['stop_t'], \
                       'pre_t':stim_data['pre_t'], 'labels': target_labels, \
                       'targets':target, 'specs':model_iterator}        



# save collated data =================
name = '{}_n{}.json'.format(cell_type,stim_data['stim_n'])
cf.save_data(OUT_avg,name) # save data
print('Saving data as {}'.format(name))



h.quit()