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

cell_type = 'dspn'
if cell_type != 'dspn' and cell_type != 'ispn':
    raise ValueError("The requested cell type is not supported.\nOnly 'dpsn' and 'ispn' are recognised.")
    
model_iterator = list(range(1)) # use model_iterator = range(specs[cell_type]['N']) for all models
# for dspn, 10 has lowest rheo, 54 has highest, 22 has mean, 41 has median; 22 is also average for experimental value
# for ispn, 8 has mean and median; 1 is average for experimental value

if pc.id() == 0:
    print('Simulating {} cell iteration(s) of type: {}'.format(len(model_iterator),cell_type), \
          flush=True)


# stimulation details
stim_info = cf.params_for_input(cell_type, 'clustered')
target = stim_info['clustered']['target']
target_labels = stim_info['clustered']['label']
stim_data = stim_info['clustered']['params']

# open library (channel distributions etc)   
with open(specs[cell_type]['lib'], 'rb') as f:
    model_sets = pickle.load(f, encoding="latin1")
    

# ===== simulate model(s) =====
    
data = {}

# model information to pass to simulations
model_data = {'specs':specs[cell_type], 'cell_type':cell_type, 'model_sets':model_sets, \
              'target':target, 'target_labels':target_labels}

start = time.time() # for timing simulations

pc.runworker() # start workers for parallelisation


if pc.nhost() == 1: # use the serial form
    
    collate = 0
    
    for cell_n, cell_index in enumerate(model_iterator): # for each model
        # simulate model
        run_info = {'curr_n':cell_n, 'tot_n':len(model_iterator)}
        data[cell_index] = sf.dpp_validation(model_data, stim_data, \
                          cell_index, run_info)
            
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
        data = pc.pyret()
        # save file to folder
        keys = list(data.keys())
        name = '{}_{}_validation'.format(data[keys[0]]['cell_type'],data[keys[0]]['id'])
        cf.save_data(data,'{}/{}.json'.format(folder,name))
  
  
pc.done() # end parallelisation


# for timing simulations
print('Simulations completed (took %.0f secs).\nNow performing calculations/collating data...' % (time.time()-start))



# =====combine and/or average data =====

data_avg = {}

if collate: # collates data if loading from files
    
    data_all = {}
    
    for n, i in enumerate(model_iterator):
        # load data
        name = '{}_{}_validation'.format(cell_type,i)
        data = cf.load_data('{}/{}.json'.format(folder,name))
        # combine data
        data_all[i] = data
        
    data = data_all
        

# averages data    
for i, lab in enumerate(target_labels):
    
    data_avg[lab] = {'vm':[], 'avg_vm':[], 'rheo':[], 'avg_rheo':[], \
                      'dur':[], 'avg_dur':[], 'amp':[], 'avg_amp':[]}
    
    for cell_index in model_iterator:
        
        data_avg[lab]['vm'].append(data[cell_index][lab]['vm'])
        data_avg[lab]['rheo'].append(data[cell_index][lab]['rheo'])
        data_avg[lab]['dur'].append(data[cell_index][lab]['dur'])
        data_avg[lab]['amp'].append(data[cell_index][lab]['amp'])
        
    data_avg[lab]['avg_vm'] = np.ndarray.tolist(np.mean(data_avg[lab]['vm'],axis=0))
    data_avg[lab]['avg_rheo'] = float(np.mean(data_avg[lab]['rheo']))
    data_avg[lab]['avg_dur'] = float(np.mean(data_avg[lab]['dur']))
    data_avg[lab]['avg_amp'] = float(np.mean(data_avg[lab]['amp']))
    
data_avg['meta'] = {'cell_type':cell_type, 'tm': data[cell_index][lab]['tm'], \
                    'dist': [data[cell_index][target_labels[0]]['dist'],data[cell_index][target_labels[1]]['dist']], \
                    'stim_n':stim_data['stim_n'], 'isi':stim_data['isi'], \
                    'stim_t':stim_data['stim_t'], 'stop_t':stim_data['stop_t'], \
                    'pre_t':stim_data['pre_t'], 'labels': target_labels, \
                    'targets':target, 'specs':model_iterator}        



# ===== save collated data =====
folder = 'Data/'
name = '{}_n{}_validation.json'.format(cell_type,stim_data['stim_n'])
cf.save_data(data_avg,folder+name)
print('Saving data as {}'.format(name))



h.quit()