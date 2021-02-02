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

if pc.id() == 0:
    print('Simulating {} cell iteration(s) of type: {}'.format(len(model_iterator),cell_type), \
          flush=True)
   
# open library (channel distributions etc)   
with open(specs[cell_type]['lib'], 'rb') as f:
    model_sets = pickle.load(f, encoding="latin1")



# ===== simulate model(s) =====
# model information to pass to simulations
model_data = {'specs':specs[cell_type], 'cell_type':cell_type, 'model_sets':model_sets}

start = time.time() # for timing simulations

pc.runworker() # start workers for parallelisation


if pc.nhost() == 1: # use the serial form
    
    collate = 0
    data = {}
    for cell_n, cell_index in enumerate(model_iterator): # for each model
        # simulate model
        run_info = {'curr_n':cell_n, 'tot_n':len(model_iterator)}
        data[cell_index] = sf.dpp_generation(model_data, cell_index, run_info, noise=True, HFI=False)
            
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

    
    
    
    