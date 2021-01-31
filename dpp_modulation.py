'''
Applies cholinergic modulation to SPNs alongside clustered inputs to examine 
effects on the plateau potential
    - WIP!!

Thomas Binns (author), 29/01/21
'''


from   neuron           import h
import numpy                as np
import pickle
import common_functions     as cf
import simulation_functions as sf
import time
import matplotlib.pyplot as plt



# for parallelisation
h.nrnmpi_init()
pc = h.ParallelContext()

# Load model mechanisms
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
cell_type         = 'ispn'    # 'dspn'/'ispn'
model_iterator    = list(range(1))  # range(specs[cell_type]['N']) gives all models; must be a list for saving data!
# for dspn, 10 has lowest rheo, 54 has highest, 22 has mean, 41 has median; 22 is also average for experimental value
# for ispn, 8 has mean and median; 1 is average for experimental value
if pc.id() == 0:
    print('Simulating {} cell(s) of type: {}'.format(len(model_iterator),cell_type))
    
# stimulation details
stim_data = cf.params_for_input(cell_type, 'ACh')

# open library (channel distributions etc)   
with open(specs[cell_type]['lib'], 'rb') as f:
    model_sets = pickle.load(f, encoding="latin1")

# model information to pass to simulations
model_data = {'specs':specs[cell_type], 'cell_type':cell_type, 'model_sets':model_sets}


# ===== gets cholinergic modulation factors =====

mod_factors = cf.draw_factors_ACh(cell_type, mode='mean')



# ===== simulate model(s) =====
    
start = time.time() # for timing simulations

pc.runworker() # start workers for parallelisation

data = {}
if pc.nhost() == 1: # use the serial form
    
    collate = 0
    
    for cell_n, cell_index in enumerate(model_iterator): # for each model
        # simulate model
        run_info = {'curr_n':cell_n, 'tot_n':len(model_iterator)}
        data[cell_index] = sf.ACh_modulation(model_data, stim_data, \
            cell_index, run_info, mod_factors)
            
else: # use the bulleting board form
    
    collate = 1
    
    # clear temp data folder
    folder = 'temp_data'
    cf.clear_folder(folder,'.json')
    
    for cell_n, cell_index in enumerate(model_iterator): # scatter processes
        # simulate model
        sim_info = {'curr_n':cell_n, 'tot_n':len(model_iterator)}
        pc.submit(sf.ACh_modulation, model_data, stim_data, cell_index, \
            run_info, mod_factors)
        
    while pc.working(): # gather results
        data = pc.pyret()
        # save file to folder
        keys_1 = data.keys()
        keys_2 = data[keys_1[0]].keys()
        name = '{}_{}'.format(cell_type,data[keys_1[0]][keys_2[0]]['id'])
        cf.save_data(data,'{}/{}.json'.format(folder,name))
    
    
pc.done() # end parallelisation


# for timing simulations
end = time.time()
print('Simulations completed (took %.0f secs).\nNow performing calculations/collating data...' % (end-start))
    
for i in data:
    for j, key_1 in enumerate(data[i]):
        plt.figure()
        for k, key_2 in enumerate(data[i][key_1]):
            plt.plot(data[i][key_1][key_2]['tm'],data[i][key_1][key_2]['vm'],label=key_2)
        plt.legend()
    
    