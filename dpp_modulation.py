'''
Applies cholinergic modulation to SPNs alongside clustered inputs to examine 
effects on the plateau potential

Thomas Binns (author), 29/01/21
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
stim_info = cf.params_for_input(cell_type, 'clustered+Ach')
target = stim_info['clustered']['target']
target_labels = stim_info['clustered']['label']

# open library (channel distributions etc)   
with open(specs[cell_type]['lib'], 'rb') as f:
    model_sets = pickle.load(f, encoding="latin1")



# simulate model(s) ======================
    
    
    
    
    
    
    
    
    
    