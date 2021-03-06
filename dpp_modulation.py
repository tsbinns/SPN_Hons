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
model_iterator    = list(range(specs[cell_type]['N']))  # range(specs[cell_type]['N']) gives all models; must be a list for saving data!
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
        run_info = {'curr_n':cell_n, 'tot_n':len(model_iterator)}
        pc.submit(sf.ACh_modulation, model_data, stim_data, cell_index, \
            run_info, mod_factors)
        
    while pc.working(): # gather results
        data = pc.pyret()
        # save file to folder
        keys_1 = list(data.keys())
        keys_2 = list(data[keys_1[0]].keys())
        name = '{}_{}_modulation'.format(cell_type,data[keys_1[0]][keys_2[0]]['id'])
        cf.save_data(data,'{}/{}.json'.format(folder,name))
    
    
pc.done() # end parallelisation


# for timing simulations
print('Simulations completed (took %.0f secs).\nNow performing calculations/collating data...' % (time.time()-start))



# ===== combine and/or average data =====

data_avg = {}

if collate: # collates data if loading from files
    
    data_all = {}
    
    for n, i in enumerate(model_iterator):
        # load data
        name = '{}_{}_modulation'.format(cell_type,i)
        data = cf.load_data('{}/{}.json'.format(folder,name))
        # combine data
        data_all[i] = data
        
    data = data_all
    

# averages data
for i, clus_lab in enumerate(stim_data['clustered']['label']): # for each clustered stimulation target
    
    
    # data for cholinergic modulation of same site as for clustered input
    data_avg[clus_lab] = {clus_lab:[]}
    data_avg[clus_lab][clus_lab] = {'vm':[], 'avg_vm':[], 'rheo':[],
        'avg_rheo':[], 'dur':[], 'avg_dur':[], 'amp':[], 'avg_amp':[]}
    
    for cell_index in model_iterator:
            
        data_avg[clus_lab][clus_lab]['vm'].append(data[cell_index][clus_lab][clus_lab]['vm'])
        data_avg[clus_lab][clus_lab]['rheo'].append(data[cell_index][clus_lab][clus_lab]['rheo'])
        data_avg[clus_lab][clus_lab]['dur'].append(data[cell_index][clus_lab][clus_lab]['dur'])
        data_avg[clus_lab][clus_lab]['amp'].append(data[cell_index][clus_lab][clus_lab]['amp'])
    
    data_avg[clus_lab][clus_lab]['avg_vm'] = np.ndarray.tolist(np.mean(data_avg[clus_lab][clus_lab]['vm'],axis=0))
    data_avg[clus_lab][clus_lab]['avg_rheo'] = float(np.mean(data_avg[clus_lab][clus_lab]['rheo']))
    data_avg[clus_lab][clus_lab]['avg_dur'] = float(np.mean(data_avg[clus_lab][clus_lab]['dur']))
    data_avg[clus_lab][clus_lab]['avg_amp'] = float(np.mean(data_avg[clus_lab][clus_lab]['amp']))
    
    
    # data for cholinergic modulation of off-site and soma
    for j, ACh_lab in enumerate(stim_data['ACh']['label']):
    
        data_avg[clus_lab][ACh_lab] = {'vm':[], 'avg_vm':[], 'rheo':[],
            'avg_rheo':[], 'dur':[], 'avg_dur':[], 'amp':[], 'avg_amp':[]}
        
        for cell_index in model_iterator: # for each simulated cell
            
            data_avg[clus_lab][ACh_lab]['vm'].append(data[cell_index][clus_lab][ACh_lab]['vm'])
            data_avg[clus_lab][ACh_lab]['rheo'].append(data[cell_index][clus_lab][ACh_lab]['rheo'])
            data_avg[clus_lab][ACh_lab]['dur'].append(data[cell_index][clus_lab][ACh_lab]['dur'])
            data_avg[clus_lab][ACh_lab]['amp'].append(data[cell_index][clus_lab][ACh_lab]['amp'])
        
        data_avg[clus_lab][ACh_lab]['avg_vm'] = np.ndarray.tolist(np.mean(data_avg[clus_lab][ACh_lab]['vm'],axis=0))
        data_avg[clus_lab][ACh_lab]['avg_rheo'] = float(np.mean(data_avg[clus_lab][ACh_lab]['rheo']))
        data_avg[clus_lab][ACh_lab]['avg_dur'] = float(np.mean(data_avg[clus_lab][ACh_lab]['dur']))
        data_avg[clus_lab][ACh_lab]['avg_amp'] = float(np.mean(data_avg[clus_lab][ACh_lab]['amp']))
    
    
# general simulation info
data_avg['meta'] = {'cell_type':cell_type, 'specs':model_iterator,
    'tm':data[cell_index][clus_lab][ACh_lab]['tm'], 'clus':[], 'ACh':[]}


# clustered input-specific info
data_avg['meta']['clus'] = {
    'dist':[data[cell_index]['proximal dend']['proximal dend']['clust_dist'],
            data[cell_index]['distal dend']['distal dend']['clust_dist']],
    'stim_n':stim_data['clustered']['params']['stim_n'],
    'isi':stim_data['clustered']['params']['isi'],
    'stim_t':stim_data['clustered']['params']['stim_t'],
    'stop_t':stim_data['clustered']['params']['stop_t'],
    'pre_t':stim_data['clustered']['params']['pre_t'],
    'labels':stim_data['clustered']['label'],
    'targets':stim_data['clustered']['target']}


# cholinergic input-specific info
data_avg['meta']['ACh'] = {
    'dist':[data[cell_index]['distal dend']['off-site']['ACh_dist'],data[cell_index]['distal dend']['soma']['ACh_dist']], 
    'stim_t':stim_data['ACh']['params']['stim_t'],
    'stop_t':stim_data['ACh']['params']['stop_t'],
    'labels':stim_data['ACh']['label'],
    'targets':stim_data['ACh']['target']}
    
    
   
# ===== save collated data =====
folder = 'Data/'
name = '{}_n{}_modulation.json'.format(cell_type,stim_data['clustered']['params']['stim_n'])
cf.save_data(data_avg,folder+name)
print('Saving data as {}'.format(name))



h.quit()



    