'''
Functions for running simulations; useful for 'bulletin board' parallelisation

Thomas Binns (author), 28/01/21
'''


from   neuron               import h
import common_functions     as cf
import MSN_builder          as build




def dpp_validation(model_data, stim_data, cell_index, sim_info):
    '''
    Provides clustered glutamatergic input to SPNs to generate the dendritic 
    plateau potential.
    
    INPUT(S):
        - model_data: model paramaters (specification, cell type, model 
            sets, and stimulation targets)
        - stim_data: stimulation parameters (number of inputs, time of inputs,
            input inter-spike intervals, time to stop simulating)
        - cell_index: cell specification being simulated
        - sim_info: information about the simulation run (this simulation and 
            total number of simulations)
        
    OUTPUT(S):
        - data: simulated data including: simulation times; simulated voltages;
            distance of stimulated target to the soma; rheobase of the model; 
            half-max width of the potential; peak amplitude of the potential;
            cell specification being simulated; cell type being simulated
            
    Thomas Binns (author), 26/01/21
    '''
    
    # print simulation info to check on process
    print('Simulating cell specification {} of {}'.format( \
          sim_info['curr_n']+1,sim_info['tot_n']),flush = True)
    
    data = {}
    for i, t in enumerate(model_data['target']): # for each simulation target
        # initiate cell
        cell = build.MSN(params=model_data['specs'][model_data['cell_type']]['par'],
                         morphology=model_data['specs'][model_data['cell_type']]['morph'],
                         variables=model_data['model_sets'][cell_index]['variables'])
        rheobase = model_data['model_sets'][cell_index]['rheobase']
        
        # record vectors
        tm = h.Vector()
        tm.record(h._ref_t)
        vm = h.Vector()
        vm.record(cell.soma(0.5)._ref_v)
        
        # add clustered inputs
        syn,stim,ncon,d2soma = cf.set_clustered_stim(cell,t,n=stim_data['stim_n'], \
                                                     act_time=stim_data['stim_t'], \
                                                     ISI=stim_data['isi'])
        # run simulation
        h.finitialize(-80)
        while h.t < stim_data['stop_t']:
            h.fadvance()
        tm = tm.to_python()
        vm = vm.to_python()
        
        # calculate dpp duration and amplitude
        base_t = tm.index(min(tm, key=lambda x:abs(x-stim_data['stim_t'])))-1
        dur = cf.dpp_dur(tm,vm,vm[base_t],stim_data['stim_t'])
        amp = cf.dpp_amp(tm,vm,vm[base_t],stim_data['stim_t'])
        
        # collate data
        data[i] = {'tm':tm, 'vm':vm, 'dist':d2soma, 'rheo':rheobase, \
                   'dur':dur, 'amp':amp, 'id':int(cell_index), \
                   'cell_type':model_data['cell_type']}
        
    return data





