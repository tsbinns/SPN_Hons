
from   neuron               import h
import common_functions     as cf
import MSN_builder          as build




def dpp_validation(model_data, stim_data, cell_index):
    '''
    Provides clustered glutamatergic input to SPNs to generate the dendritic 
    plateau potential
    '''
    
    data = {}
    tm = {}
    vm = {}
    for i, t in enumerate(model_data['target']):
        # initiate cell
        cell = build.MSN(  params=model_data['specs'][model_data['cell_type']]['par'],
                           morphology=model_data['specs'][model_data['cell_type']]['morph'],
                           variables=model_data['model_sets'][cell_index]['variables']   )
        
        rheobase        =   model_data['model_sets'][cell_index]['rheobase']
        
        # record vectors
        tm[i]  = h.Vector()
        tm[i].record(h._ref_t)
        vm[i]  = h.Vector()
        vm[i].record(cell.soma(0.5)._ref_v)
        # add clustered inputs
        syn,stim,ncon,d2soma = cf.set_clustered_stim(cell,t,n=stim_data['stim_n'], \
                                                     act_time=stim_data['stim_t'], \
                                                     ISI=stim_data['isi'])
        # run simulation
        h.finitialize(-80)
        while h.t < stim_data['stop_t']:
            h.fadvance()
        tm[i] = tm[i].to_python()
        vm[i] = vm[i].to_python()
        # calculate dpp duration and amplitude
        base_t = tm[i].index(min(tm[i], key=lambda x:abs(x-stim_data['stim_t'])))-1
        dur = cf.dpp_dur(tm[i],vm[i],vm[i][base_t],stim_data['stim_t'])
        amp = cf.dpp_amp(tm[i],vm[i],vm[i][base_t],stim_data['stim_t'])
        # store data
        data[i] = {'tm':tm[i], 'vm':vm[i], 'dist':d2soma, 'rheo':rheobase, \
                   'dur':dur, 'amp':amp, 'id':int(cell_index), \
                   'cell_type':model_data['cell_type']}
        
    return data