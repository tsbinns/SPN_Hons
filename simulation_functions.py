'''
Functions for running simulations; useful for 'bulletin board' parallelisation

Thomas Binns (author), 28/01/21
'''


from   neuron           import h
import common_functions      as cf
import MSN_builder           as build
import numpy                 as np
import modulation_lib        as modulate




def dpp_validation(model_data,
                   stim_data,
                   cell_index,
                   run_info,
                   dur_and_amp = True):
    '''
    Provides clustered glutamatergic input to SPNs to generate the dendritic 
    plateau potential in the absence of modulation.
    
    INPUT(S):
        - model_data: model paramaters (specification, cell type, model 
            sets, and stimulation targets) [dict]
        - stim_data: stimulation parameters (number of inputs, time of inputs,
            input inter-spike intervals, time to stop simulating) [dict]
        - cell_index: cell specification being simulated [int]
        - run_info: information about the simulation run (this simulation and 
            total number of simulations) [dict]
        - dur_and_amp: whether to calculate the duration and peak amplitude of
            the plateau potential (default True) [bool]
        
    OUTPUT(S):
        - data: simulated data including: simulation times; simulated voltages;
            distance of stimulated target to the soma; rheobase of the model; 
            half-max width of the potential; peak amplitude of the potential;
            cell specification being simulated; cell type being simulated
            [dict]
        - ncon: NetCon object
            
    Thomas Binns (author), 26/01/21
    '''
    
    # ===== print simulation info to monitor progress =====
    print('Simulating cell specification {} of {}'.format( \
          run_info['curr_n']+1,run_info['tot_n']),flush = True)
    
    
    # ===== simulation =====
    data = {}
    for i, t in enumerate(model_data['target']): # for each simulation target
        
        clus_lab = model_data['target_labels'][i] # label for input target
        
        # initiate cell
        cell = build.MSN(params=model_data['specs']['par'],
                         morphology=model_data['specs']['morph'],
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
        
        # collate data
        data[clus_lab] = {'tm':tm, 'vm':vm, 'dist':d2soma, 'rheo':rheobase, \
            'id':int(cell_index), 'cell_type':model_data['cell_type']}
        
        if dur_and_amp:
            # calculate dpp duration and amplitude
            base_t = tm.index(min(tm, key=lambda x:abs(x-stim_data['stim_t'])))-1
            data[clus_lab]['dur'] = \
                cf.dpp_dur(tm,vm,vm[base_t],stim_data['stim_t'])
            data[clus_lab]['amp'] = \
                cf.dpp_amp(tm,vm,vm[base_t],stim_data['stim_t'])
            
        
    return data


        

def dpp_generation(model_data,
                   cell_index,
                   run_info,
                   noise = True,
                   HFI = False,
                   HFI_delay = 0,
                   dur_and_amp = True,
                   spike = False):
    '''
    Provides clustered glutamatergic input to SPNs to generate the dendritic 
    plateau potential in the absence of modulation.
    
    INPUT(S):
        - model_data: model paramaters (specification, cell type, model 
            sets, and stimulation targets) [dict]
        - stim_data: stimulation parameters (number of inputs, time of inputs,
            input inter-spike intervals, time to stop simulating) [dict]
        - cell_index: cell specification being simulated [int]
        - run_info: information about the simulation run (this simulation and 
            total number of simulations) [dict]
        - dur_and_amp: whether to calculate the duration and peak amplitude of
            the plateau potential (default True) [bool]
        
    OUTPUT(S):
        - data: simulated data including: simulation times; simulated voltages;
            distance of stimulated target to the soma; rheobase of the model; 
            half-max width of the potential; peak amplitude of the potential;
            cell specification being simulated; cell type being simulated
            [dict]
        - ncon: NetCon object
            
    Thomas Binns (author), 26/01/21
    '''
    
    # ===== print simulation info to monitor progress =====
    print('Simulating cell specification {} of {}'.format( \
          run_info['curr_n']+1,run_info['tot_n']),flush = True)
    
    
    # ===== gets stimulation info =====
    # clustered inputs
    clus_info = cf.params_for_input(model_data['cell_type'], 'clustered')
    clus_params = clus_info['clustered']['params']
    
    # noise inputs
    if noise:
        noise_info = cf.params_for_input(model_data['cell_type'], 'noise')
        noise_params = noise_info['noise']['params']
    
    # HFIs
    if HFI:
        HFI_info = cf.params_for_input(model_data['cell_type'], 'HFI')
        HFI_params = HFI_info['HFI']['params']
        
    # ===== simulation =====
    data = {}
    
    for i, tar in enumerate(clus_info['clustered']['target']): # for each simulation target
        
        clus_lab = clus_info['clustered']['label'][i] # label for input target
        
        
        # initiate cell
        cell = build.MSN(params=model_data['specs']['par'],
                         morphology=model_data['specs']['morph'],
                         variables=model_data['model_sets'][cell_index]['variables'])
        rheobase = model_data['model_sets'][cell_index]['rheobase']
        
        
        # record vectors
        tm = h.Vector()
        tm.record(h._ref_t)
        vm = h.Vector()
        vm.record(cell.soma(0.5)._ref_v)
        
        
        # add clustered inputs
        clus_syn, clus_stim, clus_ncon, clus_d2soma = cf.set_clustered_stim(cell, tar ,n=clus_params['stim_n'],
            act_time=clus_params['stim_t'], ISI=clus_params['isi'])
        
        # add background noise
        if noise:
            '''
            noise_syn, noise_stim, noise_ncon = cf.set_noise(cell, model_data['cell_type'], 
                freq_glut=noise_params['freq_glut'], freq_gaba=noise_params['freq_gaba'],
                n_glut=noise_params['n_glut'], n_gaba=noise_params['n_gaba'], only_dend=noise_params['only dend'],
                glut_delay=noise_params['stim_t'], gaba_delay=noise_params['stim_t'])
            '''
            noise_syn, noise_stim, noise_ncon = cf.set_bg_noise(cell,model_data['cell_type'], fglut=noise_params['freq_glut'],
                fgaba=noise_params['freq_gaba'],dendOnly=noise_params['only dend'])

        
        # add high-frequency inputs
        if HFI:
            # adds HFI
            HFI_syn, HFI_stim, HFI_ncon, arrangement = cf.set_HFI(cell, model_data['cell_type'], 
                freq=HFI_params['freq'], n_inputs=HFI_params['n_inputs'],
                delay=clus_params['stim_t']+(clus_params['stim_n']*clus_params['isi'])+HFI_delay, 
                exclude=HFI_info['HFI']['exclude'])
            # collates data
            data['HFI'] = arrangement
        
        
        # run simulation
        h.finitialize(-80)
        while h.t < clus_params['stop_t']+HFI_delay+(clus_params['stim_n']*clus_params['isi']):
            h.fadvance()
        tm = tm.to_python()
        vm = vm.to_python()
        
        
        # collate data
        data[clus_lab] = {'vm':vm}
        data['meta'] = {'id':int(cell_index), 'round':run_info['round'], 'cell_type':model_data['cell_type'], 'tm':tm, 
                        'rheo':rheobase}
        
        
        # calculate dpp duration and amplitude
        if dur_and_amp:
            base_t_start = tm.index(min(tm, key=lambda x:abs(x-(clus_params['stim_t']+clus_params['pre_t']))))
            base_t_end = tm.index(min(tm, key=lambda x:abs(x-clus_params['stim_t'])))
            base_vm = np.mean(vm[base_t_start:base_t_end])
            data[clus_lab]['dur'] = cf.dpp_dur(tm, vm, base_vm, clus_params['stim_t'])
            data[clus_lab]['amp'] = cf.dpp_amp(tm, vm, base_vm, clus_params['stim_t'])
            
        # get spike-related data
        if spike:
            spiked = 0
            thresh = 0
            # checks whether spike occured
            data[clus_lab]['spiked'] = 0
            if max(vm) > thresh:
                data[clus_lab]['spiked'] = 1
                spiked = 1
            # checks time of first spike number of spikes that occured (if any)
            data[clus_lab]['first_spike'] = []
            data[clus_lab]['spike_n'] = []
            if spiked == 1:
                data[clus_lab]['first_spike'] = tm[next(i for i, x in enumerate(vm) if x > 0)]
                data[clus_lab]['spike_n'] = cf.spike_n(vm)
            
        
    return data





def dpp_ACh_modded(model_data,
                   cell_index,
                   run_info,
                   noise = True,
                   HFI = False,
                   HFI_delay = 0,
                   dur_and_amp = True,
                   spike = False,
                   mod_tar = 'all'):
    '''
    Provides clustered glutamatergic input to SPNs to generate the dendritic 
    plateau potential in the absence of modulation.
    
    INPUT(S):
        - model_data: model paramaters (specification, cell type, model 
            sets, and stimulation targets) [dict]
        - stim_data: stimulation parameters (number of inputs, time of inputs,
            input inter-spike intervals, time to stop simulating) [dict]
        - cell_index: cell specification being simulated [int]
        - run_info: information about the simulation run (this simulation and 
            total number of simulations) [dict]
        - dur_and_amp: whether to calculate the duration and peak amplitude of
            the plateau potential (default True) [bool]
        
    OUTPUT(S):
        - data: simulated data including: simulation times; simulated voltages;
            distance of stimulated target to the soma; rheobase of the model; 
            half-max width of the potential; peak amplitude of the potential;
            cell specification being simulated; cell type being simulated
            [dict]
        - ncon: NetCon object
            
    Thomas Binns (author), 26/01/21
    '''
    
    # ===== print simulation info to monitor progress =====
    print('Simulating cell specification {} of {}'.format( \
          run_info['curr_n']+1,run_info['tot_n']),flush = True)
    
    
    # ===== gets stimulation info =====
    # clustered inputs
    clus_info = cf.params_for_input(model_data['cell_type'], 'clustered')
    clus_params = clus_info['clustered']['params']
    
    
    # modulation inputs
    mod_factors = cf.draw_factors_ACh(model_data['cell_type'], mode='random')
    ACh_info = cf.params_for_input(model_data['cell_type'], 'ACh')
    ACh_params = ACh_info['ACh']['params']
    
    # gets vector for modulation timing
    state = [1 if ht >= ACh_params['stim_t'] and ht <= ACh_params['stop_t'] \
             else 0 for ht in np.arange(0,clus_params['stop_t'],h.dt)]
    if 'kaf' in mod_factors:
        kaf_state = [x * mod_factors['kaf'] for x in state]
        kaf_state = h.Vector(kaf_state)
    state = h.Vector(state)
    
    # collates time-dependent mechanism modulation scaling
    mech_scale = {}
    for key in mod_factors:
        if key == 'kaf':
            mech_scale[key] = kaf_state
        else:
            mech_scale[key] = state
    
    
    # noise inputs
    if noise:
        noise_info = cf.params_for_input(model_data['cell_type'], 'noise')
        noise_params = noise_info['noise']['params']
    
    # HFIs
    if HFI:
        HFI_info = cf.params_for_input(model_data['cell_type'], 'HFI')
        HFI_params = HFI_info['HFI']['params']
        
    # ===== simulation =====
    data = {}
    
    for i, clus_tar in enumerate(clus_info['clustered']['target']): # for each simulation target
        
        clus_lab = clus_info['clustered']['label'][i] # label for input target
        data[clus_lab] = {}
        
        # get the targets for cholinergic input
        ACh_targets = {}
        ACh_targets['target'] = [clus_tar]
        ACh_targets['label'] = [clus_lab]
        for j, ACh_t in enumerate(ACh_info['ACh']['target']):
            ACh_targets['target'].append(ACh_t)
            ACh_targets['label'].append(ACh_info['ACh']['label'][j])
        target_x = .5
        
        if mod_tar == 'all':
            ACh_targets['target'] = ['all']
            ACh_targets['label'] = ['all']
            #target_x = 'all'
        
        for j, ACh_t in enumerate(ACh_targets['target']): # for each cholinergic input target
            
            ACh_lab = ACh_targets['label'][j]
            
            # initiate cell
            cell = build.MSN(params=model_data['specs']['par'],
                             morphology=model_data['specs']['morph'],
                             variables=model_data['model_sets'][cell_index]['variables'])
            rheobase = model_data['model_sets'][cell_index]['rheobase']
            
            # record vectors
            tm = h.Vector()
            tm.record(h._ref_t)
            vm = h.Vector()
            vm.record(cell.soma(0.5)._ref_v)
            
            
            # add clustered inputs
            clus_syn, clus_stim, clus_ncon, clus_d2soma = cf.set_clustered_stim(cell, clus_tar ,n=clus_params['stim_n'],
                act_time=clus_params['stim_t'], ISI=clus_params['isi'])
            
            # add background noise
            if noise:
                noise_syn, noise_stim, noise_ncon = cf.set_bg_noise(cell, model_data['cell_type'], 
                    fglut=noise_params['freq_glut'], fgaba=noise_params['freq_gaba'], dendOnly=noise_params['only dend'])
    
            
            # add high-frequency inputs
            if HFI:
                # adds HFI
                HFI_syn, HFI_stim, HFI_ncon, arrangement = cf.set_HFI(cell, model_data['cell_type'], 
                    freq=HFI_params['freq'], n_inputs=HFI_params['n_inputs'],
                    delay=clus_params['stim_t']+(clus_params['stim_n']*clus_params['isi'])+HFI_delay, 
                    exclude=HFI_info['HFI']['exclude'])
                # collates data
                data['HFI'] = arrangement
            
            # get cholinergic modulation class
            modulation = modulate.set_ACh(cell, mod_factors, [ACh_t], target_x=target_x, play=mech_scale, dt=h.dt)
            
            # run simulation
            h.finitialize(-80)
            while h.t < clus_params['stop_t']+HFI_delay+(clus_params['stim_n']*clus_params['isi']):
                h.fadvance()
            tm = tm.to_python()
            vm = vm.to_python()
            
            
            # collate data
            data[clus_lab][ACh_lab] = {'vm':vm}
            data['meta'] = {'id':int(cell_index), 'round':run_info['round'], 'cell_type':model_data['cell_type'], 'tm':tm, 
                            'rheo':rheobase}
            
            
            # calculate dpp duration and amplitude
            if dur_and_amp:
                base_t_start = tm.index(min(tm, key=lambda x:abs(x-(clus_params['stim_t']+clus_params['pre_t']))))
                base_t_end = tm.index(min(tm, key=lambda x:abs(x-clus_params['stim_t'])))
                base_vm = np.mean(vm[base_t_start:base_t_end])
                data[clus_lab][ACh_lab]['dur'] = cf.dpp_dur(tm, vm, base_vm, clus_params['stim_t'])
                data[clus_lab][ACh_lab]['amp'] = cf.dpp_amp(tm, vm, base_vm, clus_params['stim_t'])
                
            # get spike-related data
            if spike:
                thresh = 0
                # checks whether spike occured
                data[clus_lab][ACh_lab]['spiked'] = 0
                if max(vm) > thresh:
                    data[clus_lab][ACh_lab]['spiked'] = 1
                
            
        
    return data





def dpp_DA_modded(model_data,
                   cell_index,
                   run_info,
                   noise = True,
                   HFI = False,
                   HFI_delay = 0,
                   dur_and_amp = True,
                   spike = False,
                   mod_tar = 'all'):
    '''
    Provides clustered glutamatergic input to SPNs to generate the dendritic 
    plateau potential in the absence of modulation.
    
    INPUT(S):
        - model_data: model paramaters (specification, cell type, model 
            sets, and stimulation targets) [dict]
        - stim_data: stimulation parameters (number of inputs, time of inputs,
            input inter-spike intervals, time to stop simulating) [dict]
        - cell_index: cell specification being simulated [int]
        - run_info: information about the simulation run (this simulation and 
            total number of simulations) [dict]
        - dur_and_amp: whether to calculate the duration and peak amplitude of
            the plateau potential (default True) [bool]
        
    OUTPUT(S):
        - data: simulated data including: simulation times; simulated voltages;
            distance of stimulated target to the soma; rheobase of the model; 
            half-max width of the potential; peak amplitude of the potential;
            cell specification being simulated; cell type being simulated
            [dict]
        - ncon: NetCon object
            
    Thomas Binns (author), 13/02/21
    '''
    
    # ===== print simulation info to monitor progress =====
    print('Simulating cell specification {} of {}'.format( \
          run_info['curr_n']+1,run_info['tot_n']),flush = True)
    
    
    # ===== gets stimulation info =====
    # clustered inputs
    clus_info = cf.params_for_input(model_data['cell_type'], 'clustered')
    clus_params = clus_info['clustered']['params']
    
    
    # modulation inputs
    mod_factors = cf.draw_factors_DA(model_data['cell_type'], mode='random')
    DA_info = cf.params_for_input(model_data['cell_type'], 'DA')
    DA_params = DA_info['DA']['params']
    
    # gets vector for modulation timing
    state = [1 if ht >= DA_params['stim_t'] and ht <= DA_params['stop_t'] \
             else 0 for ht in np.arange(0,clus_params['stop_t'],h.dt)]
    state = h.Vector(state)
    
    # collates time-dependent mechanism modulation scaling
    mech_scale = {}
    for key in mod_factors:
        mech_scale[key] = state
    
    
    # noise inputs
    if noise:
        noise_info = cf.params_for_input(model_data['cell_type'], 'noise')
        noise_params = noise_info['noise']['params']
    
    # HFIs
    if HFI:
        HFI_info = cf.params_for_input(model_data['cell_type'], 'HFI')
        HFI_params = HFI_info['HFI']['params']
        
    # ===== simulation =====
    data = {}
    
    for i, clus_tar in enumerate(clus_info['clustered']['target']): # for each simulation target
        
        clus_lab = clus_info['clustered']['label'][i] # label for input target
        data[clus_lab] = {}
        
        # get the targets for cholinergic input
        DA_targets = {}
        DA_targets['target'] = [clus_tar]
        DA_targets['label'] = [clus_lab]
        for j, DA_t in enumerate(DA_info['DA']['target']):
            DA_targets['target'].append(DA_t)
            DA_targets['label'].append(DA_info['DA']['label'][j])
        target_x = .5
        
        if mod_tar == 'all':
            DA_targets['target'] = ['all']
            DA_targets['label'] = ['all']
            target_x = 'all'
        
        for j, DA_t in enumerate(DA_targets['target']): # for each cholinergic input target
            
            DA_lab = DA_targets['label'][j]
            
            # initiate cell
            cell = build.MSN(params=model_data['specs']['par'],
                             morphology=model_data['specs']['morph'],
                             variables=model_data['model_sets'][cell_index]['variables'])
            rheobase = model_data['model_sets'][cell_index]['rheobase']
            
            
            # record vectors
            tm = h.Vector()
            tm.record(h._ref_t)
            vm = h.Vector()
            vm.record(cell.soma(0.5)._ref_v)
            
            
            # add clustered inputs
            clus_syn, clus_stim, clus_ncon, clus_d2soma = cf.set_clustered_stim(cell, clus_tar ,n=clus_params['stim_n'],
                act_time=clus_params['stim_t'], ISI=clus_params['isi'])
            
            # add background noise
            if noise:
                noise_syn, noise_stim, noise_ncon = cf.set_bg_noise(cell, model_data['cell_type'], 
                    fglut=noise_params['freq_glut'], fgaba=noise_params['freq_gaba'], dendOnly=noise_params['only dend'])
    
            
            # add high-frequency inputs
            if HFI:
                # adds HFI
                HFI_syn, HFI_stim, HFI_ncon, arrangement = cf.set_HFI(cell, model_data['cell_type'], 
                    freq=HFI_params['freq'], n_inputs=HFI_params['n_inputs'],
                    delay=clus_params['stim_t']+(clus_params['stim_n']*clus_params['isi'])+HFI_delay, 
                    exclude=HFI_info['HFI']['exclude'])
                # collates data
                data['HFI'] = arrangement
            
            # get cholinergic modulation class
            modulation = modulate.set_DA(cell, mod_factors, [DA_t], target_x=target_x, play=mech_scale, dt=h.dt)
            
            # run simulation
            h.finitialize(-80)
            while h.t < clus_params['stop_t']+HFI_delay+(clus_params['stim_n']*clus_params['isi']):
                h.fadvance()
            tm = tm.to_python()
            vm = vm.to_python()
            
            
            # collate data
            data[clus_lab][DA_lab] = {'vm':vm}
            data['meta'] = {'id':int(cell_index), 'round':run_info['round'], 'cell_type':model_data['cell_type'], 'tm':tm, 
                            'rheo':rheobase}
            
            
            # calculate dpp duration and amplitude
            if dur_and_amp:
                base_t_start = tm.index(min(tm, key=lambda x:abs(x-(clus_params['stim_t']+clus_params['pre_t']))))
                base_t_end = tm.index(min(tm, key=lambda x:abs(x-clus_params['stim_t'])))
                base_vm = np.mean(vm[base_t_start:base_t_end])
                data[clus_lab][DA_lab]['dur'] = cf.dpp_dur(tm, vm, base_vm, clus_params['stim_t'])
                data[clus_lab][DA_lab]['amp'] = cf.dpp_amp(tm, vm, base_vm, clus_params['stim_t'])
                
            # get spike-related data
            if spike:
                thresh = 0
                # checks whether spike occured
                data[clus_lab][DA_lab]['spiked'] = 0
                if max(vm) > thresh:
                    data[clus_lab][DA_lab]['spiked'] = 1
                    
                
            
        
    return data





def ACh_modulation(model_data,
                   stim_data,
                   cell_index,
                   run_info,
                   mod_factors,
                   dur_and_amp = True):
    '''
    Provides clustered glutamatergic input to SPNs to generate the dendritic 
    plateau potential alongside cholinergic modulation.
    
    INPUT(S):
        - model_data: model paramaters (specification, cell type, model 
            sets, and stimulation targets) [dict]
        - stim_data: stimulation parameters (number of inputs, time of inputs,
            input inter-spike intervals, time to stop simulating) [dict]
        - cell_index: cell specification being simulated [int]
        - run_info: information about the simulation run (this simulation and 
            total number of simulations) [dict]
        - mod_factors: mechanism:modulation value pairs [dict]
        - dur_and_amp: whether to calculate the duration and peak amplitude of
            the plateau potential (default True) [bool]
        
    OUTPUT(S):
        - data: simulated data including: simulation times; simulated voltages;
            distance of stimulated target to the soma; rheobase of the model; 
            half-max width of the potential; peak amplitude of the potential;
            cell specification being simulated; cell type being simulated
            [dict]
        - ncon: NetCon object
    
    Thomas Binns (author), 30/01/21
    '''
    
    # ===== print simulation info to monitor progress =====
    print('Simulating cell specification {} of {}'.format( \
          run_info['curr_n']+1,run_info['tot_n']),flush = True)
    
    
    # ===== simulation =====
    
    clus_params = stim_data['clustered']['params'] # parameters for clustered input
    ACh_params = stim_data['ACh']['params'] # parameters for cholinergic input
    
    # gets vector for modulation timing
    state = [1 if ht >= ACh_params['stim_t'] and ht <= ACh_params['stop_t'] \
             else 0 for ht in np.arange(0,clus_params['stop_t'],h.dt)]
    if 'kaf' in mod_factors:
        kaf_state = [x * mod_factors['kaf'] for x in state]
        kaf_state = h.Vector(kaf_state)
    state = h.Vector(state)
    
    # collates time-dependent mechanism modulation scaling
    mech_scale = {}
    for key in mod_factors:
        if key == 'kaf':
            mech_scale[key] = kaf_state
        else:
            mech_scale[key] = state
    
    data = {}
    
    
    for i, clus_t in enumerate(stim_data['clustered']['target']): # for each clustered input target
        
        clus_lab = stim_data['clustered']['label'][i] # label for clustered input target
        
        data[clus_lab] = {}
        
        # get the targets for cholinergic input
        ACh_targets = {}
        ACh_targets['target'] = [clus_t]
        ACh_targets['label'] = [clus_lab]
        for j, ACh_t in enumerate(stim_data['ACh']['target']):
            ACh_targets['target'].append(ACh_t)
            ACh_targets['label'].append(stim_data['ACh']['label'][j])
        
        for j, ACh_t in enumerate(ACh_targets['target']): # for each cholinergic input target
            
            ACh_lab = ACh_targets['label'][j]
            
            # initiate cell
            cell = build.MSN(params=model_data['specs']['par'],
                             morphology=model_data['specs']['morph'],
                             variables=model_data['model_sets'][cell_index]['variables'])
            rheobase = model_data['model_sets'][cell_index]['rheobase']
        
            # gets distance info for targeted regions
            dists = cf.get_dists(cell,only_sec=ACh_targets['target'])
            
            # record vectors
            tm = h.Vector()
            tm.record(h._ref_t)
            vm = h.Vector()
            vm.record(cell.soma(0.5)._ref_v)            
            
            # add clustered inputs
            syn,stim,ncon,d2soma = cf.set_clustered_stim(cell, clus_t, \
                n = clus_params['stim_n'], act_time = clus_params['stim_t'], \
                ISI = clus_params['isi'])
            
            # get cholinergic modulation class
            mod = modulate.set_ACh(cell, mod_factors, [ACh_t], target_x=.5,
                play=mech_scale, dt=h.dt)
            
            # run simulation
            h.finitialize(-80)
            while h.t < clus_params['stop_t']:
                h.fadvance()
            tm = tm.to_python()
            vm = vm.to_python()
            
            # collate data
            data[clus_lab][ACh_lab] = \
                {'tm':tm, 'vm':vm, 'clust_dist':dists[clus_t], \
                 'ACh_dist':dists[ACh_t], 'rheo':rheobase,  \
                 'id':int(cell_index), 'cell_type':model_data['cell_type']}
            
            if dur_and_amp:
                # calculate dpp duration and amplitude
                base_t = tm.index(min(tm, key=lambda x:abs(x-clus_params['stim_t'])))-1
                data[clus_lab][ACh_lab]['dur'] = cf.dpp_dur(tm,vm,vm[base_t],clus_params['stim_t'])
                data[clus_lab][ACh_lab]['amp'] = cf.dpp_amp(tm,vm,vm[base_t],clus_params['stim_t'])
            
        
    return data





