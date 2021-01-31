

#import numpy as np
from neuron import h
import common_functions as cf
    

class DA():
    ''' 
    class for setting DA modulation of cell. 
    arguments:
        cell = cell object to modulate
        mod_dict = dict with chan:factor pairs
    additional arguments carries default behavior of modulation.
    
    The 'naf', 'kas', 'kir', 'cal12','cal13','can' and 'car channels are reported to be modulated by DA in dspn. Aditionally the car channel is modulated in ispn.
    
    By default the AIS is not modulated, since this was found to increase chance of positive modulation in dspn.
        ---OBS. if uniform modulation is wanted this has to be stated!!!!
    
    Modulation has in some channel(s?) been reported to be a shift in gating. This is not
        implemented for DA (but see ACh class for example).
    
    Modulation can be set individually for gaba, glut and intrinsic.
        by default all are used.
    
    The _reset_mod method can be used to turn of modulation.
    '''
    
    def __init__(self,  cell, mod_dict,                                    
                        modulation='noAxon', 
                        intrinsic_mod=1,
                        gaba_mod=1,
                        glut_mod=1,
                        play=[],
                        syn_dict={},
                        dt=0.025
                        ):
        
        self.cell = cell
        self.mod_dict = mod_dict
        self.syn_dict = syn_dict
        self.play = play
        self.dt = dt
        
        if modulation == 'uniform':
            self.compartments = [cell.dendlist, cell.somalist, cell.axonlist]
        elif modulation == 'noAxon':
            self.compartments = [cell.dendlist, cell.somalist]
        else: raise Exception('Error, "{}" modulation, not permitted\n\tuse "uniform" or "noAxon"'.format(modulation))
        
        
        self._set_modulation(   intrinsic_mod,
                                gaba_mod,
                                glut_mod
                                )
        
    def _set_modulation(self,
                        intrinsic_mod,
                        gaba_mod,
                        glut_mod):
        for comp in self.compartments:
            for sec in comp:
                for seg in sec:
                    if intrinsic_mod:   self._update_conductance(seg)
                    if gaba_mod:        self._set_gaba(seg)
                    if glut_mod:        self._set_glut(seg)
    
    
    def _update_conductance(self, seg, reset=0):
        for mech in seg:
            if mech.name() in self.mod_dict:
                # set shift
                mech.damod = 1
                mech.maxMod = self.mod_dict[mech.name()]
                if reset: 
                    mech.level = 0
                elif len(self.play) and mech.name() in self.play:
                    self.play[mech.name()].play(mech._ref_level, self.dt)
                else:
                    mech.level = 1
                            
    def _set_gaba(self, seg, reset=0):
        for syn in seg.point_processes():
            if 'gaba' in syn.hname():
                syn.damod = 1
                syn.maxMod = self.syn_dict['GABA']
                if reset:
                    syn.level = 0 
                elif len(self.play) and 'gaba' in self.play:
                    self.play['gaba'].play(syn._ref_level, self.dt)
                else:
                    syn.level = 1
    
    def _set_glut(self, seg, reset=0):
        for syn in seg.point_processes():
            if 'glut' in syn.hname():
                syn.damod = 1
                syn.maxModNMDA = self.syn_dict['NMDA']
                syn.maxModAMPA = self.syn_dict['AMPA']
                if reset:
                    syn.l1AMPA = 0; syn.l1NMDA = 0
                elif len(self.play) and 'glut' in self.play:
                    self.play['glut'].play(syn._ref_l1NMDA, self.dt)
                    self.play['glut'].play(syn._ref_l1AMPA, self.dt)
                else:
                    syn.l1AMPA = 1; syn.l1NMDA = 1
    
    def _reset_mod(self):
        if len(self.play):
            for trans in self.play.values():
                trans.play_remove()
        for comp in self.compartments:
            for sec in comp:
                for seg in sec:
                    self._update_conductance(seg, reset=1)
                    self._set_glut(seg, reset=1)
                    self._set_gaba(seg, reset=1)
      

class ACh():
    ''' 
    class for setting ACh modulation of cell. 
    arguments:
        cell = cell object to modulate
        mod_dict = dict with chan:factor pairs
    additional arguments carries default behavior of modulation.
    
    The 'naf', 'kir','cal12','cal13','can' and 'im' channels are reported to be modulated by ACh.
    Additionally "kaf" is (by default) shifted 20 mv more negative
    
    Modulation can be set individually for gaba, glut, intrinsic and kaf shift.
        by default all are used.
    '''
    
    def __init__(self,  cell, mod_dict,                                    
                        modulation='uniform', 
                        intrinsic_mod=1,
                        shift_kaf=1,
                        gaba_mod=1,
                        glut_mod=1,
                        mv_shift_kaf=20,
                        syn_dict={},
                        play=[],
                        dt=0.025
                        ):
        
        self.cell = cell
        self.mod_dict = mod_dict
        self.mv_shift_kaf = mv_shift_kaf
        self.syn_dict = syn_dict
        self.play = play
        self.dt = dt
        
        if modulation == 'uniform':
            self.compartments = [cell.dendlist, cell.somalist, cell.axonlist]
        elif modulation == 'noAxon':
            self.compartments = [cell.dendlist, cell.somalist]
        else: raise Exception('Error, "{}" modulation, not permitted\n\tuse "uniform" or "noAxon"'.format(modulation))
        
        self._set_modulation(   intrinsic_mod,
                                shift_kaf,
                                gaba_mod,
                                glut_mod
                                )
        
    def _set_modulation(self,
                        intrinsic_mod,
                        shift_kaf,
                        gaba_mod,
                        glut_mod):
        for comp in self.compartments:
            for sec in comp:
                for seg in sec:
                    if intrinsic_mod:   self._update_conductance(seg)
                    if shift_kaf:       self._shift_kaf(seg)
                    if gaba_mod:        self._set_gaba(seg)
                    if glut_mod:        self._set_glut(seg)
        
    def _update_conductance(self, seg, reset=0):
        for mech in seg:
            if mech.name() in self.mod_dict:
                # set shift
                mech.damod = 1
                mech.max2 = self.mod_dict[mech.name()]
                if reset:
                    mech.lev2 = 0
                if len(self.play) and mech.name() in self.play:
                    self.play[mech.name()].play(mech._ref_lev2, self.dt) # ._ref_
                else:
                    mech.lev2 = 1
    
    def _shift_kaf(self, seg, reset=0):
        for mech in seg:
            if mech.name() == 'kaf':
                # set shift
                if reset:
                    mech.modShift = 0
                elif len(self.play) and 'kaf' in self.play:
                    self.play['kaf'].play(mech._ref_modShift, self.dt) # ._ref_
                else:
                    mech.modShift = self.mv_shift_kaf
                            
    def _set_gaba(self, seg, reset=0):
        for syn in seg.point_processes():
            if 'gaba' in syn.hname():
                syn.damod = 1
                syn.max2 = self.syn_dict['GABA']
                if reset:
                    syn.lev2 = 0
                elif len(self.play) and 'gaba' in self.play:
                    self.play['gaba'].play(syn._ref_lev2, self.dt) # ._ref_
                else:
                    syn.lev2 = 1
    
    def _set_glut(self, seg, reset=0):
        for syn in seg.point_processes():
            if 'glut' in syn.hname():
                syn.damod = 1
                syn.max2NMDA = self.syn_dict['NMDA']
                syn.max2AMPA = self.syn_dict['AMPA']
                if reset:
                    syn.l2AMPA = 0
                    syn.l2NMDA = 0
                elif len(self.play) and 'glut' in self.play:
                    self.play['glut'].play(syn._ref_l2NMDA, self.dt) # ._ref_
                    self.play['glut'].play(syn._ref_l2AMPA, self.dt) # ._ref_
                else:
                    syn.l2AMPA = 1
                    syn.l2NMDA = 1
                            
    
    def _reset_mod(self):
        if len(self.play):
            for trans in self.play.values():
                trans.play_remove()
        for comp in self.compartments:
            for sec in comp:
                for seg in sec:
                    self._shift_kaf(seg, reset=1)
                    self._update_conductance(seg, reset=1)
                    self._set_glut(seg, reset=1)
                    self._set_gaba(seg, reset=1)
                    
                    
                    
                    
class set_ACh():
    ''' 
    Class for setting cholinergic modulation of the cell. The 'naf', 'kir',
        'cal12','cal13','can' and 'im' channels are reported to be modulated by 
        ACh. Additionally "kaf" is shifted to a more negative voltage.
    '''
    
    def __init__(self,  cell,
                        mod_dict,                                    
                        target = ['all'],
                        target_x = 'all',
                        play=[],
                        shift_kaf = 5,
                        dt=0.025):
        '''
        Class initialisation and modulation of cell.
        
        INPUT(S):
            - cell: cell to modulate [MSN object]
            - mod_dict: mechanism:modulation value pairs [dict]
            - target: section(s) of cell to modulate (default all sections) 
                [list of section names]
            - target_x: segment of section(s) to modulate (default: all 
                segments). Can specify a specific segment as a number (e.g. 0.5
                to modulate the middle segment) [str or number]
            - play: key:value pairs containing time-specific scaling of each 
                mechanism's modulation (default: 1*scaling of modulation 
                throughout the simulation). Specify time-dependent modulation 
                as a hoc vector in which the value of vector[t/dt] == scaling 
                of modulation at time t in the simulation (see h.Vector.play 
                documentation for more information) [dict of h.Vector(s)]
            - shift_kaf: voltage by which the mechanism 'kaf' is shifted * -1
                (i.e. default: shift_kaf = 10 -> kaf shifted -10 mV). If 'play'
                is given and contains a 'kaf' key, 'shift_kaf' is ignored in
                favour of that specified by 'play' [number]
            ==== N.B. the 'kaf' scaling in 'play' (if given) should follow 
                    this same logic (i.e. for a voltage shift from 0 mV at time 
                    t to -10 mV at time t+1, play['kaf'][t/dt] == 0 and 
                    play['kaf'][t+1/dt] == 10)
            - dt: size of the time step in the simulation used for determining
                when to apply modulation scaling associated with 'play'
                (default: 0.025; units of ms) [number]
        
        OUTPUT(S):
            None
        
        Thomas Binns (modified), 30/01/21
        '''
        
        # checks for appropriate inputs
        if isinstance(target_x,str):
            if target_x != 'all':
                raise ValueError("The requested segment target '{}' is not recognised.\nThis should be 'all' or a number, e.g. 0.5.".format(target_x))
        
        # assigns values to class
        self.cell = cell
        self.mod_dict = mod_dict
        self.target = target
        self.target_x = target_x
        self.play = play
        self.shift_kaf = shift_kaf
        self.dt = dt
        
        # gets targets for modulation
        if 'all' in self.target:
            self.target = self.cell.allseclist
        
        # performs modulation of cell
        self._set_modulation()
        
        
    def _set_modulation(self):
        '''
        Modulates the requested segments of the requested sections.
        '''
        for tar in self.target:
            for sec in self.cell.allseclist:
                if sec.name() == tar:
                    if self.target_x == 'all':
                        for seg in sec:
                            self._mod_mech(seg)
                            self._mod_chan(seg)
                    else:
                        self._mod_mech(sec(self.target_x))
                        self._mod_chan(sec(self.target_x))
                
                
    def _mod_mech(self, seg, reset=False): # shift conductance
        '''
        Modulates kaf and other mechanisms (but not GABA, NMDA, or AMPA).
        '''
        for mech in seg:
            if mech.name() in self.mod_dict:
                if mech.name() == 'kaf': # shift kaf
                    if reset:
                        mech.modShift = 0
                    elif len(self.play) and 'kaf' in self.play:
                        self.play['kaf'].play(mech._ref_modShift, self.dt)
                    else:
                        mech.modShift = self.shift_kaf
                        
                else: # apply other mechanism conductance changes
                    mech.damod = 1
                    mech.max2 = self.mod_dict[mech.name()]
                    if reset:
                        mech.lev2 = 0
                    if len(self.play) and mech.name() in self.play:
                        self.play[mech.name()].play(mech._ref_lev2, self.dt)
                    else:
                        mech.lev2 = 1
    
    
    def _mod_chan(self, seg, reset=False):
        '''
        Modulates GABA, NMDA, and AMPA mechanisms.
        '''
        for syn in seg.point_processes(): # modulate channels if applicable
            
            if 'gaba' in syn.hname() and 'GABA' in self.mod_dict:
                syn.damod = 1
                syn.max2 = self.mod_dict['GABA']
                if reset:
                    syn.lev2 = 0
                elif len(self.play) and 'GABA' in self.play:
                    self.play['GABA'].play(syn._ref_lev2, self.dt)
                else:
                    syn.lev2 = 1
            
            elif 'glut' in syn.hname()and 'GLUT' in self.mod_dict:
                syn.damod = 1
                syn.max2NMDA = self.mod_dict['NMDA']
                syn.max2AMPA = self.mod_dict['AMPA']
                if reset:
                    syn.l2AMPA = 0
                    syn.l2NMDA = 0
                elif len(self.play) and 'NMDA' in self.play() and 'AMPA' in self.play:
                    self.play['NMDA'].play(syn._ref_l2NMDA, self.dt)
                    self.play['AMPA'].play(syn._ref_l2AMPA, self.dt)
                else:
                    syn.l2AMPA = 1
                    syn.l2NMDA = 1
                            
    
    def _reset_mod(self):
        '''
        Reverses modulation to return cell status to normal.
        '''
        if len(self.play):
            for trans in self.play.values():
                trans.play_remove()
                
        for tar in self.target:
            for sec in self.cell.allseclist:
                if sec.name() == tar:
                    self._mod_mech(sec(self.target_x), reset=True)
                    self._mod_chan(sec(self.target_x), reset=True)
    
