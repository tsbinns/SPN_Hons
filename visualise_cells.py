'''
Shows cell morphology, with regions of input highlighted
    - So far only visualisation of clustered inputs implemented
    
Thomas Binns (author), 29/01/21
'''



import neuron                                           as nrn
import MSN_builder                                      as build
import pickle
import common_functions                                 as cf
import numpy                                            as np
from   neuron            import h
from   matplotlib        import pyplot                  as plt
from   matplotlib.colors import LinearSegmentedColormap as lsc



# what to visualise =============

cell_type = input("What cell type should be visualised: dspn; or ispn?")
input_type = input("What input should be visualised: clustered; clustered+stim; or clustered+Ach?")

vis_info = cf.params_for_input(cell_type,input_type)


# sets up cell =============

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

with open(specs[cell_type]['lib'], 'rb') as f:
    model_sets = pickle.load(f, encoding="latin1")

# build cell
cell = build.MSN(params=specs[cell_type]['par'],
                 morphology=specs[cell_type]['morph'],
                 variables=model_sets[0]['variables'])


# changes voltage in targets (above baseline -65 mV) to highlight regions

# increases voltage of clustered inputs a large amount to highlight these 
# regions in red
for i, tar in enumerate(vis_info['clustered']['target']):
    if tar[:4] == 'dend':
        h.dend[int(tar[5:-1])](.5).v = 35
    else:
        raise ValueError("Trying to visualise input to '{}', but only visualisation of input to 'dend' currently supported".format(tar[:4]))

# increases voltage of stimulation inputs an intermediate amount to highlight 
# these regions in green
if input_type == 'clustered+stim':
    for i, tar in enumerate(vis_info['stim']['target']):
        if tar[:4] == 'dend':
            h.dend[int(tar[5:-1])](.5).v = -2.5
        else:
            raise ValueError("Trying to visualise input to '{}', but only visualisation of input to 'dend' currently supported".format(tar[:4]))

# increases voltage of stimulation inputs a small amount to highlight  these 
# regions in blue
if input_type == 'clustered+Ach':
    for i, tar in enumerate(vis_info['Ach']['target']):
        if tar[:4] == 'dend':
            h.dend[int(tar[5:-1])](.5).v = -27.5
        else:
            raise ValueError("Trying to visualise input to '{}', but only visualisation of input to 'dend' currently supported".format(tar[:4]))


# visualisation ==========

# creates colourmap
vis_cmap = lsc.from_list('my_cmap',['black','blue','green','red'])

# plots dendrites
sections = h.SectionList([sec for sec in h.allsec() if 'dend' in str(sec)])
ps = h.PlotShape(sections,False).plot(plt)
ps._do_plot(0,1,h.dend,'v',cmap=vis_cmap)

# Make sphere to mimic soma
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 12 * np.outer(np.cos(u), np.sin(v))
y = 12 * np.outer(np.sin(u), np.sin(v))
z = 12 * np.outer(np.ones(np.size(u)), np.cos(v))
ps.plot_surface(x, y+5, z+10, color='black')
ps.plot_surface(x, y+5, z+14, color='black')
ps.plot_surface(x, y+5, z+17, color='black')

# removes plot background
ps.xaxis.pane.fill = False
ps.yaxis.pane.fill = False
ps.zaxis.pane.fill = False

# removes panes
ps.xaxis.pane.set_edgecolor('w')
ps.yaxis.pane.set_edgecolor('w')
ps.zaxis.pane.set_edgecolor('w')

# removes grid
ps.grid(False)

# removes axes
ps.set_axis_off()

# adds title
ps.text2D(.4, .9, '{}, {} inputs'.format(cell_type,input_type), \
          transform=ps.transAxes)
