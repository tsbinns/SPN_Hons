

import common_functions                 as cf
import scipy.stats                      as stats
from   statsmodels.stats.diagnostic import lilliefors
import matplotlib.pyplot                as plt





dur_and_amp = 1



if dur_and_amp == 1:
    
    # load data
    data = cf.load_data('Data/dspn_HFI[0]+0_validation.json')
    clus_info = data['meta']['clustered']
    target_labels = clus_info['label']
    cell_type = data['meta']['cell type']
    
    # normality checking =====
    
    # plot histograms
    fig, axs = plt.subplots(2,2)
    fig.suptitle(cell_type)
    for i, lab in enumerate(target_labels):
        
        # duration data
        axs[i,0].hist(data['all'][lab]['dur'])
        axs[i,0].set_title('dur, ' + lab)
        
        # amplitude data
        axs[i,1].hist(data['all'][lab]['amp'])
        axs[i,1].set_title('amp, ' + lab)
        
        plt.tight_layout()
        
    
    # gets stats
    data['stats'] = {}
    for lab in target_labels:
        data['stats'][lab] = {'dur':{}, 'amp':{}}
        # kurtosis
        data['stats'][lab]['dur']['kurt'] = stats.kurtosistest(data['all'][lab]['dur'])
        data['stats'][lab]['amp']['kurt'] = stats.kurtosistest(data['all'][lab]['amp'])
        # skewness
        data['stats'][lab]['dur']['skew'] = stats.skewtest(data['all'][lab]['dur'])
        data['stats'][lab]['amp']['skew'] = stats.skewtest(data['all'][lab]['amp'])
        # ks test
        data['stats'][lab]['dur']['ks'] = lilliefors(data['all'][lab]['dur'])
        data['stats'][lab]['amp']['ks'] = lilliefors(data['all'][lab]['amp'])
    
    
    # tests on data =====
    '''
    for dspn, dur and amp data was normally distributed; for ispn, a bit iffy but still within reason for normal distribution
    '''
    
    data['stats']['diff'] = {'dur':{}, 'amp':{}}
    
    # duration
    data['stats']['diff']['dur']['result'] = stats.ttest_rel(data['all'][target_labels[1]]['dur'], \
        data['all'][target_labels[0]]['dur'])
    data['stats']['diff']['dur']['label'] = 'paired t-test, two-sided'
    # ttest_rel() doesn't have one-sided implemented, so to see if duration of one greater than other, 
        # t-statistic must be > 0 and p-value/2 < .05
        
    # amplitude
    data['stats']['diff']['amp']['result'] = stats.ttest_rel(data['all'][target_labels[1]]['amp'], \
        data['all'][target_labels[0]]['amp'])
    data['stats']['diff']['amp']['label'] = 'paired t-test, two-sided'
    
    
    
    