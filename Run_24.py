import Helper_Functions as hf

# Set up the constant dict
constant_dict ={
    'At':4.16, # mM
    'K_app':4.4 * 10**(-6),  #
    'asp':1.6, # mM
    'glu':5.3, # mM
    'Pi':2.44*10**(-3), # mM
    'Nt':1.07, # mM
    'C':6.75 * 10**(-6), #mMmV^-1
    'a':0.1, # mV^-1
    'b':0.004, # mM^-1
    'psi_m':150, # mV
    'K':2, #mM
    'F':96485*10**(-3) , # C mMol^-1
    'R':8.314 , # mMolK ^-1
    'T':298, # K
    'K_6_eq':0.12,
}

# Set up the conentration dict with the steady state concentrations

# Use the inverted steady state when running inv_factor = -1
#Inverted Steady State
# concentration_dict ={
#     'Pyr_bar': 0.13956603,
#     'NAD_bar': 0.82738712,
#     'OAA_bar': 0.0050072,
#     'AcCoA_bar': 0.06980711,
#     'Cit_bar': 0.39991786,
#     'KG_bar': 0.24883437,
#     'ATP_bar': 2.50625812,
#     'psi_bar': 149.82949285,
# }

# Use the steady state when running inv_factor = 1
# Steady State
concentration_dict ={
    'Pyr_bar': 0.13956603,
    'NAD_bar': 0.82738712,
    'OAA_bar': 0.0050072,
    'AcCoA_bar': 0.06980711,
    'Cit_bar': 0.39991786,
    'KG_bar': 0.24883437,
    'ATP_bar': 2.50625812,
    'psi_bar': 149.82949285,
}

# Calculate the starting concentration values, which are +/- 10% from the steady state concentration
concentration_dict['Pyr'] = concentration_dict['Pyr_bar']*1.1
concentration_dict['NAD'] = concentration_dict['NAD_bar']*0.9
concentration_dict['OAA'] = concentration_dict['OAA_bar']*.9
concentration_dict['AcCoA'] = concentration_dict['AcCoA_bar']*.9
concentration_dict['Cit'] = concentration_dict['Cit_bar']*1.1
concentration_dict['KG'] = concentration_dict['KG_bar']*.9
concentration_dict['ATP'] = concentration_dict['ATP_bar']*1.1
concentration_dict['psi'] = concentration_dict['psi_bar']*.9

# Set the rate dict
rate_dict ={
    'k1':38, # mMs^-1
    'k2':152, #
    'k3':57142, #
    'k4':53, #
    'k5':82361*10**(-3), # mM^-2 s^-1
    'k6':3.2, #
    'k7':40, # mMS ^-1
    'k8':3.6*10**3, # s^-1
    'k_resp':2.5*10**3, # 2.5 mMs^-1
    'k_atp':131.9*10**3, # mMs^-1
    'k_ant':0.05387 *10**3, #mMs^-1
    'k_leak':0.426 #mM(mVs)^-1
}

# Set to -1 to run under inverted assumption
inv_factor = 1

# Modify rate constants based off predictions from the machine learning model results at 24 HPI
rate_dict['k2'] = rate_dict['k2'] + rate_dict['k2']*-0.13251471026541703*inv_factor # PD Complex
rate_dict['k3'] = rate_dict['k3'] + rate_dict['k3']*0.5514591191595479*inv_factor # Citrate Synthase
rate_dict['k4'] = rate_dict['k4'] + rate_dict['k4']*((-0.1310524382650617+0.5969921753499483)/2) *inv_factor# ID complex + Aconitase
rate_dict['k5'] = rate_dict['k5'] + rate_dict['k5']*-0.021764402776499214*inv_factor # AK Complex
rate_dict['k6'] = rate_dict['k6'] + rate_dict['k6']*((-0.5645056313414636+-0.878429548494093+0.5966979833763801+-0.42810565855100235)/4)*inv_factor # Succinyl-CoA synthase, Succinic dehydrogenase, Fumarase, Malate dehydrogenase

# Create the TCA model, Run it, and then save the concentration matrix and tau list
tca = hf.TCA_Cycle(constant_dict, concentration_dict, rate_dict)
tca.run(0.0000001, 40, method='RK4')
tca.save('./Classes/APC523/Run_24_40_V2_steady')