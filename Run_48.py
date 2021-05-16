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
#     'Pyr_bar': 0.13113873,
#     'NAD_bar': 0.80521847,
#     'OAA_bar': 0.00475478,
#     'AcCoA_bar': 0.17232445,
#     'Cit_bar': 0.56493739,
#     'KG_bar': 0.27041412,
#     'ATP_bar': 2.561009,
#     'psi_bar': 152.60864354,
# }

# Use the steady state when running inv_factor = 1
#Steady State
concentration_dict ={
    'Pyr_bar': 0.15084048,
    'NAD_bar': 0.84400219,
    'OAA_bar': 0.00525139,
    'AcCoA_bar': 0.04102625,
    'Cit_bar': 0.3037923,
    'KG_bar': 0.22458203,
    'ATP_bar': 2.43072101,
    'psi_bar': 147.13135702,
}

#Calculate the starting concentration values, which are +/- 10% from the steady state concentration
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

# Modify rate constants based off predictions from the machine learning model results at 48 HPI
rate_dict['k2'] = rate_dict['k2'] + rate_dict['k2']*0.118632524021445*inv_factor # PD Complex
rate_dict['k3'] = rate_dict['k3'] + rate_dict['k3']*0.5214117541190887*inv_factor # Citrate Synthase
rate_dict['k4'] = rate_dict['k4'] + rate_dict['k4']*((-0.1485026292034054+0.9442442393155126)/2)*inv_factor # ID complex + Aconitase
rate_dict['k5'] = rate_dict['k5'] + rate_dict['k5']*-0.03375129144961664*inv_factor # AK Complex
rate_dict['k6'] = rate_dict['k6'] + rate_dict['k6']*((-0.7624722973789384+1.9236571856095588+-0.45380219458984217+0.4717007310085516)/4)*inv_factor # Succinyl-CoA synthetase, Succinic dehydrogenase, Fumarase, Malate dehydrogenase


# Create the TCA model, Run it, and then save the concentration matrix and tau list
tca = hf.TCA_Cycle(constant_dict, concentration_dict, rate_dict)
tca.run(0.0000001, 40, method='RK4')
tca.save('./Classes/APC523/Run_48_40_V2_steady')