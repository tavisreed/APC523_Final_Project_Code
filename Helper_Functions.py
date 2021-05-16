import pickle
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv


# This Function is used to save any object as a pickle file
def save_object(obj, save_name):
    f = open(save_name + ".pkl", "wb")
    pickle.dump(obj, f, -1)
    f.close()

# This Function is used to load pickle files
def load_object(file_name):
    f = open(file_name + ".pkl", "rb")
    obj = pickle.load(f)
    f.close()
    return obj


# This function is to plot the concentration matrix obtained after running the TCA model
# Provided as seperate function for when concentration matrix is loaded from pickle file
def plot_concentration_matrix(matrix, t, title=''):
    # Set the figure Size
    f, ax = plt.subplots(figsize=(10,7))

    # Extract the row corresponding to each of the Metabolites
    Pyr = matrix[0]
    NAD = matrix[1]
    OAA = matrix[2]
    AcCoA = matrix[3]
    Cit = matrix[4]
    KG = matrix[5]
    ATP = matrix[6]
    psi = matrix[7]

    # Plot each Metabolite
    plt.plot(t, Pyr, label='Pyr')
    plt.plot(t, NAD, label='NAD')
    plt.plot(t, OAA, label='OAA')
    plt.plot(t, AcCoA, label='AcCoA')
    plt.plot(t, Cit, label='Cit')
    plt.plot(t, KG, label='KG')
    plt.plot(t, ATP, label='ATP')
    plt.plot(t, psi, label='psi')
    plt.title(title)
    plt.legend()
    plt.xlabel('Tau')
    plt.show()

# This Function is for plotting the concentration of a particular metabolite across different HPIs (Model Runs)
def plot_hpi_concentration(y, t=[0, 24, 48, 72, 96], title=''):
    # Set the Plot Size
    f, ax = plt.subplots(figsize=(5, 3))

    # Make all the values realtive to the value at 0 HPI
    y = y / y[0]

    # Plot the curve
    plt.plot(t, y)
    plt.xlabel('Hours Post Infection (HPI)')
    plt.ylabel('Relative Concentration')
    plt.title(title)
    plt.show()

# This class contains the TCA model, and all the functions needed to run the model using different algorithms,
# as well as some auxiliary functions
class TCA_Cycle(object):
    # Upon initialization, take in a constant, concentration, and rate dict. Generally, the constant dict should
    # always be the same, but doing this allows for change
    def __init__(self, constant_dict, concentration_dict, rate_dict):
        # Store the constants
        self.At = constant_dict['At']
        self.K_app = constant_dict['K_app']
        self.asp = constant_dict['asp']
        self.glu = constant_dict['glu']
        self.Pi = constant_dict['Pi']
        self.Nt = constant_dict['Nt']
        self.C = constant_dict['C']
        self.a = constant_dict['a']
        self.b = constant_dict['b']
        self.psi_m = constant_dict['psi_m']
        self.K = constant_dict['K']
        self.F = constant_dict['F']
        self.R = constant_dict['R']
        self.T = constant_dict['T']
        self.K_6_eq = constant_dict['K_6_eq']
        # Calculate and store K_eq
        self.K_eq = self.K_6_eq * (self.glu / self.asp)

        # Load in the Concentrations
        Pyr = concentration_dict['Pyr']
        NAD = concentration_dict['NAD']
        OAA = concentration_dict['OAA']
        AcCoA = concentration_dict['AcCoA']
        Cit = concentration_dict['Cit']
        KG = concentration_dict['KG']
        ATP = concentration_dict['ATP']
        psi = concentration_dict['psi']

        # Load in the steady state concentrations
        Pyr_bar = concentration_dict['Pyr_bar']
        NAD_bar = concentration_dict['NAD_bar']
        OAA_bar = concentration_dict['OAA_bar']
        AcCoA_bar = concentration_dict['AcCoA_bar']
        Cit_bar = concentration_dict['Cit_bar']
        KG_bar = concentration_dict['KG_bar']
        ATP_bar = concentration_dict['ATP_bar']
        psi_bar = concentration_dict['psi_bar']

        # Calculate the scaled concentrations
        p = Pyr / Pyr_bar
        n = NAD / self.Nt
        o = OAA / OAA_bar
        a = AcCoA / AcCoA_bar
        c = Cit / Cit_bar
        k = KG / KG_bar
        e = ATP / self.At
        s = psi / self.psi_m

        # Store the scaled concentrations in the concentrations matrix
        self.concentration_matrix = np.array([[p], [n], [o], [a],
                                              [c], [k], [e], [s]])

        # Calculate the ettas
        self.etta_1 = AcCoA_bar / Pyr_bar
        self.etta_2 = Cit_bar / Pyr_bar
        self.etta_3 = KG_bar / Pyr_bar
        self.etta_4 = OAA_bar / Pyr_bar
        self.etta_5 = self.Nt / Pyr_bar
        self.etta_6 = self.At / Pyr_bar
        self.etta_7 = (self.psi_m / Pyr_bar) * self.C

        # Load in the Rate constants, commeted are the associated protein/process
        self.k1 = rate_dict['k1']  # Incoming Pyruvate
        self.k2 = rate_dict['k2']  # Pyruvate dehydrogenase
        self.k3 = rate_dict['k3']  # Citrate Synthase
        self.k4 = rate_dict['k4']  # Aconitase and Isocitrate dehydrogenase
        self.k5 = rate_dict['k5']  # Reversible a-ketoglutarate dehydrogenase
        self.k6 = rate_dict['k6']  # Reversible Succinyl-CoA synthetase, Succinic dehydrogenase, Fumarase, Malate dehydrogenase
        self.k7 = rate_dict['k7']  # Pyruvate carboxylase
        self.k8 = rate_dict['k8']  # OAA decay
        self.k_resp = rate_dict['k_resp'] # Resperation
        self.k_atp = rate_dict['k_atp'] # ATP production
        self.k_ant = rate_dict['k_ant'] # Exchange od ATP molecules across membrane
        self.k_leak = rate_dict['k_leak'] # Proton leak across membrane

        # Calculate the scaled rate constants
        self.B2 = (self.k2 / self.k1) * self.Nt * Pyr_bar
        self.B3 = (self.k3 / self.k1) * OAA_bar * AcCoA_bar
        self.B4 = (self.k4 / self.k1) * self.Nt * Cit_bar
        self.B5 = (self.k5 / self.k1) * self.Nt * self.At * KG_bar
        self.B6 = (self.k6 / self.k1) * OAA_bar
        self.B7 = (self.k7 / self.k1) * self.At * Pyr_bar
        self.B8 = (self.k8 / self.k1) * OAA_bar
        self.B_ant = (self.k_ant / self.k1) * self.At
        self.B_leak = (self.k_leak / self.k1) * self.psi_m
        self.B_resp = (self.k_resp / self.k1)
        self.B_atp = (self.k_atp / self.k1)
        self.delta_6 = KG_bar / (OAA_bar * self.K_eq)
        self.delta_r1 = self.K / self.Nt
        self.delta_r2 = self.a * self.psi_m
        self.delta_atp = self.b * self.At
        self.delta_crit = 3 * (1.2 * self.F * self.psi_m / (self.R * self.T))
        self.K_app_prime = self.K_app * self.Pi

        # Calculate a proto tau that will be used to calculate tau later on
        self.proto_tau = self.k1 / Pyr_bar

        # Calculate the reaction rates
        self.update_reaction_rates(self.concentration_matrix)

        # Set inital tau value
        self.tau = [0]

    # A function for clearing the concentration matrix and tau list
    def restart(self):
        self.concentration_matrix = self.concentration_matrix[:, 0:1]
        self.tau = self.tau[0]

    # Append a new column to the concentration matrix
    def update_concentration_matrix(self, new_matrix):
        self.concentration_matrix = np.concatenate((self.concentration_matrix, new_matrix), axis=-1)

    # Calculate the reaction rates using a concentration matrix
    def update_reaction_rates(self, concentration_matrix):
        # Get each current concentration value
        p = concentration_matrix[0][-1]
        n = concentration_matrix[1][-1]
        o = concentration_matrix[2][-1]
        a = concentration_matrix[3][-1]
        c = concentration_matrix[4][-1]
        k = concentration_matrix[5][-1]
        e = concentration_matrix[6][-1]
        s = concentration_matrix[7][-1]

        # Solve for each reaction rate
        self.v1 = 1
        self.v2 = self.B2 * p * n
        self.v3 = self.B3 * o * a
        self.v4 = self.B4 * c * n
        self.v5 = self.B5 * k * n * (1 - e)
        self.v6 = self.B6 * (o - self.delta_6 * k)
        self.v7 = self.B7 * p * e
        self.v8 = self.B8 * o
        self.v_ant = self.B_ant * e
        self.v_leak = self.B_leak * s
        self.v_resp = self.B_resp * ((1 - n) / (self.delta_r1 + 1 - n)) * (1 / (1 + np.exp(self.delta_r2 * (s - 1))))
        self.e_crit = (self.K_app_prime) / (self.K_app_prime + np.exp(-self.delta_crit * s))
        self.v_atp = self.B_atp * (2 / (1 + np.exp(self.delta_atp * (e - self.e_crit))) - 1)

    # Compute dy/dt for each concentration using the reaction rate
    def compute_dt(self):
        dp_dt = self.v1 - self.v2 - self.v7
        dn_dt = (-1 * self.v2 - self.v4 - 2 * self.v5 + self.v_resp) / self.etta_5
        do_dt = (self.v5 + self.v7 - self.v3 - self.v8 - self.v6) / self.etta_4
        da_dt = (self.v2 - self.v3) / self.etta_1
        dc_dt = (self.v3 - self.v4) / self.etta_2
        dk_dt = (self.v4 + self.v6 - self.v5) / self.etta_3
        de_dt = (self.v_atp - self.v_ant + self.v5 - self.v7) / self.etta_6
        ds_dt = (10 * self.v_resp - 3 * self.v_atp - self.v_leak - self.v_ant) / self.etta_7

        # Store values in matrix and return that matrix
        dt_matrix = np.array([[dp_dt], [dn_dt], [do_dt], [da_dt],
                              [dc_dt], [dk_dt], [de_dt], [ds_dt]])

        return dt_matrix

    # Implementation of Forward Euler algorithm
    def Forward_Euler(self, dt):
        # Calculate dtau
        dtau = dt * self.proto_tau
        # Get the current concentration column
        concentrations = self.concentration_matrix[:, -1:]
        # Compute the new concentration column
        new_concentration_matrix = concentrations + dtau * self.compute_dt()
        # Update the concentration matrix, reaction rates, and tau list
        self.update_concentration_matrix(new_concentration_matrix)
        self.update_reaction_rates(self.concentration_matrix)
        self.tau.append(self.tau[-1] + dt * self.proto_tau)

    # Helper function for RK4 algorithm that computes dy/dt at a particular point
    def compute_F(self, concentration_matrix, dt):
        # update the reaction rates for the current concentration values
        self.update_reaction_rates(concentration_matrix)

        # calculate concentration at the new time point
        concentration_matrix = np.array(concentration_matrix + dt * self.compute_dt())

        # update the reaction rates for the new concentration values
        self.update_reaction_rates(concentration_matrix)

        # Calculate d/dt for the new concentration values
        d_dt = self.compute_dt()

        return d_dt

    #  Implementation of the RK4 algorithm
    def RK4(self, dt):
        # Calculate dtau
        dtau = np.array(dt * self.proto_tau)
        # Get the current concentration column
        concentrations = np.array(self.concentration_matrix[:, -1:])
        # Calculate the xi values
        xi_1 = concentrations

        xi_2 = concentrations + np.array(dt / 2 * self.compute_F(xi_1, 0))

        xi_3 = concentrations + np.array(dt / 2 * self.compute_F(xi_2, dtau / 2))

        xi_4 = concentrations + np.array(dt * self.compute_F(xi_3, dtau / 2))

        # Use the xi values to compute the new concentrations
        new_concentration_matrix = np.array(concentrations + dtau / 6 * (self.compute_F(xi_1, 0)
                                                                         + 2 * self.compute_F(xi_2, dtau / 2)
                                                                         + 2 * self.compute_F(xi_3, dtau / 2)
                                                                         + self.compute_F(xi_4, dtau)))
        # Update the concentration matrix, reaction rates, and tau list
        self.update_concentration_matrix(new_concentration_matrix)
        self.update_reaction_rates(self.concentration_matrix)
        self.tau.append(self.tau[-1] + dt * self.proto_tau)

    # A helper function for the Implicit algorithm that computes the Jacobian
    def compute_Jacobian(self, concentration_matrix):
        # Get each current concentration value
        p = concentration_matrix[0][0]
        n = concentration_matrix[1][0]
        o = concentration_matrix[2][0]
        a = concentration_matrix[3][0]
        c = concentration_matrix[4][0]
        k = concentration_matrix[5][0]
        e = concentration_matrix[6][0]
        s = concentration_matrix[7][0]

        # Calculate each row of the Jacobian
        dp = np.array([[-self.B2 * n - self.B7 * e, -self.B2 * p, 0, 0, 0, 0, -self.B7 * p, 0]], dtype=float)
        dn = np.array([[-self.B2 * n / self.etta_5, (-self.B2 * p - self.B4 * c - 2 * self.B5 * k * (1 - e) - (
                    self.B_resp / ((-n + self.delta_r1 + 1) * (np.exp((s - 1) * self.delta_r2)) + 1)) + (
                                                                 ((1 - n) * self.B_resp) / (
                                                                     ((-n + self.delta_r1 + 1) ** 2) * (
                                                                 np.exp((s - 1) * self.delta_r2)) + 1))) / self.etta_5,
                        0, 0, -self.B4 * n / self.etta_5, 2 * self.B5 * n * (e - 1) / self.etta_5,
                        2 * self.B5 * k * n / self.etta_5,
                        -((1 - n) * self.delta_r2 * self.B_resp * np.exp((s - 1) * self.delta_r2)) / (
                                    self.etta_5 * (-n + self.delta_r1 + 1) * (
                                        np.exp((s - 1) * self.delta_r2) + 1) ** 2)]], dtype=float)
        do = np.array([[self.B7 * e / self.etta_4, -self.B5 * k * (e - 1) / self.etta_4,
                        (-self.B5 * self.B3 * a + self.B6 + self.B8) / self.etta_4,
                        -self.B5 * self.B3 * o / self.etta_4, 0,
                        (self.B5 * n * e + self.B5 * n + self.B6 * self.delta_6) / self.etta_4,
                        (self.B7 * p - self.B5 * k * n) / self.etta_4, 0]], dtype=float)
        da = np.array([[self.B2 * n / self.etta_1, self.B2 * p / self.etta_1, -self.B3 * a / self.etta_1,
                        -self.B3 * o / self.etta_1, 0, 0, 0, 0]], dtype=float)
        dc = np.array([[0, -self.B4 * c / self.etta_2, self.B3 * a / self.etta_2, self.B3 * o / self.etta_2,
                        -self.B4 * n / self.etta_2, 0, 0, 0]], dtype=float)
        dk = np.array([[0, (self.B4 * c + (e - 1) * self.B5 * k) / self.etta_3, self.B6 / self.etta_3, 0,
                        self.B4 * n / self.etta_3, ((e - 1) * self.B5 * n - self.B6 * self.delta_6) / self.etta_3,
                        self.B5 * k * n / self.etta_3, 0]], dtype=float)
        de = np.array([[-self.B7 * e / self.etta_6, -self.B5 * k * (e - 1) / self.etta_6, 0, 0, 0,
                        -self.B5 * n * (e - 1) / self.etta_6, (-(
                        (2 * self.B_atp * self.delta_atp * np.exp(self.delta_atp * (e - self.e_crit))) / (np.exp(
                    self.delta_atp * (
                                e - self.e_crit)) + 1) ** 2) - self.B5 * k * n - self.B7 * p - self.B_ant) / self.etta_6,
                        (2 * self.B_atp * self.delta_crit * self.delta_atp * self.K_app_prime * np.exp(
                            self.delta_atp * (e - self.e_crit) - self.delta_crit * s)) / (
                                    self.etta_6 * ((np.exp(-self.delta_crit * s) + self.K_app_prime) ** 2) * (
                                        np.exp(self.delta_atp * (e - self.e_crit)) + 1) ** 2)]], dtype=float)
        ds = np.array([[0, -10 * self.B_resp * self.delta_r1 / (
                    self.etta_7 * ((n - self.delta_r1 - 1) ** 2) * (np.exp(self.delta_r2 * (s - 1)) + 1)), 0, 0, 0, 0, (
                                    ((6 * self.delta_atp * self.B_atp * np.exp(self.delta_atp * (e - self.e_crit))) / (
                                                np.exp(self.delta_atp * (
                                                            e - self.e_crit)) + 1) ** 2) - self.B_ant) / self.etta_7, -(
                    (10 * self.v_resp * (1 - n) * self.delta_r2 * np.exp(self.delta_r2 * (s - 1))) / (
                        (-n + self.delta_r1 + 1) * (np.exp(self.delta_r2 * (s - 1)) + 1) ** 2) - self.B_leak - (
                                6 * self.delta_crit * self.delta_atp * self.B_atp * self.K_app_prime * np.exp(
                            self.delta_atp * (e - self.e_crit) - self.delta_crit * s)) / (
                                ((np.exp(-self.delta_crit * s) + self.K_app_prime) ** 2) * (
                                    np.exp(self.delta_atp * (e - self.e_crit)) + 1) ** 2)) / self.etta_7]], dtype=float)

        # Combine the rows into the Jacobian matrix then return it
        J = np.concatenate((dp, dn, do, da, dc, dk, de, ds), axis=0)

        return J

    #Implementation of the Implicit method algorithm
    def implicit(self, dt):
        # Get each current concentration value
        dtau = np.array(dt * self.proto_tau)
        # Compute the current d/dt matrix
        dt_matrix = np.array(self.compute_dt())
        # Get the current concentration column
        concentrations = np.array(self.concentration_matrix[:, -1:])
        # Compute the Jacobian
        J = self.compute_Jacobian(concentrations)
        # Compute the inverted Jacobian
        J_inv = inv(np.identity(8) - dtau * J)
        # Compute the new concentrations
        inner_sum = np.array(concentrations + dtau * dt_matrix - np.matmul(J * dtau, concentrations))
        new_concentration_matrix = np.matmul(J_inv, inner_sum)
        # Update the concentration matrix, reaction rates, and tau list
        self.update_concentration_matrix(new_concentration_matrix)
        self.update_reaction_rates(self.concentration_matrix)
        self.tau.append(self.tau[-1] + dt * self.proto_tau)

    # A function for running the model for any of the algorithms
    def run(self, dt, tf, method='RK4', drop_freq=2):
        # Keep track of iterations to know when to drop values to save on storage space and run faster
        iteration = 0
        # Keep running until the final tau value is reached or passed
        while self.tau[-1] < tf:
            # Run the desired method with the given dt
            if method == 'RK4':
                self.RK4(dt)
            elif method == 'implicit':
                self.implicit(dt)
            else:
                self.Forward_Euler(dt)

            # Update the iteration counter
            iteration = iteration + 1
            # If we have reached the drop frequency, drop the second to last column of the concentration matrix
            # and the second to last value in the tau list
            if iteration % drop_freq == 0:
                self.concentration_matrix = np.delete(self.concentration_matrix, obj=-2, axis=1)
                del self.tau[-2]



    # A function to save the concentration matrix and tau list as a pickle file
    def save(self, save_name):
        dicty = {
            'concentration_matrix': self.concentration_matrix,
            'tau': self.tau
        }
        save_object(dicty, save_name)

    # A function to plot the concentration matrix
    def plot(self, title=''):
        # Set the figure size
        f, ax = plt.subplots(figsize=(10, 7))
        # Extract the row corresponding to each of the Metabolites
        Pyr = self.concentration_matrix[0]
        NAD = self.concentration_matrix[1]
        OAA = self.concentration_matrix[2]
        AcCoA = self.concentration_matrix[3]
        Cit = self.concentration_matrix[4]
        KG = self.concentration_matrix[5]
        ATP = self.concentration_matrix[6]
        psi = self.concentration_matrix[7]
        # Plot each Metabolite
        plt.plot(self.tau, Pyr, label='Pyr')
        plt.plot(self.tau, NAD, label='NAD')
        plt.plot(self.tau, OAA, label='OAA')
        plt.plot(self.tau, AcCoA, label='AcCoA')
        plt.plot(self.tau, Cit, label='Cit')
        plt.plot(self.tau, KG, label='KG')
        plt.plot(self.tau, ATP, label='ATP')
        plt.plot(self.tau, psi, label='psi')
        plt.legend()
        plt.title(title)
        plt.xlabel('Tau')
        plt.show()

    # Calculate the error curve for psi for a given method and range
    def calculate_error_curve(self, dts, dtf, dt_inc, tf, method='RK4', p=4):
        # Keep track of error and dt
        eh_list = []
        dt_list = []
        dt = dts
        while dt < dtf:
            # Clear out old concentraion matrix and tau list
            self.restart()
            # Run with grid spacing 2h
            self.run(dt, tf, method=method, drop_freq=999999999)
            # Get value of psi at tf
            h = self.concentration_matrix[7, -1]


            # Clear out old concentraion matrix and tau list
            self.restart()
            # Run with grid spacing h
            self.run(dt / 2, tf, method=method, drop_freq=999999999)
            # Get value of psi at tf
            h_half = self.concentration_matrix[7, -1]

            # Calculate the empirical error estimate
            eh = (h - h_half) / (2 ** p - 1)

            # Keep track of grid spacing and error
            eh_list.append(eh)
            dt_list.append(dt)

            # update the current grid spacing
            dt = dt + dt_inc

        return eh_list, dt_list

    # A function for plotting the error curves for each of the three algorithms
    def plot_error_curves(self, dts, dtf, dt_inc, tf):
        # Calculate error curves for each algorithm
        eh_list_f, dt_list_f = self.calculate_error_curve(dts, dtf, dt_inc, tf, method='F', p=1)
        eh_list_RK4, dt_list_RK4 = self.calculate_error_curve(dts, dtf, dt_inc, tf, method='RK4', p=4)
        eh_list_implicit, dt_list_implicit = self.calculate_error_curve(dts, dtf, dt_inc, tf, method='implicit', p=1)

        # Plot the curves
        plt.figure()
        plt.plot(dt_list_f, eh_list_f, label='Forward Euler')
        plt.plot(dt_list_RK4, eh_list_RK4, label='RK4')
        plt.plot(dt_list_implicit, eh_list_implicit, label='Implicit')
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel('dt')
        plt.ylabel('Empirical Estimate of Error of delta_psi')
        plt.title('Emperical Estimate of Error tf=' + str(tf))
        plt.legend()
        plt.show()