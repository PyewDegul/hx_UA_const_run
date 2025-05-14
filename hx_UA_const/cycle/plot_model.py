import numpy as np
import matplotlib.pyplot as plt

class PlotModel:
    def __init__(self, sim, params, res):
        self.sim = sim
        self.params = params
        self.res = res

        '''TS_diagram, PH_diagram'''
        # self.sim에 dependent

        '''찍혀야 하는 점'''
        # P_con_sol, P_eva_sol
        # 1 Point 1 //  (P_eva_sol, h_comp_in, s_comp_in, T_comp_in = T_eva_elem[-1] + DSH)
        # 2.Point 1 // (P_cond_sol, h_comp_out, s_comp_out, T_comp_out)
        # 3.Point N_cond // (P_cond_sol, h_cond_elem, s_cond_elem, T_cond_elem)
        # 4.Point 1 // (P_eva_sol, h_exp_out, s_exp_out, T_exp_out)
        # 5.Point N_eva // (P_eva_sol, h_eva_elem, s_eva_elem, T_eva_elem)

        
    # def plot_TS(self):
        # Plot Ts diagram

    # def plot_PH(self):
        # Plot Ph diagram