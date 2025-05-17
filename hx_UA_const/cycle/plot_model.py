import numpy as np
import matplotlib.pyplot as plt

import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(__file__, '..', '..', '..')))
from hx_UA_const.core.plot_tool import Graph_class

class PlotModel():
    def __init__(self, sim, params, res):
        self.sim = sim
        self.state = self.sim.get_state()
        self.params = params
        self.res = res

        self.backend = self.state.backend_name()
        self.fluid = self.state.fluid_names()

        '''TS_diagram, PH_diagram'''
        # self.sim에 dependent

        '''찍혀야 하는 점'''
        # P_con_sol, P_eva_sol
        # 1 Point 1 //  (P_eva_sol, h_comp_in, s_comp_in, T_comp_in = T_eva_elem[-1] + DSH)
        # 2.Point 1 // (P_cond_sol, h_comp_out, s_comp_out, T_comp_out)
        # 3.Point N_cond // (P_cond_sol, h_cond_elem, s_cond_elem, T_cond_elem)
        # 4.Point 1 // (P_eva_sol, h_exp_out, s_exp_out, T_exp_out)
        # 5.Point N_eva // (P_eva_sol, h_eva_elem, s_eva_elem, T_eva_elem)

    ''' Ts, Ph diagram 더 정확하게(critical 인근에서 구하는 법) 구현 // 작동 유체 온도 역시 함께 표시하기 '''
    def plot_TS(self):
        T_min = self.sim.T_TP + 1e-6
        n_points = 500
        T_vals = np.linspace(T_min, self.sim.T_C * 0.999, n_points)
        sat_liq_s = np.zeros(n_points)
        sat_vap_s = np.zeros(n_points)

        for i in range(n_points):
            sat_liq_s[i] = self.sim.get_single('QT_inputs', 0, T_vals[i], ('S')) / 1000.0 # entropy [kJ/kgK]
            sat_vap_s[i] = self.sim.get_single('QT_inputs', 1, T_vals[i], ('S')) / 1000.0 # entropy [kJ/kgK]
        
        T_vals_C = T_vals # [K]
        '''
        # Critical point properties for marking
        s_crit = self.sim.get_single('PT_inputs', self.sim.P_C, self.sim.T_C, ('S'))  / 1000.0 # entropy [kJ/kgK]  
        T_crit_C = self.sim.T_C  # Critical temperature in K
        '''
        # Prepare cycle data points (convert to proper units)
        # Evaporator exit / Compressor inlet (Point 1)
        T_comp_in_C = self.res.T_eva_elem[-1]  # [K]
        s_comp_in = self.res.s_eva_elem[-1] / 1000.0    # kJ/kg-K
        # Compressor outlet (Point 2)
        T_comp_out_C = self.res.T_comp_out     # [K]
        s_comp_out = self.res.s_comp_out / 1000.0       # kJ/kg-K
        # Expansion valve outlet (Point 4)
        T_exp_out_C = self.res.T_exp_out      # [K]
        s_exp_out = self.res.s_exp_out / 1000.0         # kJ/kg-K
        # Combine points for the condenser segment (Point 2 -> 3) and evaporator segment (Point 4 -> 5)
        # Convert entire arrays to lists in proper units
        T_cond_C_list = list(self.res.T_cond_elem)   # [K]
        s_cond_list = list(self.res.s_cond_elem / 1000.0)     # kJ/kg-K
        T_eva_C_list = list(self.res.T_eva_elem)     # [K]
        s_eva_list = list(self.res.s_eva_elem / 1000.0)       # kJ/kg-K

        # Build the complete cycle path in order: comp inlet -> comp outlet -> condenser -> exp valve outlet -> evaporator -> back to comp inlet
        cycle_T_C = [T_comp_in_C, T_comp_out_C] + T_cond_C_list + [T_exp_out_C] + T_eva_C_list
        cycle_s   = [s_comp_in,   s_comp_out]   + s_cond_list   + [s_exp_out]   + s_eva_list

        # Plotting
        Graph = Graph_class()
        fig, ax = plt.subplots(figsize=(6,4))
        Graph.setPlotStyle(ax=ax, title='T-s Diagram', 
                        x_label='Entropy [kJ/(kg-K)]', 
                        y_label='Temperature [°C]')  # apply styling:contentReference[oaicite:12]{index=12}
        # Plot saturation lines
        ax.plot(sat_liq_s, T_vals_C, color='black', linestyle='-', label='Sat. Liquid')
        ax.plot(sat_vap_s, T_vals_C, color='black',  linestyle='-', label='Sat. Vapor')
        # Plot critical point
        # ax.scatter(s_crit, T_crit_C, color='red', zorder=3, label='Critical Point')
        # Plot cycle process line
        ax.plot(cycle_s, cycle_T_C, color='black', markersize=1, label='Cycle Path')
        # Add legend (optional, to identify lines)
        ax.legend(loc='best')
        # Show and save figure
        plt.show()

        # 이미지 저장
        script_dir = os.path.dirname(os.path.abspath(__file__))
        save_dir = os.path.join(script_dir, "Diagrams")
        os.makedirs(save_dir, exist_ok=True)               # 폴더 없으면 생성
        fig_filename = os.path.join(save_dir, f"{self.backend}_{self.fluid}_TS_diagram.png")
        fig.savefig(fig_filename)

    def plot_PH(self):
        T_min = self.sim.T_TP + 1e-6
        n_points = 500
        T_vals = np.linspace(T_min, self.sim.T_C * 0.999, n_points)
        sat_liq_h= np.zeros(n_points)
        sat_liq_p = np.zeros(n_points)
        sat_vap_h = np.zeros(n_points)
        sat_vap_p = np.zeros(n_points)

        for i in range(n_points):
            sat_liq_h[i], sat_liq_p[i] = self.sim.get_multiple('QT_inputs', 0, T_vals[i], ('H', 'P')) # [J/kg], [Pa]
            sat_vap_h[i], sat_vap_p[i] = self.sim.get_multiple('QT_inputs', 1, T_vals[i], ('H', 'P')) # [J/kg], [Pa]
        
        sat_liq_h = sat_liq_h / 1000.0 # [kJ/kg]
        sat_liq_p = sat_liq_p / 1e6 # [kJ/kgK]
        sat_vap_h = sat_vap_h / 1000.0 # [kJ/kg]
        sat_vap_p = sat_vap_p / 1e6 # [kJ/kgK]

        '''
        # Critical point properties for marking
        h_crit = self.sim.get_single('PT_inputs', self.sim.P_C, self.sim.T_C, ('H')) / 1000.0 # entropy [kJ/kg]  
        P_crit = self.sim.P_C / 1e6 # critical pressure [MPa]
        '''
        # Cycle points (enthalpy in kJ/kg, pressure in MPa)
        # Evaporator exit / Compressor inlet
        h_comp_in = self.res.h_eva_elem[-1] / 1000.0    # kJ/kg
        P_comp_in = self.res.P_eva_sol / 1e6            # MPa
        # Compressor outlet
        h_comp_out = self.res.h_comp_out / 1000.0       # kJ/kg
        P_comp_out = self.res.P_cond_sol / 1e6          # MPa
        # Expansion valve outlet
        h_exp_out = self.res.h_exp_out / 1000.0         # kJ/kg
        P_exp_out = self.res.P_eva_sol / 1e6            # MPa (back to evaporator pressure)
        # Prepare condenser and evaporator segment arrays/lists
        h_cond_list = list(self.res.h_cond_elem / 1000.0)   # kJ/kg
        P_cond_list = [self.res.P_cond_sol / 1e6] * len(h_cond_list)  # MPa (constant for all condenser points)
        h_eva_list = list(self.res.h_eva_elem / 1000.0)     # kJ/kg
        P_eva_list = [self.res.P_eva_sol / 1e6] * len(h_eva_list)     # MPa (constant for all evaporator points)
        
        # Build combined cycle path: comp in -> comp out -> condenser -> exp valve -> evaporator -> back to comp in
        cycle_h = [h_comp_in, h_comp_out] + h_cond_list + [h_exp_out] + h_eva_list
        cycle_P = [P_comp_in, P_comp_out] + P_cond_list + [P_exp_out] + P_eva_list

        # Plotting
        Graph = Graph_class()
        fig, ax = plt.subplots(figsize=(6,4))
        Graph.setPlotStyle(ax=ax, title='P-h Diagram',
                        x_label='Enthalpy [kJ/kg]', 
                        y_label='Pressure [MPa]')
        # Plot saturation curves in P-h coordinates
        ax.plot(sat_liq_h, sat_liq_p, color='black', linestyle='-', label='Sat. Liquid')
        ax.plot(sat_vap_h, sat_vap_p, color='black',  linestyle='-', label='Sat. Vapor')
        # Plot critical point
        # ax.scatter(h_crit, P_crit, color='red', zorder=3, label='Critical Point')
        # Plot cycle path
        ax.plot(cycle_h, cycle_P, color='black', markersize=1, label='Cycle Path')
        # Set y-axis to log scale for better visualization (optional but common for P-h):contentReference[oaicite:16]{index=16}
        ax.set_yscale('log')
        ax.legend(loc='best')
        plt.show()

        script_dir = os.path.dirname(os.path.abspath(__file__))
        save_dir = os.path.join(script_dir, "Diagrams")
        os.makedirs(save_dir, exist_ok=True)               # 폴더 없으면 생성
        fig_filename = os.path.join(save_dir, f"{self.backend}_{self.fluid}_PH_diagram.png")
        fig.savefig(fig_filename)

