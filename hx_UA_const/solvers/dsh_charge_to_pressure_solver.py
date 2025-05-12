import scipy.optimize as opt
from ..core.sim_cycle import SimCycle

from hx_UA_const.components.compressor import Compressor
from hx_UA_const.components.expansion_valve import ExpansionValve
from hx_UA_const.components.heat_exchanger import Condenser, Evaporator
from hx_UA_const.components.connector import Connector

from hx_UA_const.metrics.dsh_dsc_cal import DSHCalculator
from hx_UA_const.metrics.charge_cal import ChargeCalculator

class PressureSolver_charge:
    def __init__(self,
                 sim:SimCycle,
                 params):
        self.sim = sim
        self.params = params

        self.comp = Compressor(sim, self.params)
        self.cond = Condenser(sim, self.params.N_cond,
                              self.params.UA_cond / self.params.N_cond,
                              self.params.T_cond_air, self.params.V_elem_cond)
        self.exp = ExpansionValve(sim, self.params)
        self.eva = Evaporator(sim, self.params.N_eva,
                              self.params.UA_eva / self.params.N_eva,
                              self.params.T_eva_air, self.params.V_elem_eva)
        self.conn = Connector(sim, self.params)
        self.dsh = DSHCalculator(sim, self.params.DSH_target)
        self.charge = ChargeCalculator(sim, self.params.charge_target)

        self.tol = self.params.tol


    def solve_evap(self, P_cond: float, T_eva_air: float):
        def cycle_dsh(P_eva):
            h_comp_out, s_comp_out, T_comp_out, mdot = self.comp.process(P_eva, P_cond)
            h_cond_out, s_cond_out, T_cond_out, m_cond = self.cond.exchange(mdot, P_cond, h_comp_out)
            # h_exp_out, s_exp_out, T_exp_out = self.exp.process(P_eva, P_cond, h_cond_out)
            h_eva_out, s_eva_out, T_eva_out, m_eva = self.eva.exchange(mdot, P_eva, h_cond_out)
            
            return self.dsh.error(T_eva_out, P_eva)
        
        # bisect or brentq or toms748
        
        P_eva_low = self.sim.get_single('QT_inputs', 1, T_eva_air - 20, ('P'))
        P_eva_high = self.sim.get_single('QT_inputs', 1, T_eva_air, ('P'))
        P_eva_sol = opt.brentq(cycle_dsh, P_eva_low, P_eva_high, xtol=self.tol)

        h_comp_out, s_comp_out, T_comp_out, mdot = self.comp.process(P_eva_sol, P_cond)
        h_cond_out, s_cond_out, T_cond_out, m_cond = self.cond.exchange(mdot, P_cond, h_comp_out)
        # h_exp_out, s_exp_out, T_exp_out = self.exp.process(P_eva_sol, P_cond, h_cond_out)
        h_eva_out, s_eva_out, T_eva_out, m_eva = self.eva.exchange(mdot, P_eva_sol, h_cond_out)

        return P_eva_sol, h_comp_out, s_comp_out, T_comp_out, h_cond_out, s_cond_out, T_cond_out, h_eva_out, s_eva_out, T_eva_out, m_cond, m_eva, mdot
    # return P_eva_sol, h_comp_out, s_comp_out, T_comp_out, h_cond_out, s_cond_out, T_cond_out, h_exp_out, s_exp_out, T_exp_out, h_eva_out, s_eva_out, T_eva_out, m_cond, m_eva, mdot
        
    def solve_cond(self, T_cond_air: float, T_eva_air: float):
        def cycle_charge(P_cond):
            P_eva_sol, h_comp_out, s_comp_out, T_comp_out, h_cond_out, s_cond_out, T_cond_out, h_eva_out, s_eva_out, T_eva_out, m_cond, m_eva, mdot = self.solve_evap(P_cond, T_eva_air)
            # P_eva_sol, h_comp_out, s_comp_out, T_comp_out, h_cond_out, s_cond_out, T_cond_out, h_exp_out, s_exp_out, T_exp_out, h_eva_out, s_eva_out, T_eva_out, m_cond, m_eva, mdot  = self.solve_evap(P_cond, T_eva_air)
            m_conn = self.conn.process(P_eva_sol, h_cond_out, h_eva_out)
            mtot = m_cond + m_eva + m_conn
            return self.charge.error(mtot)
        
        # bisect or brentq or toms748
        P_cond_low = self.sim.get_single('QT_inputs', 0, T_cond_air, ('P'))
        P_cond_high = self.sim.get_single('QT_inputs', 0, T_cond_air + 30, ('P'))

        P_cond_sol = opt.brentq(cycle_charge, P_cond_low, P_cond_high, xtol=self.tol)
        P_eva_sol, h_comp_out, s_comp_out, T_comp_out, h_cond_out, s_cond_out, T_cond_out, h_eva_out, s_eva_out, T_eva_out, m_cond, m_eva, mdot = self.solve_evap(P_cond_sol, T_eva_air)
        # P_eva_sol, h_comp_out, s_comp_out, T_comp_out, h_cond_out, s_cond_out, T_cond_out, h_exp_out, s_exp_out, T_exp_out, h_eva_out, s_eva_out, T_eva_out, m_cond, m_eva, mdot  = self.solve_evap(P_cond, T_eva_air)
        m_conn = self.conn.process(P_eva_sol, h_cond_out, h_eva_out)
        mtot = m_cond + m_eva + m_conn
        
        return P_cond_sol, P_eva_sol
        # return P_cond_sol, P_eva_sol, h_comp_out, s_comp_out, T_comp_out, h_cond_out, s_cond_out, T_cond_out, h_eva_out, s_eva_out, T_eva_out, mtot, mdot