import scipy.optimize as opt
from ..core.sim_cycle import SimCycle

from hx_UA_const.components.compressor import Compressor
from hx_UA_const.components.expansion_valve import ExpansionValve
from hx_UA_const.components.heat_exchanger import Condenser, Evaporator
from hx_UA_const.components.connector import Connector

from hx_UA_const.metrics.dsh_dsc_cal import DSHCalculator
from hx_UA_const.metrics.charge_cal import ChargeCalculator
from hx_UA_const.metrics.mdot_cal import mdotCalculator

from dataclasses import dataclass

@dataclass
class solved_cond_results:
    P_cond_sol: float
    h_comp_out: float
    s_comp_out: float
    T_comp_out: float
    h_cond_out: float
    s_cond_out: float
    T_cond_out: float
    h_exp_out: float
    s_exp_out: float
    T_exp_out: float
    mdot: float
    m_cond: float

@dataclass
class solved_eva_results(solved_cond_results):
    P_eva_sol: float
    h_eva_out: float
    s_eva_out: float
    T_eva_out: float
    m_eva: float

@dataclass
class solved_results(solved_eva_results):
    mtot: float
    DSH_sol: float



class PressureSolver_mdot_DSH_charge:
    def __init__(self,
                 sim:SimCycle,
                 params):
        self.sim = sim
        self.params = params

        self.comp = Compressor(sim, self.params)
        self.cond = Condenser(sim, self.params)
        self.exp = ExpansionValve(sim, self.params)
        self.eva = Evaporator(sim, self.params)
        self.conn = Connector(sim, self.params)
        self.dsh = DSHCalculator(sim, self.params.DSH_target)
        self.charge = ChargeCalculator(sim, self.params.charge_target)
        self.Mdot = mdotCalculator(sim)

        self.tol = params.tol
 

    def solve_cond(self, P_eva: float, T_cond_air: float):
        def cycle_mdot(P_cond):
            h_comp_out, s_comp_out, T_comp_out, mdot_comp = self.comp.process(P_eva, P_cond)
            h_cond_out, s_cond_out, T_cond_out, m_cond = self.cond.exchange(mdot_comp, P_cond, h_comp_out)
            h_exp_out, s_exp_out, T_exp_out, mdot_exp = self.exp.process(P_eva, P_cond, h_cond_out)
            return self.Mdot.error(mdot_comp, mdot_exp)
        
        # bisect or brentq or toms748
        P_cond_low = self.sim.get_single('QT_inputs', 0, T_cond_air, ('P'))
        P_cond_high = self.sim.get_single('QT_inputs', 0, T_cond_air + 30, ('P'))
        P_cond_sol = opt.brentq(cycle_mdot, P_cond_low, P_cond_high, xtol=self.tol)

        h_comp_out, s_comp_out, T_comp_out, mdot = self.comp.process(P_eva, P_cond_sol)
        h_cond_out, s_cond_out, T_cond_out, m_cond = self.cond.exchange(mdot, P_cond_sol, h_comp_out)
        h_exp_out, s_exp_out, T_exp_out, _ = self.exp.process(P_eva, P_cond_sol, h_cond_out)

        return solved_cond_results(
            P_cond_sol=P_cond_sol,
            h_comp_out=h_comp_out,
            s_comp_out=s_comp_out,
            T_comp_out=T_comp_out,
            h_cond_out=h_cond_out,
            s_cond_out=s_cond_out,
            T_cond_out=T_cond_out,
            h_exp_out=h_exp_out,
            s_exp_out=s_exp_out,
            T_exp_out=T_exp_out,
            mdot=mdot,
            m_cond=m_cond
        )
            
    def solve_evap(self, T_cond_air: float, T_eva_air: float):
        def cycle_DSH(P_eva):
            solved_cond = self.solve_cond(P_eva, T_cond_air)
            h_eva_out, s_eva_out, T_eva_out, m_eva = self.eva.exchange(solved_cond.mdot, P_eva, solved_cond.h_exp_out)
            return self.dsh.error(T_eva_out, P_eva)
        
        # bisect or brentq or toms748
        P_eva_low = self.sim.get_single('QT_inputs', 1, T_eva_air - 30, ('P'))
        P_eva_high = self.sim.get_single('QT_inputs', 1, T_eva_air - 10, ('P'))
        P_eva_sol = opt.brentq(cycle_DSH, P_eva_low, P_eva_high, xtol=self.tol)

        solved_cond = self.solve_cond(P_eva_sol, T_cond_air)
        h_eva_out, s_eva_out, T_eva_out, m_eva = self.eva.exchange(solved_cond.mdot, P_eva_sol, solved_cond.h_exp_out)

        return solved_eva_results(
            **vars(solved_cond),
            P_eva_sol=P_eva_sol,
            h_eva_out=h_eva_out,
            s_eva_out=s_eva_out,
            T_eva_out=T_eva_out,
            m_eva=m_eva
        )
    
    def solve_charge(self, T_cond_air: float, T_eva_air: float):
        def cycle_charge(DSH_ass):
            self.params.DSH_target = DSH_ass
            self.comp.set_target(DSH_ass)
            self.dsh.set_target(DSH_ass)

            solved_eva = self.solve_evap(T_cond_air, T_eva_air)
            m_conn = self.conn.process(solved_eva.P_eva_sol, solved_eva.h_cond_out, solved_eva.h_eva_out)
            mtot = solved_eva.m_cond + solved_eva.m_eva + m_conn
            return self.charge.error(mtot)
        
        # bisect or brentq or toms748
        DSH_low = 1e-3
        DSH_high = 20
        DSH_sol = opt.brentq(cycle_charge, DSH_low, DSH_high, xtol=self.tol)
        
        self.params.DSH_target = DSH_sol
        self.comp.set_target(DSH_sol)
        self.dsh.set_target(DSH_sol)

        solved_eva = self.solve_evap(T_cond_air, T_eva_air)
        m_conn = self.conn.process(solved_eva.P_eva_sol, solved_eva.h_cond_out, solved_eva.h_eva_out)
        mtot = solved_eva.m_cond + solved_eva.m_eva + m_conn

        return solved_results(
            **vars(solved_eva),
            mtot=mtot,
            DSH_sol=DSH_sol
        )
    

    