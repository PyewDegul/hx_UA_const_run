import scipy.optimize as opt
from ..core.sim_cycle import SimCycle

class PressureSolver_charge:
    def __init__(self,
                 sim:SimCycle,
                 compressor,
                 condenser,
                 evaporator,
                 connector,
                 dsh_calc,
                 charge_calc,
                 tol: float):
        self.comp = compressor
        self.cond = condenser
        self.eva = evaporator
        self.conn = connector
        self.dsh = dsh_calc
        self.charge = charge_calc
        self.tol = tol
        self.sim = sim

    def solve_evap(self, P_cond: float, T_eva_air: float):
        def cycle_dsh(P_eva):
            # 1) 압력 출력 (Pa → bar)
            # print(f"P_eva = {P_eva/1e6:.3f} MPa, P_cond = {P_cond/1e6:.3f} MPa")

            # 2) 압축기 출구 엔탈피 & 유량 (h: J/kg → kJ/kg, mdot: kg/s)
            mdot, h_comp_out = self.comp.process(P_eva, P_cond)
            # print(f"mdot = {mdot:.4f} kg/s, h_comp_out = {h_comp_out/1e3:.2f} kJ/kg")

            # 3) 응축기 출구 상태 (h: J/kg → kJ/kg, T: K → °C, m_cond: kg/s)
            h_cond_out, T_cond_out, m_cond = self.cond.exchange(mdot, P_cond, h_comp_out)
            '''print(f"h_cond_out = {h_cond_out/1e3:.2f} kJ/kg, "
                f"T_cond_out = {T_cond_out-273.15:.2f} °C, "
                f"m_cond = {m_cond:.4f} kg")
            '''
            # 4) 증발기 출구 상태 (h: J/kg → kJ/kg, T: K → °C, m_eva: kg/s)
            h_eva_out, T_eva_out, m_eva = self.eva.exchange(mdot, P_eva, h_cond_out)
            '''print(f"h_eva_out = {h_eva_out/1e3:.2f} kJ/kg, "
                f"T_eva_out = {T_eva_out-273.15:.2f} °C, "
                f"m_eva = {m_eva:.4f} kg")'''
            # 최종 에러 리턴
            return self.dsh.error(T_eva_out, P_eva)
        
        # bisect or brentq or toms748
        
        P_eva_low = self.sim.get_single('QT_inputs', 1, T_eva_air - 20, ('P'))
        P_eva_high = self.sim.get_single('QT_inputs', 1, T_eva_air, ('P'))
        P_eva_sol = opt.brentq(cycle_dsh, P_eva_low, P_eva_high, xtol=self.tol)

        mdot, h_comp_out = self.comp.process(P_eva_sol, P_cond)
        h_cond_out, T_cond_out, m_cond = self.cond.exchange(mdot, P_cond, h_comp_out)
        h_eva_out, T_eva_out, m_eva = self.eva.exchange(mdot, P_eva_sol, h_cond_out)

        return P_eva_sol, h_comp_out, h_cond_out, T_cond_out, m_cond, h_eva_out, T_eva_out, m_eva, mdot
        
    def solve_cond(self, T_cond_air: float, T_eva_air: float):
        def cycle_charge(P_cond):
            P_eva_sol, h_comp_out, h_cond_out, T_cond_out, m_cond, h_eva_out, T_eva_out, m_eva, mdot = self.solve_evap(P_cond, T_eva_air)
            m_conn = self.conn.process(P_eva_sol, h_cond_out, h_eva_out)
            mtot = m_cond + m_eva + m_conn
            return self.charge.error(mtot)
        
        # bisect or brentq or toms748
        P_cond_low = self.sim.get_single('QT_inputs', 0, T_cond_air, ('P'))
        P_cond_high = self.sim.get_single('QT_inputs', 0, T_cond_air + 30, ('P'))

        # print(cycle_charge(P_cond_low), cycle_charge(P_cond_high))
        
        P_cond_sol = opt.brentq(cycle_charge, P_cond_low, P_cond_high, xtol=self.tol)

        P_eva_sol, h_comp_out, h_cond_out, T_cond_out, m_cond, h_eva_out, T_eva_out, m_eva, mdot = self.solve_evap(P_cond_sol, T_eva_air)
        
        return P_cond_sol, P_eva_sol
        # return P_cond_sol, P_eva_sol, h_comp_out, h_cond_out, T_cond_out, h_eva_out, T_eva_out, mdot