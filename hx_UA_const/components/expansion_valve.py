class ExpansionValve:
    def __init__(self, sim, params):
        self.sim = sim
        self.p = params

    def process(self, P_eva: float, P_cond:float, h_cond_out) -> tuple[float, float]:
        # Calculate the enthalpy and entropy at the evaporator inlet and outlet
        h_exp_in = h_cond_out
        h_exp_out = h_exp_in
        s_exp_out, T_exp_out = self.sim.get_multiple('PH_inputs', P_eva, h_exp_out, ('S', 'T'))

        # Calculate the mass flow rate(Orifice area)
        # mdot = f(h_exp_in, h_exp_out, P_eva, P_cond, self.Area)
        
        return h_exp_out, s_exp_out, T_exp_out, # mdot