class ExpansionValve:
    def __init__(self, sim, params):
        self.sim = sim
        self.p = params

    def process(self, P_eva: float, P_cond:float, h_cond_out) -> tuple[float, float]:
        # Calculate the enthalpy and entropy at the evaporator inlet and outlet
        h_exp_in = h_cond_out
        rho_exp_in = self.sim.get_single('HP_inputs', h_exp_in, P_cond, 'D')
        h_exp_out = h_exp_in
        s_exp_out, T_exp_out = self.sim.get_multiple('HP_inputs', h_exp_out, P_eva, ('S', 'T'))

        if self.p.CA is not None:
            CA = self.p.CA
            mdot = CA * (2 * rho_exp_in * (P_cond - P_eva)) ** 0.5
            return h_exp_out, s_exp_out, T_exp_out, mdot
        else:
            return h_exp_out, s_exp_out, T_exp_out

        
        
        