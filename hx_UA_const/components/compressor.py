class Compressor:
    def __init__(self, sim, params):
        self.sim = sim
        self.p = params
        self.DSH = params.DSH_target

    def set_DSH(self, target: float):
        self.DSH = target

    def process(self, P_eva: float, P_cond:float) -> tuple[float, float]:
        # Initialize the simulation cycle(압축기는 inlet은 기체! Coolprop 헷갈리지 않게 하기)
        self.sim.specify_phase('G')
        
        # 압축기 입구 온도 및 DSH 계산산
        T_sat = self.sim.get_single('PQ_inputs', P_eva, 1, ('T'))
        T_comp_in = T_sat + self.DSH

        # 압축기 입구 엔탈피, 엔트로피, 밀도, 유량의 계산 ()
        h_comp_in, s_comp_in, rho_comp_in = self.sim.get_multiple('PT_inputs', P_eva, T_comp_in, ('H', 'S', 'D'))
        mdot = rho_comp_in * self.p.V_comp * self.p.f_comp

        # 압축기 출구 엔탈피, 유량 계산
        h_comp_out_iso = self.sim.get_single('PS_inputs', P_cond, s_comp_in, ('H'))
        h_comp_out = h_comp_in + (h_comp_out_iso - h_comp_in) / self.p.isen_eff
        s_comp_out, T_comp_out = self.sim.get_multiple('HP_inputs', h_comp_out, P_cond, ('S', 'T'))

        self.sim.specify_phase('None')

        return h_comp_out, s_comp_out, T_comp_out, mdot
