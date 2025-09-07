# run.py
from hx_UA_const.cycle.cycle_model import CycleModel, CycleModel_charge, CycleModel_mdot_charge
from hx_UA_const.core.params import SystemParams_DSH_DSC, SystemParams_DSH_charge, SystemParams_mdot_DSH_charge # 예시 경로, 실제 구조에 따라 조정
from dataclasses import asdict 

class SingleRunner:
    def __init__(self,
                 backend: str,
                 fluid: str,
                 base_kwargs_map: dict[str, dict] | None = None):
        self.base_kwargs_map = base_kwargs_map
        self.backend = backend
        self.fluid = fluid

    def _run_sweep(self,
                   param_cls,
                   model_cls,
                   base_kwargs: dict):
        """
        param_cls   : SystemParams_* 클래스
        model_cls   : CycleModel_* 클래스
        base_kwargs : 파라미터 (dict)
        """
        # 매번 깨끗한 리스트로 시작
        # param 생성 (기본 파라미터 + sweep 축 파라미터)
        params = param_cls(**base_kwargs)
        
        try:
            model = model_cls(params,
                            backend=self.backend,
                            fluid=self.fluid)
            result = model.calculate()
            rec = result.to_dict(key_with_unit = True)
            return model, rec
        except ValueError:
            # 실패한 케이스도 축 값은 남겨둠
            print(f"Error in calculation for {base_kwargs}")

    def run_DSH_DSC(self):
        # 1) DSH vs DSC
        model, rec= self._run_sweep(
            param_cls   = SystemParams_DSH_DSC,
            model_cls   = CycleModel,
            base_kwargs = self.base_kwargs_map["DSH_DSC"],
        )
        return model, rec

    def run_DSH_charge(self):
        # 1) DSH vs DSC
        model, rec = self._run_sweep(
            param_cls   = SystemParams_DSH_charge,
            model_cls   = CycleModel_charge,
            base_kwargs = self.base_kwargs_map["DSH_charge"],
        )
        return model, rec

    def run_mdot_charge(self):
        # 1) DSH vs DSC
        model, rec = self._run_sweep(
            param_cls   = SystemParams_mdot_DSH_charge,
            model_cls   = CycleModel_mdot_charge,
            base_kwargs = self.base_kwargs_map["mdot_charge"],
        )
        return model, rec


if __name__ == "__main__":
    # Example usage
    base_kwargs_1 = {
                "UA_total":8000, "N_cond":20, "N_eva":10,
                "T_cond_air":45+273.15, "T_eva_air":40+273.15,
                "isen_eff":0.9, "V_comp":2e-5, "f_comp":50,
                "DSH_target":5, "DSC_target":5,
                "CA":None, "tol":0.01
    }
    base_kwargs_2 = dict(
            U_cond=1000, U_eva=1000, N_cond=200, N_eva=50,
            D_cond=8e-3, L_cond=30, D_eva=6e-3, L_eva=30, L_connect=5, 
            T_cond_air=35 + 273.15, T_eva_air=27+ 273.15,
            isen_eff = 0.7, V_comp =2e-5, f_comp=50, 
            CA = None,
            DSH_target=5, charge_target=0.45, tol=0.01
        )
    base_kwargs_3 = dict(
            U_cond=1000, U_eva=1000, N_cond=200, N_eva=50,
            D_cond=8e-3, L_cond=30, D_eva=6e-3, L_eva=30, L_connect=5, 
            T_cond_air=35 + 273.15, T_eva_air=27+ 273.15,
            isen_eff = 0.7, V_comp =2e-5, f_comp=50,
            DSH_target= 0.0, # DSH_target is dummpy for this model
            CA=8e-7, charge_target = 0.32, tol=0.01
        )
    
    runner = SingleRunner(base_kwargs_map = {"DSH_DSC": base_kwargs_1, 
                           "DSH_charge": base_kwargs_2, 
                           "mdot_charge": base_kwargs_3},
                           backend = "REFPROP", 
                           fluid = "R410a")

    # Run the models
    model1, rec1 = runner.run_DSH_DSC()
    # model2, rec2 = runner.run_DSH_charge()
    # model3, rec3 = runner.run_mdot_charge()
    print(rec1)
    # print(rec2)
    # print(rec3)
    '''    
    # Plot the model results
    model1.plot()
    model2.plot()
    model3.plot()
    '''
