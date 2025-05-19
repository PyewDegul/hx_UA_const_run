# batch_run.py
import os
import numpy as np
import pandas as pd
from itertools import product
from hx_UA_const.cycle.cycle_model import CycleModel, CycleModel_charge ,CycleModel_mdot_charge
from hx_UA_const.core.params import SystemParams_DSH_DSC, SystemParams_DSH_charge, SystemParams_mdot_DSH_charge  # 예시

from time import time

class BatchRunner:
    def __init__(self,
                 sweeps,
                 backend: str, fluid: str,
                 res_dir: str = None):
        
        self.x1_vals = sweeps["x1_vals"]
        self.x2_vals = sweeps["x2_vals"]
        self.base_kwargs_map = sweeps["base_kwargs_map"] or {
            "DSH_DSC": {},
            "DSH_charge": {},
            "mdot_charge": {}
            }
        self.backend = backend
        self.fluid = fluid

        # 결과 디렉토리
        script_dir = os.path.dirname(os.path.abspath(__file__))
        self.res_dir = res_dir or os.path.join(script_dir, "results")
        os.makedirs(self.res_dir, exist_ok=True)

    # 아직 batch의 rec을 csv로 저장하는 것까지 있음
    # 데이터 분석 및 특정 히트맵 찾기 / 최적 모델 plot은 없음
    def _run_sweep(self,
                   param_cls,
                   model_cls,
                   sweep_axes: dict[str, np.ndarray],
                   base_kwargs: dict,
                   filename: str):
        """
        param_cls   : SystemParams_* 클래스
        model_cls   : CycleModel_* 클래스
        sweep_axes  : {'DSH': array1, 'DSC': array2} 같은 모양
        base_kwargs : sweep 축 외 공통 파라미터 (dict)
        filename    : 결과 CSV 파일명
        """
        # 매번 깨끗한 리스트로 시작
        records: list[dict] = []

        axis_names = list(sweep_axes.keys())
        axis_vals  = list(sweep_axes.values())

        for vals in product(*axis_vals):
            # ex) {'DSH': 3.2, 'DSC': 5.8}
            sweep_kwargs = dict(zip(axis_names, vals))
            # param 생성 (기본 파라미터 + sweep 축 파라미터)
            params = param_cls(**base_kwargs, **sweep_kwargs)

            # 모델 생성
            try:
                model = model_cls(params,
                              backend=self.backend,
                              fluid=self.fluid)
                result = model.calculate()
                rec = {**sweep_kwargs}
                rec.update(result.to_dict(key_with_unit=True))
            except ValueError:
                # 실패한 케이스도 축 값은 남겨둠
                rec = {**sweep_kwargs}
            records.append(rec)

        # DataFrame → CSV
        df = pd.DataFrame(records)
        out_path = os.path.join(self.res_dir, filename)
        df.to_csv(out_path, index=False)
        print(f"완료: {filename}  ({len(records)} cases)")

    def run_DSH_DSC(self):
        # 1) DSH vs DSC
        self._run_sweep(
            param_cls   = SystemParams_DSH_DSC,
            model_cls   = CycleModel,
            sweep_axes  = {'DSH_target': self.x1_vals, 'DSC_target': self.x2_vals},
            base_kwargs = self.base_kwargs_map["DSH_DSC"],
            filename    = f"{self.backend}_{self.fluid}_{SystemParams_DSH_DSC.__name__}.csv"
        )
        
    def run_DSH_charge(self):
        # 2) DSH vs Charge
        self._run_sweep(
            param_cls   = SystemParams_DSH_charge,
            model_cls   = CycleModel_charge,
            sweep_axes  = {'DSH_target': self.x1_vals, 'charge_target': self.x2_vals},
            base_kwargs = self.base_kwargs_map["DSH_charge"],
            filename    = f"{self.backend}_{self.fluid}_{SystemParams_DSH_charge.__name__}.csv"
        )

    def run_mdot_charge(self):
        # 3) CA vs Charge
        self._run_sweep(
            param_cls   = SystemParams_mdot_DSH_charge,
            model_cls   = CycleModel_mdot_charge,
            sweep_axes  = {'CA': self.x1_vals, 'charge_target': self.x2_vals},
            base_kwargs = self.base_kwargs_map["mdot_charge"],
            filename    = f"{self.backend}_{self.fluid}_{SystemParams_mdot_DSH_charge.__name__}.csv"
        )

if __name__ == "__main__":
    '''x1, x2 축의 개수'''
    x1_N = 20
    x2_N = 20
    '''기본 parameter'''
    # DSH_DSC
    base_kwarg_1 = dict(
                UA_total=1000, N_cond=200, N_eva=50,
                T_cond_air=35+273.15, T_eva_air=27+273.15,
                isen_eff=0.9, V_comp=2e-5, f_comp=50,
                CA = None, tol=0.01
    )
    # DSH_charge
    base_kwarg_2 = dict(
                U_cond=1000, U_eva=1000, N_cond=200, N_eva=50,
                D_cond=8e-3, L_cond=30, D_eva=6e-3, L_eva=30, L_connect=5,
                T_cond_air=35+273.15, T_eva_air=27+273.15,
                isen_eff=0.7, V_comp=2e-5, f_comp=50,
                CA=None, tol=0.01
    )
    # CA_charge
    base_kwarg_3 = dict(
                U_cond=1000, U_eva=1000, N_cond=200, N_eva=50,
                D_cond=8e-3, L_cond=30, D_eva=6e-3, L_eva=30, L_connect=5,
                T_cond_air=35+273.15, T_eva_air=27+273.15,
                isen_eff=0.7, V_comp=2e-5, f_comp=50,
                DSH_target=0.0, tol=0.01
    )
    # DSH_DSC
    sweeps_1 = {
        "x1_vals": np.linspace(1, 10, x1_N),
        "x2_vals": np.linspace(1, 10, x2_N),
        "base_kwargs_map": {
            "DSH_DSC": base_kwarg_1,
            "DSH_charge": {},
            "mdot_charge": {}
        }
    }
    # DSH_charge
    sweeps_2 = {
        "x1_vals": np.linspace(1, 10, x1_N),
        "x2_vals": np.linspace(0.1, 1, x2_N),
        "base_kwargs_map": {
            "DSH_DSC": {},
            "DSH_charge": base_kwarg_2,
            "mdot_charge": {}
        }
    }
    # CA_charge
    sweeps_3 = {
        "x1_vals": np.linspace(8e-7, 1e-6, x1_N),
        "x2_vals": np.linspace(0.2, 0.5, x2_N),
        "base_kwargs_map": {
            "DSH_DSC": {},
            "DSH_charge": {},
            "mdot_charge": base_kwarg_3
        }
    }


    '''2 loop - DSH/DSC'''
    run_batch = BatchRunner(
        sweeps_1,
        backend = "BICUBIC&HEOS",
        fluid   = "R32"
    )
    run_batch.run_DSH_DSC()
    
    '''2 loop - DSH/Charge'''
    run_batch = BatchRunner(
        sweeps_2,
        backend = "BICUBIC&HEOS",
        fluid   = "R32"
    )
    run_batch.run_DSH_charge()

    '''3 loop - CA/Charge'''
    run_batch = BatchRunner(
        sweeps_3,
        backend = "BICUBIC&HEOS",
        fluid   = "R32"
    )
    run_batch.run_mdot_charge()
