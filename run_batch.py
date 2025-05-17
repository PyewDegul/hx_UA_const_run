# batch_run.py
import os
import numpy as np
import pandas as pd
from itertools import product
from hx_UA_const.cycle.cycle_model import CycleModel, CycleModel_charge ,CycleModel_mdot_charge
from hx_UA_const.core.params import SystemParams_DSH_DSC, SystemParams_DSH_charge, SystemParams_mdot_DSH_charge  # 예시
from dataclasses import asdict

x1_N = 100
x2_N = 100

'''1. DSH, DSC to Pressure Solver(2-loop)'''
# 1) 파라미터 범위
DSH_vals = np.linspace(1, 10, x1_N)
DSC_vals = np.linspace(1, 10, x2_N) 

# 2) 결과 수집
records = []
for dsh, dsc in product(DSH_vals, DSC_vals):
    params = SystemParams_DSH_DSC(
        # Cond, Evap UA, N 개별 설정
        UA_total=1000, N_cond=200, N_eva=50,
        # 이차 유체 물성 설정
        T_cond_air=35 + 273.15, T_eva_air=27+ 273.15,
        # 압축기 작동 설정
        isen_eff = 0.9, V_comp =2e-5, f_comp=50,
        # DSH, DSC 설정
        CA = None, DSH_target=dsh, DSC_target=dsc,
        # solver tol 설정
        tol=0.01
    )
    try:
        model = CycleModel(params, backend="BICUBIC&HEOS", fluid="R32")
        res = model.calculate()  # dataclass 인스턴스
        res_dict = asdict(res)
    except ValueError:
        res_dict  = {}
    record = {"DSH[K]": dsh, "DSC[K]": dsc}
    record.update(res_dict)  # 결과 추가
    records.append(record)

# 3) 저장할 디렉터리와 파일명 정의
script_dir = os.path.dirname(os.path.abspath(__file__))
res_dir = os.path.join(script_dir, "results")
os.makedirs(res_dir, exist_ok=True)               # 폴더 없으면 생성
res_filename = os.path.join(res_dir, f"{model.backend}_{model.fluid}_N1_N2_{x1_N}&{x2_N}_DSH_DSC.csv")

# 4) DataFrame 변환 및 저장
df = pd.DataFrame(records) 
df.to_csv(res_filename, index=False)
print("배치 완료, 결과: DSH_DSC_table")

'''2. DSH, Charge to Pressure Solver(2-loop)'''

# 1) 파라미터 범위
DSH_vals = np.linspace(1, 10, x1_N)
charge_vals = np.linspace(0.1, 1, x2_N)

# 2) 결과 수집
records = []
for dsh, charge_val in product(DSH_vals, charge_vals):
    params = SystemParams_DSH_charge(
        # Cond, Evap U, N 개별 설정
        U_cond=1000, U_eva=1000, N_cond=200, N_eva=50,
        # Cond, Evap, Conn 기하 형상 설정
        D_cond=8e-3, L_cond=30, D_eva=6e-3, L_eva=30, L_connect=5,
        # 이차 유체 물성 설정
        T_cond_air=35 + 273.15, T_eva_air=27+ 273.15,
        # 압축기 작동 설정
        isen_eff = 0.7, V_comp =2e-5, f_comp=50, 
        # DSH, charge 설정
        CA = None, DSH_target=dsh, charge_target=charge_val, 
        # solver tol 설정
        tol=0.01
        )
    try:
        model = CycleModel_charge(params, backend="BICUBIC&HEOS", fluid="R32")
        res = model.calculate()  # dataclass 인스턴스
        res_dict = asdict(res)
    except ValueError:
        res_dict  = {}
    record = {"DSH[K]": dsh, "charge[kg]": charge_val}
    record.update(res_dict)  # 결과 추가
    records.append(record)

# 3) 저장할 디렉터리와 파일명 정의
script_dir = os.path.dirname(os.path.abspath(__file__))
res_dir = os.path.join(script_dir, "results")
os.makedirs(res_dir, exist_ok=True)               # 폴더 없으면 생성
res_filename = os.path.join(res_dir, f"{model.backend}_{model.fluid}_N1_N2_{x1_N}&{x2_N}_DSH_charge.csv")

# 4) DataFrame 변환 및 저장
df = pd.DataFrame(records) 
df.to_csv(res_filename, index=False)
print("배치 완료, 결과: DSH_charge_table")


'''3. mdot, DSH, Charge to Pressure Solver(3-loop) // CA, m_charge'''
# 1) 파라미터 범위
CA_vals = np.linspace(1e-8, 1e-5, x1_N)  # 예시로 추가한 CA 범위
charge_vals = np.linspace(0.1, 0.5, x2_N)  # 예시로 추가한 charge 범위

# 2) 결과 수집
records = []
for CA_val, charge_val in product(CA_vals, charge_vals):
    params = SystemParams_mdot_DSH_charge(
        # Cond, Evap U, N 개별 설정
        U_cond=1000, U_eva=1000, N_cond=200, N_eva=50,
        # Cond, Evap, Conn 기하 형상 설정
        D_cond=8e-3, L_cond=30, D_eva=6e-3, L_eva=30, L_connect=5,
        # 이차 유체 물성 설정
        T_cond_air=35 + 273.15, T_eva_air=27+ 273.15,
        # 압축기 작동 설정
        isen_eff = 0.7, V_comp =2e-5, f_comp=50,
        # CA(mdot) 개도, charge 설정
        CA = CA_val, # CA 값은 linspace를 통해 구간 설정
        DSH_target=0.0,  # DSH_target을 0으로 설정(dummy)
        charge_target=charge_val,  # charge 값은 linspace를 통해 구간 설정
        # solver tol 설정
        tol = 0.01
    )
    try:
        model = CycleModel_mdot_charge(params, backend="BICUBIC&HEOS", fluid="R32")
        res = model.calculate()  # dataclass 인스턴스
        res_dict = asdict(res)
    except ValueError:
        res_dict  = {}
    record = {"CA[m^2]": CA_val, "charge[kg]": charge_val}
    record.update(res_dict)  # 결과 추가
    records.append(record)

# 3) 저장할 디렉터리와 파일명 정의
script_dir = os.path.dirname(os.path.abspath(__file__))
res_dir = os.path.join(script_dir, "results")
os.makedirs(res_dir, exist_ok=True)               # 폴더 없으면 생성
res_filename = os.path.join(res_dir, f"{model.backend}_{model.fluid}_N1_N2_{x1_N}&{x2_N}_CA_charge_table.csv")

# 4) DataFrame 변환 및 저장
df = pd.DataFrame(records) 
df.to_csv(res_filename, index=False)
print("배치 완료, 결과: CA_charge_table")  