# core/sim_cycle.py
import CoolProp.CoolProp as CP
from operator import itemgetter, attrgetter
import numpy as np
from scipy.interpolate import RectBivariateSpline
import os
import cProfile

class SimCycle:
    # Backend 설정(with Coolprop wrapper), ("REFPROP" or "HEOS", ...)
    # Fluid 설정 (ex. "R32", "R134a", "R410A", ...)
    # REFPROP 백엔드 : https://www.coolprop.org/_static/doxygen/html/class_cool_prop_1_1_r_e_f_p_r_o_p_backend.html
    # 필요하면 여기서 가져올 수 있음 : 
    def __init__(self, backend_name: str, fluid_name: str):
        self.state = CP.AbstractState(backend_name, fluid_name)
        self.P_C = self.state.p_critical()
        self.T_C = self.state.T_critical()
        self.S_C = self.state.trivial_keyed_output(CP.iP_triple)

        self.P_TP = self.state.trivial_keyed_output(CP.iP_triple)
        self.T_TP = self.state.Ttriple()

        # input_keys: CoolProp에서 사용하는 입력값들
        # S, H, D는 모두 질량 기준
        # input_keys : https://www.coolprop.org/_static/doxygen/html/namespace_cool_prop.html#a58e7d98861406dedb48e07f551a61efb
        self.input_keys = {
            'QT_inputs': CP.QT_INPUTS,
            'PQ_inputs': CP.PQ_INPUTS,
            'PT_inputs': CP.PT_INPUTS,
            'PS_inputs': CP.PSmass_INPUTS,
            'HQ_inputs': CP.HmassQ_INPUTS, 	
            'HP_inputs': CP.HmassP_INPUTS,
            'HT_inputs': CP.HmassT_INPUTS,   
        }
        # output_pairs: CoolProp에서 사용하는 출력값들
        # output_pairs : https://www.coolprop.org/_static/doxygen/html/namespace_cool_prop.html#a4b49eeb37210a720b188f493955d8364
        self.output_pairs = {
            'T': CP.iT,
            'P': CP.iP,
            'H': CP.iHmass,
            'S': CP.iSmass,
            'D': CP.iDmass,
            'Q': CP.iQ,
            'C': CP.iCpmass,
        }
        # phase_pairs: CoolProp에서 사용하는 유체의 상
        self.phase_pairs = {
            'L': CP.iphase_liquid,
            'G': CP.iphase_gas,
            'S': CP.iphase_supercritical,
            'TP': CP.iphase_twophase,
            'None': CP.iphase_not_imposed
        }
    
    def get_state(self):
        return self.state
    
    # Specify the phase of the fluid
    def specify_phase(self, phase):
        self.state.specify_phase(self.phase_pairs[phase])

    # Update the state with the given inputs
    def update(self, arg, input1, input2):
        self.state.update(self.input_keys[arg], input1, input2)

    
    # Get the value of a single property(float) or multiple properties(ndarray)
    def get_single_no_update(self, prop : str)-> dict:
        return self.state.keyed_output(self.output_pairs[prop])

    def get_multiple_no_update(self, props : tuple):
        vals = itemgetter(*props) (self.output_pairs)
        return tuple(map(self.state.keyed_output, vals))

    def get_single(self, arg, input1, input2, props):
        self.update(arg, input1, input2)
        return self.get_single_no_update(props)
    
    def get_multiple(self, arg, input1, input2, props):
        self.update(arg, input1, input2)
        return self.get_multiple_no_update(props)


class SimCycle_Interp:
    def __init__(self, backend_name: str, fluid_name: str, 
                 n_points_p: int = 200, n_points_h: int = 200):
        """압력-엔탈피 보간 테이블 초기화 (CoolProp REFPROP 사용)"""
        self.backend = backend_name
        self.fluid = fluid_name
        self.state = CP.AbstractState(backend_name, fluid_name)
        # 임계압력 및 삼중점 (필요 시 사용)
        self.P_c = self.state.p_critical()
        self.T_c = self.state.T_critical()
        self.P_triple = self.state.trivial_keyed_output(CP.iP_triple)
        self.T_triple = self.state.Ttriple()
        
        # 압력 및 엔탈피 범위 설정
        P_min = 1e5  # 0.1 MPa 
        if P_min < self.P_triple:
            P_min = self.P_triple * 1.01  # 삼중점보다 약간 높게 (안전상)
        P_max = 0.95 * self.P_c
        # 포화 엔탈피 한계 (최소 압력에서)
        self.state.update(CP.PQ_INPUTS, P_min, 0)   # 포화 액체
        h_liq_min = self.state.keyed_output(CP.iHmass)
        self.state.update(CP.PQ_INPUTS, P_min, 1)   # 포화 증기
        h_vap_min = self.state.keyed_output(CP.iHmass)
        h_min, h_max = h_liq_min, h_vap_min
        
        # 1) 저장할 디렉터리와 파일명 정의
        script_dir = os.path.dirname(os.path.abspath(__file__))
        cache_dir = os.path.join(script_dir, "cache")
        os.makedirs(cache_dir, exist_ok=True)               # 폴더 없으면 생성
        cache_filename = os.path.join(cache_dir, f"{self.backend}_{self.fluid}_{n_points_p}_{n_points_h}_PH_table.npz")
        
        try:
            # 캐시된 테이블 불러오기 시도
            data = np.load(cache_filename)
            P_vals = data['P_vals']; H_vals = data['H_vals']
            T_table = data['T_table']; S_table = data['S_table']
            D_table = data['D_table']; Q_table = data['Q_table']
            C_table = data['C_table']  # cp table
        except FileNotFoundError:
            # 테이블 생성 수행
            P_vals = np.linspace(P_min, P_max, n_points_p)
            H_vals = np.linspace(h_min, h_max, n_points_h)
            # 테이블 메모리 할당
            T_table = np.empty((n_points_p, n_points_h))
            S_table = np.empty((n_points_p, n_points_h))
            D_table = np.empty((n_points_p, n_points_h))
            Q_table = np.empty((n_points_p, n_points_h))
            C_table = np.empty((n_points_p, n_points_h))
            # 각 압력 단계별로 루프
            for i, P in enumerate(P_vals):
                # 현재 압력에서 포화 엔탈피 (액/증기) 계산
                if P < self.P_c:
                    self.state.update(CP.PQ_INPUTS, P, 0)
                    h_liq = self.state.keyed_output(CP.iHmass)
                    self.state.update(CP.PQ_INPUTS, P, 1)
                    h_vap = self.state.keyed_output(CP.iHmass)
                else:
                    # 임계압력 이상에서는 포화 개념 없음
                    h_liq = h_vap = None
                # 엔탈피 방향으로 루프
                for j, H in enumerate(H_vals):
                    self.state.update(CP.HmassP_INPUTS, H, P)
                    # CoolProp에서 필요한 출력 추출
                    T_val = self.state.keyed_output(CP.iT)
                    S_val = self.state.keyed_output(CP.iSmass)
                    D_val = self.state.keyed_output(CP.iDmass)
                    Q_val = self.state.keyed_output(CP.iQ)
                    C_val = self.state.keyed_output(CP.iCpmass)
                    # 혹시 Q가 NaN이거나 (<0, >1) 범위를 벗어나면 보정
                    if not np.isfinite(Q_val):
                        # 두상 영역 밖: 포화선 이내 값을 벗어나면 Q를 0 또는 1로 설정
                        if h_liq is not None and H < h_liq:
                            Q_val = 0.0
                        elif h_vap is not None and H > h_vap:
                            Q_val = 1.0
                        else:
                            # 그 외 예외적인 경우는 0으로
                            Q_val = 0.0
                    else:
                        # Q 값이 부정확하게 <0 또는 >1인 경우 클램핑
                        if Q_val < 0.0: Q_val = 0.0
                        if Q_val > 1.0: Q_val = 1.0
                    # cp 값이 비정상 (예: 두상 혼합에서 매우 큼)인 경우 NaN 처리
                    if not np.isfinite(C_val):
                        C_val = np.nan
                    T_table[i, j] = T_val
                    S_table[i, j] = S_val
                    D_table[i, j] = D_val
                    Q_table[i, j] = Q_val
                    C_table[i, j] = C_val
            # 계산 완료 후 캐시 파일로 저장
            np.savez_compressed(cache_filename,
                                 P_vals=P_vals, H_vals=H_vals,
                                 T_table=T_table, S_table=S_table,
                                 D_table=D_table, Q_table=Q_table,
                                 C_table=C_table)
        # 2D 보간 스플라인 생성 (건도 Q는 1차, 나머지 3차 보간)
        self.T_spline   = RectBivariateSpline(P_vals, H_vals, T_table, kx=3, ky=3, s=0)
        self.S_spline   = RectBivariateSpline(P_vals, H_vals, S_table, kx=3, ky=3, s=0)
        self.D_spline   = RectBivariateSpline(P_vals, H_vals, D_table, kx=3, ky=3, s=0)
        self.Q_spline   = RectBivariateSpline(P_vals, H_vals, Q_table, kx=1, ky=1, s=0)
        # cp 테이블에 NaN이 있을 경우 보간기에서 자동 처리되지만, 
        # 필요시 NaN을 0으로 대체하여 보간할 수도 있음 (여기서는 그대로 둠)
        self.C_spline   = RectBivariateSpline(P_vals, H_vals, C_table, kx=3, ky=3, s=0)
        # 마지막 업데이트된 상태 저장용 딕셔너리
        self.last_vals = {}

    def specify_phase(self, phase: str):
        """CoolProp 상 계산 시 상(specify phase) 지정 (보간 계산에는 영향 없음)"""
        # SimCycle의 phase_pairs 맵 활용
        phase_pairs = {
            'L': CP.iphase_liquid,
            'G': CP.iphase_gas,
            'S': CP.iphase_supercritical,
            'TP': CP.iphase_twophase,
            'None': CP.iphase_not_imposed
        }
        if phase in phase_pairs:
            self.state.specify_phase(phase_pairs[phase])
        else:
            raise ValueError(f"Unknown phase specification: {phase}")

    def update(self, arg: str, input1: float, input2: float):
        """상태 업데이트 - 현재는 Hmass-P 입력만 지원"""
        if arg != 'HP_inputs':
            raise NotImplementedError("SimCycle_Interp only supports HP_inputs (enthalpy-pressure) updates.")
        # input1 = enthalpy [J/kg], input2 = pressure [Pa]
        h = input1; P = input2
        # 보간기로 값 계산
        T = float(self.T_spline.ev(P, h))
        S = float(self.S_spline.ev(P, h))
        D = float(self.D_spline.ev(P, h))
        Q = float(self.Q_spline.ev(P, h))
        C = float(self.C_spline.ev(P, h))
        # last_vals 딕셔너리에 저장
        self.last_vals = {'P': P, 'H': h, 'T': T, 'S': S, 'D': D, 'Q': Q, 'C': C}

    def get_single_no_update(self, prop: str):
        """마지막 상태에서 단일 속성 반환"""
        return self.last_vals[prop]

    def get_multiple_no_update(self, props: tuple):
        """마지막 상태에서 여러 속성 값을 튜플로 반환"""
        return tuple(self.last_vals[p] for p in props)

    def get_single(self, arg: str, input1: float, input2: float, prop: str):
        """입력 갱신 후 단일 속성 반환"""
        self.update(arg, input1, input2)
        return self.get_single_no_update(prop)

    def get_multiple(self, arg: str, input1: float, input2: float, props: tuple):
        """입력 갱신 후 다중 속성 반환"""
        self.update(arg, input1, input2)
        return self.get_multiple_no_update(props)



def test_sim():
    sim = SimCycle("REFPROP", "R32")
    N = int(1e5)

    T_elem = np.empty(N)
    s_elem = np.empty(N)
    q_elem = np.empty(N)
    rho_elem = np.empty(N)
    h_elem = np.linspace(280*1e3, 320*1e3, N)
    for i in range(N):
        T_elem[i], s_elem[i], q_elem[i], rho_elem[i] = sim.get_multiple('HP_inputs', h_elem[i], 3e6, ('T', 'S', 'Q', 'D'))

def test_sim_interp():
    sim = SimCycle_Interp("REFPROP", "R32", n_points_p =  200, n_points_h = 200)
    N = int(1e5)

    T_elem = np.empty(N)
    s_elem = np.empty(N)
    q_elem = np.empty(N)
    rho_elem = np.empty(N)
    h_elem = np.linspace(280*1e3, 320*1e3, N)
    for i in range(N):
        T_elem[i], s_elem[i], q_elem[i], rho_elem[i] = sim.get_multiple('HP_inputs', h_elem[i], 3e6, ('T', 'S', 'Q', 'D'))

cProfile.run('test_sim()')
cProfile.run('test_sim_interp()')
