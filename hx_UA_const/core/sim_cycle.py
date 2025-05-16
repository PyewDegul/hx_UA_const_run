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
        self.backend_name = backend_name
        self.fluid_name = fluid_name
        # critical point
        self.P_C = self.state.p_critical()
        self.T_C = self.state.T_critical()
        
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


class SimCycle_Interp(SimCycle):
    def __init__(self, backend_name: str, fluid_name: str, 
                 n_points_p: int = 200, n_points_h: int = 200, 
                 degree_x: int = 3, degree_y: int = 3):
        super().__init__(backend_name, fluid_name)
        T_max = self.state.Tmax()

        # 1) 저장할 디렉터리와 파일명 정의
        script_dir = os.path.dirname(os.path.abspath(__file__))
        cache_dir = os.path.join(script_dir, "cache")
        os.makedirs(cache_dir, exist_ok=True)               # 폴더 없으면 생성
        cache_filename = os.path.join(cache_dir, f"{self.backend_name}_{self.fluid_name}_n_ph_{n_points_p}&{n_points_h}_deg_{degree_x}&{degree_y}_PH_table.npz")

        # 압력 및 엔탈피 범위 설정        
        # 포화 엔탈피 한계 (최소 압력에서)
        """압력-엔탈피 보간 테이블 초기화 (CoolProp REFPROP 사용)"""
        P_min = 0.1*1e6  # 삼중점보다 약간 높게 (안전상)
        P_max = 0.95 * self.P_C  # 최대 압력 (최대 압력)
        h_min = self.get_single('PQ_inputs', P_min, 0, 'H')
        h_max = self.get_single('PT_inputs', P_min, T_max, 'H')

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
            
            # 1) n_points_p×n_points_h×5 배열(raw) 구축
            raw = np.array([
                [self.get_multiple('HP_inputs', H, P, ('T','S','D','Q','C'))
                for H in H_vals ]
                for P in P_vals
            ], dtype=np.float64)   # shape == (nP, nH, 5)

            # 2) C (cp) 컬럼의 NaN 처리 (in-place)
            C = raw[..., 4]
            C[~np.isfinite(C)] = np.nan

            # 3) Q (quality) 컬럼 클램핑 + NaN 보정
            Q = raw[..., 3]
            # 3-1) 0–1 범위 clamp
            np.clip(Q, 0.0, 1.0, out=Q)

            # h_liq/vap 벡터화
            h_liq_arr = np.full_like(P_vals, np.nan)
            h_vap_arr = np.full_like(P_vals, np.nan)
            # h_liq, h_vap 도 벡터로 계산 (P_vals 별)
            mask_P = P_vals < self.P_C
            h_liq_arr[mask_P] = np.array(
                [self.get_single('PQ_inputs', p, 0, 'H') for p in P_vals[mask_P]],
                dtype=np.float64
            )
            h_vap_arr[mask_P] = np.array(
                [self.get_single('PQ_inputs', p, 1, 'H') for p in P_vals[mask_P]],
                dtype=np.float64
            )
            
            # 3-2) broadcast 마스크
            mask = ~np.isfinite(Q)
            mask_liq = (H_vals[None, :] < h_liq_arr[:, None]) & mask
            mask_vap = (H_vals[None, :] > h_vap_arr[:, None]) & mask
            Q[mask_liq] = 0.0;  Q[mask_vap] = 1.0;  Q[mask & ~(mask_liq|mask_vap)] = 0.0
            
            # 4) 최종 T, S, D, Q, C 테이블로 언패킹
            T_table, S_table, D_table, Q_table, C_table = (
                raw[..., k] for k in range(5)
            )   

            # 계산 완료 후 캐시 파일로 저장
            np.savez_compressed(cache_filename,
                                 P_vals=P_vals, H_vals=H_vals,
                                 T_table=T_table, S_table=S_table,
                                 D_table=D_table, Q_table=Q_table,
                                 C_table=C_table)
            
        # 2D 보간 스플라인 생성 (건도 Q는 1차, 나머지 3차 보간)
        self.T_spline   = RectBivariateSpline(P_vals, H_vals, T_table, kx=degree_x, ky=degree_y, s=0)
        self.S_spline   = RectBivariateSpline(P_vals, H_vals, S_table, kx=degree_x, ky=degree_y, s=0)
        self.D_spline   = RectBivariateSpline(P_vals, H_vals, D_table, kx=degree_x, ky=degree_y, s=0)
        self.Q_spline   = RectBivariateSpline(P_vals, H_vals, Q_table, kx=1, ky=1, s=0)
        # cp 테이블에 NaN이 있을 경우 보간기에서 자동 처리되지만, 
        # 필요시 NaN을 0으로 대체하여 보간할 수도 있음 (여기서는 그대로 둠)
        self.C_spline   = RectBivariateSpline(P_vals, H_vals, C_table, kx=degree_x, ky=degree_y, s=0)

        # --- in SimCycle_Interp.__init__() 끝부분에 -------------------
        self._idx = {'P': 0, 'H': 1, 'T': 2, 'S': 3,
                    'D': 4, 'Q': 5, 'C': 6}
        # 7개의 최근 상태 값을 담을 1-D ndarray
        self.last_vals = np.empty(7, dtype=np.float64)

    def specify_phase_itp(self, phase: str):
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

    def update_itp(self, arg: str, input1: float, input2: float):
        """상태 업데이트 - 현재는 Hmass-P 입력만 지원"""
        if arg != 'HP_inputs':
            raise NotImplementedError("SimCycle_Interp only supports HP_inputs (enthalpy-pressure) updates.")
        # input1 = enthalpy [J/kg], input2 = pressure [Pa]
        h = input1; P = input2
        # spline은 여전히 1 스칼라씩 계산해야 하지만
        # 결과를 바로 ndarray에 저장한다
        self.last_vals[:] = (
            P, h,
            self.T_spline.ev(P, h),
            self.S_spline.ev(P, h),
            self.D_spline.ev(P, h),
            self.Q_spline.ev(P, h),
            self.C_spline.ev(P, h)
        )

    def get_single_no_update_itp(self, prop: str):
        """마지막 상태에서 단일 속성 반환"""
        return float(self.last_vals[self._idx[prop]])

    def get_multiple_no_update_itp(self, props: tuple):
        """마지막 상태에서 여러 속성 값을 튜플로 반환"""
        idxs = itemgetter(*props) (self._idx)
        return tuple(self.last_vals.take(idxs))

    def get_single_itp(self, arg: str, input1: float, input2: float, prop: str):
        """입력 갱신 후 단일 속성 반환"""
        self.update_itp(arg, input1, input2)
        return self.get_single_no_update_itp(prop)

    def get_multiple_itp(self, arg: str, input1: float, input2: float, props: tuple):
        """입력 갱신 후 다중 속성 반환"""
        self.update_itp(arg, input1, input2)
        return self.get_multiple_no_update_itp(props)



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
    sim = SimCycle_Interp("REFPROP", "R32", n_points_p =  1000, n_points_h = 1000, degree_x = 1, degree_y = 1)
    N = int(1e5)

    T_elem = np.empty(N)
    s_elem = np.empty(N)
    q_elem = np.empty(N)
    rho_elem = np.empty(N)
    h_elem = np.linspace(280*1e3, 320*1e3, N)
    for i in range(N):
        T_elem[i], s_elem[i], q_elem[i], rho_elem[i] = sim.get_multiple_itp('HP_inputs', h_elem[i], 3e6, ('T', 'S', 'Q', 'D'))

test_sim_interp()
'''
cProfile.run('test_sim()')
cProfile.run('test_sim_interp()')
'''