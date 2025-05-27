import numpy as np
import pandas as pd
import lightgbm as lgb
import time
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score

# 1) 데이터 로드 및 전처리
df = pd.read_csv('results/REFPROP_R32_SystemParams_DSH_DSC.csv')
features = ['DSH_target', 'DSC_target', 'f_comp']
df = df.dropna(subset=features + ['Q_cond [kW]', 'COP_H'])

X = df[features]
y_q    = df['Q_cond [kW]']
y_cop  = df['COP_H']

# 2) Train/Test 분할
X_train, X_test, yq_train, yq_test, ycop_train, ycop_test = train_test_split(
    X, y_q, y_cop, test_size=0.2, random_state=42
)

# 3) LightGBM Dataset 생성
dtrain_q   = lgb.Dataset(X_train, label=yq_train)
dtrain_cop = lgb.Dataset(X_train, label=ycop_train)

# 4) 하이퍼파라미터 설정
params = {
    'objective': 'regression',
    'metric': 'rmse',
    'learning_rate': 0.05,
    'num_leaves': 31,
    'verbose': -1,
    'seed': 42
}

# 5) 모델 학습
model_q   = lgb.train(params, dtrain_q,   num_boost_round=150)
model_cop = lgb.train(params, dtrain_cop, num_boost_round=150)

# 4) 단일 포인트 추론 시간 측정
sample_point = X_test.iloc[[0]]  # DataFrame 형태로 전달
n_runs = 1000

# Q_cond
start = time.perf_counter()
for _ in range(n_runs):
    _ = model_q.predict(sample_point)
end = time.perf_counter()
print(f"Q_cond 예측: {(end - start)/n_runs * 1e3:.3f} ms/점")

# COP_H
start = time.perf_counter()
for _ in range(n_runs):
    _ = model_cop.predict(sample_point)
end = time.perf_counter()
print(f"COP_H 예측:  {(end - start)/n_runs * 1e3:.3f} ms/점")