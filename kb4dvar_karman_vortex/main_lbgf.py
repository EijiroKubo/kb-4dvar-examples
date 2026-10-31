import os
import yaml
import numpy as np
import pandas as pd
from scipy.optimize import minimize

import calc_ns_field
import calc_edmd_adj as EDMD_ADJ

import obsdata_controller

import calc_dct as DCT2


def read_physic_model():
    """osに応じたsharedobjectを指定して読み込む"""
    if os.name == "nt":
        f = os.path.join("ns_model", "calc_forward_ns.dll")
    else:
        f = os.path.join("ns_model", "calc_forward_ns.so")

    module = calc_ns_field.NSFWD(f)

    return module

def read_config():
    """計算条件はCONFIG.ymlファイルから読みとる"""

    with open("set_config.yml", "r", encoding="utf-8") as yml:
        config = yaml.safe_load(yml)
    return config

def pool_state_to_edmd_snapshots(pool) -> np.ndarray:
    """Fortran側の (x, y, var, t) を EDMD側の (var, t, y, x) に揃える"""

    state = pool.qin_time_sr[0 : pool.imx1, 0 : pool.jmx1, 0:N_VARIABLE, :]
    return np.transpose(state, (2, 3, 1, 0))


class ForwardCalc:
    def __init__(self):

        # 設定ファイルの読み込み
        self.config = read_config()

        # クラス継承
        self.pool = read_physic_model()

        # 初期化
        self.pool.read_initial_condition()

        imax = np.empty([], dtype=np.int32)
        imx1 = np.empty([], dtype=np.int32)
        jmax = np.empty([], dtype=np.int32)
        jmx1 = np.empty([], dtype=np.int32)
        t_no = np.empty([], dtype=np.int32)

        self.pool.return_mesh_num(imax, imx1, jmax, jmx1, t_no)

        self.pool.initialize_condition(t_no)

    def setting_init_state(self, p1_external: np.ndarray = None):
        """初期浸水深を直接いじる"""

        if p1_external is None:
            # ここが同化前の初期計算条件となる。
            self.pool.set_validation_initial_state()  # 検証用の初期状態の設定

        else:
            # 初期条件
            # EDMD側は (var, y, x)、Fortran側は (x, y) で受ける。
            p1_real = np.real(p1_external)
            u0 = p1_real[0].T[0 : self.pool.imx1, 0 : self.pool.jmax]
            v0 = p1_real[1].T[0 : self.pool.imax, 0 : self.pool.jmx1]
            p0 = p1_real[2].T[0 : self.pool.imax, 0 : self.pool.jmax]
            self.pool.set_init_state(u0, v0, p0)

    def calculation(self):
        self.pool.proceed_forward()




def exec_simulation(save_anim: bool = False):

    # 物理計算による時間発展計算クラス
    forward = ForwardCalc()

    # ex-dmdの設定
    edmd = EDMD_ADJ.EDMD(
        forward.pool.t_no, forward.pool.imx1, forward.pool.jmx1, N_VARIABLE
    )

    # 離散コサイン変換
    # dct = DCT2.DCT(forward.pool.jmx1, forward.pool.imx1)
    # elim_ny = int((1 - CUT_RATE) * forward.pool.imx1)
    # elim_nx = int((1 - CUT_RATE) * forward.pool.jmx1)

    # 観測データ
    obs = obsdata_controller.ObsDataProcessor()
    obs.read_observation(OBS_SOURCES)

    # 観測誤差分散行列 -> 一旦単位行列で設定
    Rk = np.eye(obs.obs_num)

    # 更新対象とする初期条件
    p1_temp = None

    # 以下、同化計算
    for t_i in range(int(CONFIG["DA"]["global_loops"])):
        # 物理計算
        forward.setting_init_state(p1_temp)
        forward.calculation()

        calc_dat = pool_state_to_edmd_snapshots(forward.pool)
        edmd.exec_mode_decomposition_from_snapshops(calc_dat)

        # 随伴コード
        edmd.prepare_adj()


        # 擬似随伴コードによる評価関数・勾配計算
        # 1回目: 従来どおり固定学習率 alpha の勾配法
        # 2回目以降: SciPy の L-BFGS-B（履歴数 maxcor=20）
        local_loops = int(CONFIG["DA"]["local_loops"])
        assimilate_steps = int(CONFIG["DA"]["assimilate_steps"])
        alpha = float(CONFIG["DA"]["alpha"])

        obs_matrix_initialized = False
        evaluation_cache = {}

        def objective_and_gradient(x):
            """EDMD予測から評価関数 J と随伴勾配 dJ/dx を返す。"""
            nonlocal obs_matrix_initialized, evaluation_cache

            # SciPyから渡される最適化変数をEDMD初期値へ反映
            edmd.p1_init = np.asarray(x, dtype=np.float64).copy()

            # EDMDによる時間発展
            VV = np.asarray(edmd.exec_dmd())
            obs_calc_mat = obs.obs_function_for_adj(VV)

            if not obs_matrix_initialized:
                nv, nt, ny, nx = np.shape(VV)
                obs.set_obsmatrix_for_adj(nv, ny, nx, nt, OBS_STEP)
                obs_matrix_initialized = True

            # 評価関数 J = 1/2 * sum((H(x)-y)^2)
            # 現行コードと同様、obs_mask を掛けた観測誤差を用いる。
            func_j = 0.0
            p_t_prev = np.zeros_like(edmd.p1_init)

            for t_k in range(assimilate_steps):
                k_step = assimilate_steps - (t_k + 1)

                # アドジョイントコードの実行
                p_t_temp = edmd.adjoint_code(p_t_prev)

                error_temp = (
                    obs_calc_mat[:, :, :, k_step]
                    - obs.obs_mat[:, :, :, k_step]
                ) * obs.obs_mask[:, :, :, k_step]

                error_real = np.asarray(error_temp.real, dtype=np.float64)
                func_j += 0.5 * float(np.dot(error_real.ravel(), error_real.ravel()))

                # 観測誤差をリフト空間と同じサイズへ埋め込む
                error_mat = np.zeros((np.size(p_t_prev), 1), dtype=np.float64)
                state_size = edmd.nvar * edmd.ny * edmd.nx
                error_mat[:state_size, 0] = error_real.reshape(state_size)

                p_t_prev = p_t_temp + error_mat

            gradient = np.asarray(
                np.real(p_t_prev).reshape(-1), dtype=np.float64
            )

            if gradient.size != edmd.p1_init.size:
                raise ValueError(
                    "Gradient size does not match edmd.p1_init: "
                    f"{gradient.size} != {edmd.p1_init.size}"
                )

            evaluation_cache = {
                "x": edmd.p1_init.copy(),
                "fun": func_j,
                "grad": gradient.copy(),
                "VV": VV,
                "obs_calc_mat": obs_calc_mat,
            }
            return func_j, gradient

        if local_loops > 0:
            # --------------------------------------------------------------
            # Step 1: 従来の固定学習率による勾配法
            # --------------------------------------------------------------
            x0 = np.asarray(np.real(edmd.p1_init), dtype=np.float64).copy()
            f0, g0 = objective_and_gradient(x0)
            x1 = x0 - alpha * g0
            edmd.p1_init = x1.copy()

            print(
                "gradient step:",
                "global=", t_i,
                "J=", f0,
                "|g|=", np.linalg.norm(g0),
            )

            # --------------------------------------------------------------
            # Step 2以降: L-BFGS-B
            # maxcor=20 が保存する (s, y) 履歴数に対応する。
            # --------------------------------------------------------------
            if local_loops > 1:
                lbfgs_iteration = 0

                def lbfgs_callback(xk):
                    nonlocal lbfgs_iteration
                    lbfgs_iteration += 1

                    print(
                        "L-BFGS-B step:",
                        "global=", t_i,
                        "local=", lbfgs_iteration + 1,
                        "J=", evaluation_cache.get("fun", np.nan),
                        "|g|=", np.linalg.norm(
                            evaluation_cache.get("grad", np.array([np.nan]))
                        ),
                    )

                result = minimize(
                    objective_and_gradient,
                    x1,
                    method="L-BFGS-B",
                    jac=True,
                    callback=lbfgs_callback,
                    options={
                        "maxiter": local_loops - 1,
                        "maxcor": 20,
                        "maxls": 20,
                        "ftol": 0.4 , # 1.0e-12,  # ここで、発散を制限
                        "gtol": 1.0e-6,
                    },
                )
                edmd.p1_init = np.asarray(result.x, dtype=np.float64).copy()


            # 最終状態で評価し直し、診断出力との整合を取る
            final_j, final_g = objective_and_gradient(edmd.p1_init)
            VV = evaluation_cache["VV"]
            obs_calc_mat = evaluation_cache["obs_calc_mat"]
            t_j = local_loops - 1

            if save_anim == True:
                header_id_tmp = int((t_i + 1) * 100) + int(local_loops)
                edmd.show_result(VV, header_id_tmp)

        else:
            t_j = -1

        p1_temp = np.reshape(
            edmd.p1_init[: forward.pool.imx1 * forward.pool.jmx1 * edmd.nvar],
            (edmd.nvar, forward.pool.jmx1, forward.pool.imx1),
        )

        # グローバルループ終了
        print("end g-loop...", "iter=" + str(t_i))

    # 補正後予測計算
    forward.setting_init_state(p1_temp)
    forward.calculation()
    calc_dat = pool_state_to_edmd_snapshots(forward.pool)


    if save_anim == True:
        # アニメーション表示用データ
        edmd.show_result(calc_dat, 10001)

        # edmdでの計算 -> 最後の計算結果
        edmd.show_result(VV, 10000)


CONFIG = read_config()

# 計算時間の計測時はFalseにする
SAVE_ANIM = False


N_VARIABLE = 3
OBS_STEP = int(CONFIG["OBS"]["step"])
OBS_SOURCES = [
    (0, "observe_u.zip"),  # u
    (1, "observe_v.zip"),  # v
]


if __name__ == "__main__":
    # カレントディレクトリの位置を変更
    os.chdir(os.path.dirname(__file__))

    # 計算
    exec_simulation(SAVE_ANIM)
