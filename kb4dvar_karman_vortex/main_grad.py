import os
import yaml
import numpy as np

import calc_ns_field
import calc_edmd_adj as EDMD_ADJ

import obsdata_controller


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
            self.pool.set_validation_initial_state()  

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


def exec_simulation():

    # 物理計算による時間発展計算クラス
    forward = ForwardCalc()

    # ex-dmdの設定
    edmd = EDMD_ADJ.EDMD(
        forward.pool.t_no, forward.pool.imx1, forward.pool.jmx1, N_VARIABLE
    )

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

        #　
        calc_dat = pool_state_to_edmd_snapshots(forward.pool)
        edmd.exec_mode_decomposition_from_snapshops(calc_dat)

        # 随伴コード
        edmd.prepare_adj()


        # 擬似随伴コードによる計算
        for t_j in range(int(CONFIG["DA"]["local_loops"])):
            # print(t_j, "...start-localloop")

            # dmdでのfeedfoward
            VV = edmd.exec_dmd()

            calc_dat = []
            for i in range(edmd.nvar):
                calc_dat.append(VV[i])
            calc_dat = np.array(calc_dat)
            obs_calc_mat = obs.obs_function_for_adj(VV)

            p_t_prev = np.zeros_like(edmd.p1_init)

            if t_j == 0:
                nv, nt, ny, nx = np.shape(VV)
                obs.set_obsmatrix_for_adj(nv, ny, nx, nt, OBS_STEP)

            for t_k in range(int(CONFIG["DA"]["assimilate_steps"])):
                k_step = int(CONFIG["DA"]["assimilate_steps"]) - (t_k + 1)

                # アドジョイントコードの実行
                p_t_temp = edmd.adjoint_code(p_t_prev)

                # 誤差評価
                error_temp = (
                    obs_calc_mat[:, :, :, k_step] - obs.obs_mat[:, :, :, k_step]
                ) * obs.obs_mask[:, :, :, k_step]

                # 観測誤差行列の穴埋め用のマトリクス
                error_mat = np.zeros((np.size(p_t_prev), 1))
                state_size = edmd.nvar * edmd.ny * edmd.nx
                error_mat[:state_size, :] = np.reshape(error_temp.real, (state_size, 1))

                p_t = p_t_temp + error_mat
                p_t_prev = p_t

            # 更新
            edmd.p1_init -= float(CONFIG["DA"]["alpha"]) * np.real(p_t.flatten())


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


CONFIG = read_config()

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
    exec_simulation()

