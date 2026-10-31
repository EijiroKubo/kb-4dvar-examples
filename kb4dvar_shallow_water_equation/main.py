import os
import yaml
import numpy as np
import pandas as pd

import calc_sea_field as SEA
import calc_edmd_adj as EDMD_ADJ

import obsdata_controller

import calc_dct as DCT2

import matplotlib.pyplot as plt
from matplotlib import animation


def read_sea_field():
    """osに応じたsharedobjectを指定して読み込む"""
    if os.name == "nt":
        f = os.path.join("seafield", "CalcSeaFieldMod.dll")
    else:
        f = os.path.join("seafield", "CalcSeaFieldMod.so")

    module = SEA.SEAFIELD(f)

    return module


def read_config():
    """計算条件はCONFIG.ymlファイルから読みとる"""

    with open("config.yml", "r", encoding="utf-8") as yml:
        config = yaml.safe_load(yml)
    return config




class ForwardCalc:
    def __init__(self):

        # 設定ファイルの読み込み
        self.config = read_config()

        # クラス継承
        self.sea = read_sea_field()

        self.sea.common_settings(
            self.config["MODEL"]["dt"],
            self.config["MODEL"]["time_proceed"],
            self.config["MODEL"]["t_update"],
            self.config["MODEL"]["number_of_threads"],
        )

        # 一応初期化
        self.sea.initialize_condition(
            self.config["MODEL"]["initial_depth"],
            self.config["MODEL"]["read_initial_file"],
            self.config["MODEL"]["i_num"],
            self.config["MODEL"]["j_num"],
            self.config["MODEL"]["dy"],
            self.config["MODEL"]["dx"],
        )

    def setting_init_state(self, p1_external: np.ndarray = None):
        """初期浸水深を直接いじる"""

        if p1_external is None:
            # ここが同化前の初期計算条件となる。
            self.sea.h_time_sr[
                0,
                X_IMPACT - X_IMP_RANGE : X_IMPACT + X_IMP_RANGE,
                Y_IMPACT - Y_IMP_RANGE : Y_IMPACT + Y_IMP_RANGE,
            ] += IMPACT_DEPTH

        else:
            # 初期条件
            # self.sea.h_time_sr[0,:,:] = p1_external
            self.sea.h_time_sr[0, :, :] = p1_external[0]
            self.sea.u_time_sr[0, :, :] = p1_external[1]
            self.sea.v_time_sr[0, :, :] = p1_external[2]

    def calculation(self):

        # 物理計算を実行
        self.sea.proceed_sea_field_with_velocity()

        # 注意：計算結果のh_time_sr は [時間/横/縦] のインデックスを有している
        self.sea.allocate_reshaped_res()



def exec_simulation(save_anim: bool = False):

    # 物理計算による時間発展計算クラス
    forward = ForwardCalc()

    # ex-dmdの設定
    edmd = EDMD_ADJ.EDMD(
        forward.sea.t_no, forward.sea.j_no, forward.sea.i_no, N_VARIABLE
    )

    # 離散コサイン変換
    dct = DCT2.DCT(forward.sea.j_no, forward.sea.i_no)
    elim_ny = int((1 - CUT_RATE) * forward.sea.i_no)
    elim_nx = int((1 - CUT_RATE) * forward.sea.j_no)

    # 観測データ
    obs = obsdata_controller.ObsDataProcessor()
    obs.read_observation()

    # 観測誤差分散行列 -> 一旦単位行列で設定
    Rk = np.eye(obs.obs_num)

    # 更新対象とする初期条件
    p1_temp = None

    # 以下、同化計算
    for t_i in range(int(CONFIG["DA"]["global_loops"])):
        # 物理計算
        forward.setting_init_state(p1_temp)
        forward.calculation()

        # INPOOL2の対応
        if DIMENSION_STANDARDIZE == True:
            calc_dat = []

            # Falseの場合と同じ3変数(h, u, v)のみを標準化してEDMDに入力する
            h_mean = forward.sea.h_time_sr.mean()
            u_mean = forward.sea.u_time_sr.mean()
            v_mean = forward.sea.v_time_sr.mean()

            h_std = max(forward.sea.h_time_sr.std(), 0.001)
            u_std = max(forward.sea.u_time_sr.std(), 0.001)
            v_std = max(forward.sea.v_time_sr.std(), 0.001)

            calc_dat.append((forward.sea.h_time_sr - h_mean) / h_std)
            calc_dat.append((forward.sea.u_time_sr - u_mean) / u_std)
            calc_dat.append((forward.sea.v_time_sr - v_mean) / v_std)

            calc_dat = np.array(calc_dat)
            edmd.exec_mode_decomposition_from_snapshops(calc_dat)

            '''
            hu_temp = forward.sea.h_time_sr * forward.sea.u_time_sr
            hv_temp = forward.sea.h_time_sr * forward.sea.v_time_sr
            uv_temp = forward.sea.u_time_sr * forward.sea.v_time_sr
            vv_temp = forward.sea.v_time_sr * forward.sea.v_time_sr
            uu_temp = forward.sea.u_time_sr * forward.sea.u_time_sr

            # 標準化処理（ x - μ ) / ∑
            h_mean = forward.sea.h_time_sr.mean()
            u_mean = forward.sea.u_time_sr.mean()
            v_mean = forward.sea.v_time_sr.mean()
            hu_mean = hu_temp.mean()
            hv_mean = hv_temp.mean()
            uv_mean = uv_temp.mean()
            vv_mean = vv_temp.mean()
            uu_mean = uu_temp.mean()

            h_std = max(forward.sea.h_time_sr.std(), 0.001)
            u_std = max(forward.sea.u_time_sr.std(), 0.001)
            v_std = max(forward.sea.v_time_sr.std(), 0.001)
            hu_std = max(hu_temp.std(), 0.001)
            hv_std = max(hv_temp.std(), 0.001)
            uv_std = max(uv_temp.std(), 0.001)
            vv_std = max(vv_temp.std(), 0.001)
            uu_std = max(uu_temp.std(), 0.001)

            h_norm = (forward.sea.h_time_sr - h_mean) / h_std
            u_norm = (forward.sea.u_time_sr - u_mean) / u_std
            v_norm = (forward.sea.v_time_sr - v_mean) / v_std
            hu_norm = (hu_temp - hu_mean) / hu_std
            hv_norm = (hv_temp - hv_mean) / hv_std
            uv_norm = (uv_temp - uv_mean) / uv_std
            vv_norm = (vv_temp - vv_mean) / vv_std
            uu_norm = (uu_temp - uu_mean) / uu_std

            calc_dat.append(h_norm)
            calc_dat.append(u_norm)
            calc_dat.append(v_norm)
            calc_dat.append(uu_norm)
            calc_dat.append(vv_norm)
            calc_dat.append(hu_norm)
            calc_dat.append(hv_norm)
            '''

            calc_dat = np.array(calc_dat)

            edmd.exec_mode_decomposition_from_snapshops(calc_dat)

        else:
            calc_dat = []
            calc_dat.append(forward.sea.h_time_sr)
            calc_dat.append(forward.sea.u_time_sr)
            calc_dat.append(forward.sea.v_time_sr)
            calc_dat = np.array(calc_dat)

            edmd.exec_mode_decomposition_from_snapshops(calc_dat)

        # 随伴コード
        edmd.prepare_adj()

        #
        if save_anim == True:
            if DIMENSION_STANDARDIZE == True:
                params_st = np.array(
                    [
                        [h_mean, h_std],
                        [u_mean, u_std],
                        [v_mean, v_std],
                    ],
                    dtype=float,
                )
            else:
                params_st = None

            edmd.show_result(calc_dat, int((t_i + 1) * 100), param_standards=params_st)

        # 擬似随伴コードによる計算
        for t_j in range(int(CONFIG["DA"]["local_loops"])):
            VV = edmd.exec_dmd()

            calc_dat = []
            for i in range(edmd.nvar):
                calc_dat.append(VV[i])
            calc_dat = np.array(calc_dat)

            # check
            if save_anim == True:
                header_id_tmp = int((t_i + 1) * 100) + int(t_j + 1)
                edmd.show_result(VV, header_id_tmp, param_standards=params_st)

            obs_calc_mat = obs.obs_function_for_adj(VV[0])  # 水位だけを対象にする
            p_t_prev = np.zeros_like(edmd.p1_init)

            if t_j == 0:
                nv, nt, ny, nx = np.shape(VV)

                if DIMENSION_STANDARDIZE == True:
                    # 観測データが無いところも標準化すると、誤差が拡大するため
                    # 観測データのあるメッシュをmaskとして作成しておく
                    obs_mask = np.zeros((ny, nx), dtype=bool)

                    for ip in range(np.shape(obs.obs_point)[0]):
                        jj = int(obs.obs_point[ip][0].item())  # X
                        ii = int(obs.obs_point[ip][1].item())  # Y
                        obs_mask[ii, jj] = 1  # True

                obs.set_obsmatrix_for_adj(ny, nx, nt)

            if DIMENSION_STANDARDIZE == True:
                obs_mat_norm = (obs.obs_mat - h_mean) / h_std
                obs_mat_norm *= obs_mask[:, :, None]

            for t_k in range(int(CONFIG["DA"]["assimilate_steps"])):
                k_step = int(CONFIG["DA"]["assimilate_steps"]) - (t_k + 1)

                # アドジョイントコードの実行
                p_t_temp = edmd.adjoint_code(p_t_prev)

                # 誤差評価
                if DIMENSION_STANDARDIZE == True:
                    error_temp = (
                        obs_calc_mat[:, :, k_step] - obs_mat_norm[:, :, k_step]
                    ) * obs_mask
                else:
                    error_temp = obs_calc_mat[:, :, k_step] - obs.obs_mat[:, :, k_step]

                # 観測誤差行列の穴埋め用のマトリクス
                error_mat = np.zeros((np.size(p_t_prev), 1))
                error_mat[: np.size(error_temp), :] = np.reshape(
                    error_temp.real, (np.size(error_temp), 1)
                )

                p_t = p_t_temp + error_mat
                p_t_prev = p_t

            # 更新
            edmd.p1_init -= float(CONFIG["DA"]["alpha"]) * np.real(p_t.flatten())
            pt_phys = np.real(p_t.flatten()[: edmd.xsize]).reshape(
                edmd.nvar, forward.sea.j_no, forward.sea.i_no
            )
            for ivar, name in enumerate(["h", "u", "v"]):
                print(
                    f"p_t {name}: "
                    f"rms={np.sqrt(np.mean(pt_phys[ivar] ** 2)):.6e}, "
                    f"maxabs={np.max(np.abs(pt_phys[ivar])):.6e}"
                )

        p1_temp = np.reshape(
            edmd.p1_init[: forward.sea.i_no * forward.sea.j_no * edmd.nvar],
            (edmd.nvar, forward.sea.j_no, forward.sea.i_no),
        )

        # 以下、chatGPTより
        nxy = forward.sea.i_no * forward.sea.j_no

        # h成分だけに仮想誤差を入れる
        pk = np.zeros_like(edmd.p1_init)
        pk[:nxy] = 1.0

        adj = edmd.adjoint_code(pk).real.flatten()

        adj_phys = adj[: edmd.xsize].reshape(
            edmd.nvar,
            forward.sea.j_no,
            forward.sea.i_no,
        )
        for ivar, name in enumerate(["h", "u", "v"]):
            rms = np.sqrt(np.mean(adj_phys[ivar] ** 2))
            maxabs = np.max(np.abs(adj_phys[ivar]))


        if DIMENSION_STANDARDIZE == True:
            # 標準化量→物理量　に戻してから物理モデルの初期状態を更新する
            p1_temp[0] = p1_temp[0] * h_std + h_mean
            p1_temp[1] = p1_temp[1] * u_std + u_mean
            p1_temp[2] = p1_temp[2] * v_std + v_mean

        # dctで高周波成分を除去
        if CUT_RATE < 1:
            # for i in range(edmd.nvar):
            p1_dct = dct.dct_2d(p1_temp[0])
            p1_dct[elim_ny:, :] = 0
            p1_dct[:, elim_nx:] = 0
            p1_idct = dct.inv_dct_2d(p1_dct)
            p1_temp[0] = p1_idct


    # 補正後予測計算
    forward.setting_init_state(p1_temp)
    forward.calculation()

    # アニメーション表示用データ
    calc_dat = []
    calc_dat.append(forward.sea.h_time_sr)
    calc_dat.append(forward.sea.u_time_sr)
    calc_dat.append(forward.sea.v_time_sr)
    calc_dat = np.array(calc_dat)

    edmd.show_result(calc_dat, 10001)

    # edmdでの計算 -> 最後の計算結果
    edmd.show_result(VV, 10000)


Y_IMPACT = 100
X_IMPACT = 100
Y_IMP_RANGE = 8
X_IMP_RANGE = 8
IMPACT_DEPTH = 1.0
CONFIG = read_config()
SAVE_ANIM = False
CUT_RATE = 1 - float(CONFIG["DMD"]["dct_threshold"])  # 0.2

N_VARIABLE = 3
DIMENSION_STANDARDIZE = False 


if __name__ == "__main__":

    # カレントディレクトリの位置を変更
    os.chdir(os.path.dirname(__file__))

    # 計算
    exec_simulation(SAVE_ANIM)
