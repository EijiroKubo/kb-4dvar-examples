import os
import yaml
import numpy as np
from numpy import dot, multiply, diag, power
from numpy.linalg import inv, eig, pinv
from scipy.linalg import svd

import matplotlib.pyplot as plt
from matplotlib import animation


class EDMD:
    def __init__(self, n_snaps: int, nx: int, ny: int, n_variables: int = 1):

        # 基本設定
        with open("config.yml", "r", encoding="utf-8") as yml:
            config = yaml.safe_load(yml)

        self.rank = config["DMD"]["rank"]  # dmdモードの行列ランク
        self.dt = config["DMD"]["dt"]  # dmdの⊿t
        self.nlift = config["DMD"]["NumberOfLifts"]  # Liftの数
        self.substract_ave = config["DMD"]["substract_ave"]

        self.tt = 0.0
        self.N = n_snaps - 1  # dmdモード抽出用の元計算結果のスナップ数
        self.nx = nx
        self.ny = ny
        self.nvar = n_variables
        self.xsize = ny * nx * n_variables
        self.t_save = self.N  # 出力するスナップショット数


    def createmat(self, X, N, h=1):
        # モード分解を行うためのinputデータを作成

        print("create snapshots")
        p1 = np.array([])
        p2 = np.array([])

        x_size = np.size(X[0])

        for i in range(N - h + 1):
            for j in range(h):
                p1 = np.append(p1, X[i + j])
                p2 = np.append(p2, X[i + j + 1])

        p1 = np.reshape(p1, (N - h + 1, x_size * h))
        p1 = np.transpose(p1)
        p2 = np.reshape(p2, (N - h + 1, x_size * h))
        p2 = np.transpose(p2)

        if self.substract_ave == True:
            self.x_ave = np.mean(X)
            p1 -= self.x_ave
            p2 -= self.x_ave

        return p1, p2

    def createmat_nval(self, Xs: np.ndarray, N: int):
        # モード分解を行うためのinputデータを作成

        print("create snapshots with integrate each variables")
        p1 = np.zeros((self.xsize, N))
        p2 = np.zeros((self.xsize, N))

        for i in range(N):
            for j in range(self.nvar):
                p1_temp = Xs[j][i].flatten()
                p2_temp = Xs[j][i + 1].flatten()

                p1[j * self.ny * self.nx : (j + 1) * self.ny * self.nx, i] = p1_temp
                p2[j * self.ny * self.nx : (j + 1) * self.ny * self.nx, i] = p2_temp

        nxy = self.nx * self.ny

        if self.substract_ave == True:
            self.x_ave = np.zeros((self.xsize, 1))
            for j in range(self.nvar):
                st = j * nxy
                ed = (j + 1) * nxy
                m_j = p1[st:ed, :].mean()
                self.x_ave[st:ed, 0] = m_j
            p1 -= self.x_ave
            p2 -= self.x_ave

        return p1, p2

    def extend_mat_with_lifting(self, p1, p2):
        """拡張モード分解のためにリフト分を考慮して対となる状態空間行列を作成"""

        p1_ext = []
        p2_ext = []

        for k in range(self.N - 1):
            p1_ext_temp = []
            p2_ext_temp = []

            for m in range(self.nlift):
                ext1_temp = p1[:, k] ** (m + 1)
                ext2_temp = p2[:, k] ** (m + 1)

                p1_ext_temp.append(ext1_temp)
                p2_ext_temp.append(ext2_temp)

            p1_ext_temp = np.array(p1_ext_temp)
            p2_ext_temp = np.array(p2_ext_temp)
            p1_ext.append(p1_ext_temp.flatten())
            p2_ext.append(p2_ext_temp.flatten())

        p1_ext = np.reshape(
            np.transpose(p1_ext),
            (self.nx * self.ny * self.nlift * self.nvar, self.N - 1),
        )
        p2_ext = np.reshape(
            np.transpose(p2_ext),
            (self.nx * self.ny * self.nlift * self.nvar, self.N - 1),
        )

        return p1_ext, p2_ext

    def lifting_one_state(self, p1):
        """拡張モード分解のためのリフト考慮
        Args:
            - p1: 初期の浸水深分布
        Note: 初期条件を変更する場合に使用する
        """

        p1_ext = []

        for m in range(self.nlift):
            ext1_temp = p1 ** (m + 1)
            p1_ext.append(ext1_temp)

        p1_ext = np.array(p1_ext)
        return p1_ext.flatten()

    def mode_decomposing(self, p1, p2):
        # モード分解
        print("mode decomposition")

        self.U2, self.Sig2, self.Vh2 = svd(p1, False)
        U = self.U2[:, : self.rank]
        Sig = np.diag(self.Sig2)[: self.rank, : self.rank]
        V = self.Vh2.conj().T[:, : self.rank]

        F = np.dot(np.dot(np.dot(U.conj().T, p2), V), np.linalg.pinv(Sig))
        mu, W = eig(F)

        phi = np.dot(U, W)
        # phi = dot(dot(dot( p2, V), inv(Sig)), W) #記事2に詳述

        mu = np.nan_to_num(mu, nan=0)
        phi = np.nan_to_num(phi, nan=0)

        return mu, phi

    def exec_mode_decomposition_from_snapshops(self, X):
        """計算結果の時系列変化（スナップショット）からdmdを実行"""

        if self.nvar > 1:
            # スナップショットに複数の変数が含まれる場合
            # input:Xは[ val, time , x , y ]となる
            _p1_tmp, _p2_tmp = self.createmat_nval(X, self.N)

        else:
            _p1_tmp, _p2_tmp = self.createmat(X, self.N)

        _p1, _p2 = self.extend_mat_with_lifting(_p1_tmp, _p2_tmp)
        self.mu, self.phi = self.mode_decomposing(_p1, _p2)
        self.p1_init = _p1[:, 0]

    def apply_dmdmode(self, d_length: float = None):
        # dmdモデルを適用して浸水深時系列の元データを構築する

        if d_length is None:
            tt = np.linspace(0, self.N * self.dt, self.t_save)
        else:
            tt = np.linspace(0, d_length * self.dt, self.t_save)

        UU = dot(pinv(self.phi), self.p1_init)
        Psi = np.zeros([self.rank, len(tt)], dtype="complex")

        for i, _t in enumerate(tt):
            Psi[:, i] = multiply(power(self.mu, _t / self.dt), UU)

        return np.dot(self.phi, Psi)

    def exec_dmd(self, p1_external=None):

        if p1_external is None:
            pass
        else:
            self.p1 = self.lifting_one_state(p1_external)

        VV = self.apply_dmdmode()  # モードからデータを再現
        VV = VV[: self.nx * self.ny * self.nvar, :]  # 下半分は無視

        if self.nvar > 1:
            v_reshaped = np.reshape(VV, (self.nvar, self.ny, self.nx, self.t_save))
            v_reshaped = np.transpose(v_reshaped, (0, 3, 1, 2))
            return v_reshaped

        else:
            v_reshaped = np.reshape(VV, (self.ny, self.nx, self.t_save))
            v_reshaped = np.transpose(v_reshaped, (2, 0, 1))

            if self.substract_ave == True:
                return v_reshaped + self.x_ave
            else:
                return v_reshaped

    def prepare_adj(self, d_length: float = None):
        """adj実施前に時間軸を設定"""

        if d_length is None:
            self.tt = np.linspace(0, self.N * self.dt, self.t_save)
        else:
            self.tt = np.linspace(0, d_length * self.dt, self.t_save)

        # 複素数行列ならconjをconjで随伴行列にしないと、との指摘
        self.phi_t = self.phi.conj().T
        self.phi_t_inv = pinv(self.phi_t)
        self.mu_conj = np.reshape(self.mu.conj(), (self.rank, 1))

    def adjoint_code(self, pk):
        """各モードを使用して、アドジョイント計算を行う
        """

        uu_k_t = dot(self.phi_t, pk).reshape((self.rank, 1))  # -> size=[rank]
        psi_k_t = multiply(self.mu_conj, uu_k_t)  # -> size=[rank]
        adj_res = dot(self.phi_t_inv, psi_k_t)  # -> size=[st-space]

        return adj_res

    def show_result(
        self,
        dat,
        header_id: int = None,
        layout: str = "row",
        param_standards: np.ndarray = None,
    ):
        """
        dat: shape = (nvar, nt, ny, nx)
            nvar: 変数数 (self.nvar)
            nt  : 時間ステップ数
            ny  : 格子数 (self.ny)
            nx  : 格子数 (self.nx)
        """
        # ---- 形状チェック（最低限） ----
        if dat.ndim != 4:
            raise ValueError(
                f"dat.ndim = {dat.ndim} ですが、(nvar, t, ny, nx) の4次元を期待しています。"
            )

        nvar, nt, ny, nx = dat.shape

        if (nvar != self.nvar) or (ny != self.ny) or (nx != self.nx):
            raise ValueError(
                f"dat.shape = {dat.shape} ですが、(nvar={self.nvar}, t, ny={self.ny}, nx={self.nx}) を期待しています。"
            )

        # ---- subplot レイアウト ----
        if layout == "row":
            nrows, ncols = 1, nvar
        elif layout == "grid":
            ncols = int(np.ceil(np.sqrt(nvar)))
            nrows = int(np.ceil(nvar / ncols))
        else:
            raise ValueError("layout は 'row' か 'grid' を指定してください。")

        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(4.8 * ncols, 4 * nrows),
            squeeze=False,
        )

        # ---- 変数ごとのカラースケール・カラーマップ ----
        vmins = []
        vmaxs = []
        cmaps = []

        for v in range(nvar):
            if v == 0:
                vmins.append(9.5)
                vmaxs.append(10.5)
                cmaps.append("viridis")
            else:
                vmins.append(-0.05)
                vmaxs.append(0.05)
                cmaps.append("seismic")

        # ---- カラーバーを1回だけ作るための初期 imshow ----
        base_images = []

        for v in range(nvar):
            r, c = divmod(v, ncols)
            ax = axes[r, c]

            ax.set_title(f"var {v}")
            ax.set_xticks([])
            ax.set_yticks([])

            img0 = np.real(dat[v, 0, :, :])

            if param_standards is not None:
                # 標準化を戻す
                img0 = img0 * param_standards[v][1] + param_standards[v][0]

            im0 = ax.imshow(
                img0,
                animated=True,
                vmin=vmins[v],
                vmax=vmaxs[v],
                cmap=cmaps[v],
            )

            base_images.append(im0)

            # 各 subplot にカラーバーを付ける
            fig.colorbar(
                im0,
                ax=ax,
                fraction=0.046,
                pad=0.04,
            )

        # 未使用 axes を非表示にする
        for v in range(nvar, nrows * ncols):
            r, c = divmod(v, ncols)
            axes[r, c].axis("off")

        # ---- アニメフレーム作成 ----
        ims = []

        for it in range(nt):
            artists = []

            for v in range(nvar):
                r, c = divmod(v, ncols)
                ax = axes[r, c]

                img = np.real(dat[v, it, :, :])

                if param_standards is not None:
                    # 標準化を戻す
                    img = img * param_standards[v][1] + param_standards[v][0]

                im = ax.imshow(
                    img,
                    animated=True,
                    vmin=vmins[v],
                    vmax=vmaxs[v],
                    cmap=cmaps[v],
                )

                artists.append(im)

            ims.append(artists)

        anim = animation.ArtistAnimation(
            fig,
            ims,
            interval=10,
            repeat_delay=100,
            repeat=False,
        )

        out_dir = "result_anim"
        os.makedirs(out_dir, exist_ok=True)

        if header_id is None:
            fname = "03_edmd_anim.gif"
        else:
            fname = f"{header_id}_edmd_anim.gif"

        fpath = os.path.join(out_dir, fname)
        anim.save(fpath, writer="pillow", fps=10)

        plt.close(fig)
