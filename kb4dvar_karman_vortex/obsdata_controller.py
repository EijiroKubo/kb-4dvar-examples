import os
import zipfile
import numpy as np
import pandas as pd


class observation:
    def __init__(self):

        # データ同化流体科学 表5.14を参照に4メッシュ間隔で設定
        x_tmp = np.arange(0, 240, 4)
        y_tmp = np.arange(0, 80, 4)

        data_tmp = []
        for j in x_tmp:
            for i in y_tmp:
                data_tmp.append([j, i])

        self.data = np.array(data_tmp)


class ObsDataProcessor:
    def __init__(self):
        self.obs_pass = os.path.join(os.getcwd(), "observe")
        self.obs_point = np.empty((0, 3), dtype=int)
        self.obs = np.empty((0, 0), dtype=float)
        self.obs_num = 0
        self.nvar = 1

    def set_obs_place(self, obs_point):
        self.obs_point = np.asarray(obs_point, dtype=int)

    @staticmethod
    def _observation_index(filename: str) -> int:
        stem = os.path.splitext(filename)[0]
        return int(stem.split("_")[-1])

    @staticmethod
    def _point_parts(point):
        point_arr = np.asarray(point).reshape(-1)
        if point_arr.size == 2:
            return 0, int(point_arr[0]), int(point_arr[1])
        return int(point_arr[0]), int(point_arr[1]), int(point_arr[2])

    def obs_function(self, X):
        """X[y,x,t] or X[var,t,y,x]"""
        x_temp = []
        for i in range(np.shape(self.obs_point)[0]):
            var, jj, ii = self._point_parts(self.obs_point[i])
            if X.ndim == 4:
                x_temp.append(X[var, :, ii, jj])
            else:
                x_temp.append(X[ii, jj, :])

        x_temp = np.array(x_temp)
        return x_temp

    def obs_function_for_adj(self, X):
        """観測行列のマトリクスを生成する
        Input: X ( t, ny, nx ) or ( nvar, t, ny, nx )
        Output: ( ny, nx, t ) or ( nvar, ny, nx, t )
        """

        if X.ndim == 4:
            nvar, t, ny, nx = np.shape(X)
            x_temp = np.zeros((nvar, ny, nx, t), dtype=X.dtype)

            for i in range(np.shape(self.obs_point)[0]):
                var, jj, ii = self._point_parts(self.obs_point[i])
                x_temp[var, ii, jj, :] = X[var, :, ii, jj].real

            return x_temp

        t, ny, nx = np.shape(X)
        x_temp = np.zeros((ny, nx, t), dtype=X.dtype)

        for i in range(np.shape(self.obs_point)[0]):
            _, jj, ii = self._point_parts(self.obs_point[i])
            x_temp[ii, jj, :] = X[:, ii, jj].real

        return x_temp

    def set_obsmatrix_for_adj(self, *shape_args):
        """実観測データを観測行列の形に整形しておく"""

        obs_step = 1
        if len(shape_args) == 5:
            nvar, ny, nx, nt, obs_step = shape_args
        elif len(shape_args) == 4:
            nvar, ny, nx, nt = shape_args
        else:
            ny, nx, nt = shape_args

        obs_step = int(obs_step)

        if len(shape_args) >= 4:
            self.obs_mat = np.zeros((nvar, ny, nx, nt))
            self.obs_mask = np.zeros((nvar, ny, nx, nt))

            for i in range(np.shape(self.obs_point)[0]):
                var, jj, ii = self._point_parts(self.obs_point[i])
                t_ids = np.arange(np.size(self.obs[i])) * obs_step
                t_ids = t_ids[t_ids < nt]
                iter_num = np.size(t_ids)
                self.obs_mat[var, ii, jj, t_ids] = self.obs[i][:iter_num]
                self.obs_mask[var, ii, jj, t_ids] = 1.0

            return

        self.obs_mat = np.zeros((ny, nx, nt))
        self.obs_mask = np.zeros((ny, nx, nt))

        for i in range(np.shape(self.obs_point)[0]):
            _, jj, ii = self._point_parts(self.obs_point[i])
            t_ids = np.arange(np.size(self.obs[i])) * obs_step
            t_ids = t_ids[t_ids < nt]
            iter_num = np.size(t_ids)
            self.obs_mat[ii, jj, t_ids] = self.obs[i][:iter_num]
            self.obs_mask[ii, jj, t_ids] = 1.0

    def create_test_observation(self, X):
        obs = self.obs_function(X)

        t = 0
        for i in range(np.shape(self.obs_point)[0]):
            f_name = "observations_" + str(i) + ".csv"
            f = os.path.join(self.obs_pass, f_name)
            df = pd.DataFrame(obs[i])
            df.to_csv(f, sep=",", index=None, header=None)

            f_place_name = "obs_place_" + str(i) + ".csv"
            f_place = os.path.join(self.obs_pass, f_place_name)
            _, jj, ii = self._point_parts(self.obs_point[i])
            d = np.array([jj, ii])
            df_place = pd.DataFrame(d)
            df_place.to_csv(f_place, sep=",", index=None, header=None)

    def read_observation(self, obs_sources=None):
        """観測データを取得"""

        if obs_sources is None:
            obs_sources = [(0, "observe")]

        obs = []
        obs_point_temp = []
        default_obs_points = observation().data

        for var_idx, obs_dir_name in obs_sources:
            obs_dir = obs_dir_name
            if not os.path.isabs(obs_dir):
                obs_dir = os.path.join(os.getcwd(), obs_dir)

            source_observations = []
            if zipfile.is_zipfile(obs_dir):
                with zipfile.ZipFile(obs_dir) as archive:
                    obs_files = [
                        f
                        for f in archive.namelist()
                        if os.path.basename(f).startswith("observations_")
                        and f.endswith(".csv")
                    ]
                    obs_files = sorted(obs_files, key=self._observation_index)
                    for obs_file in obs_files:
                        with archive.open(obs_file) as f_obs:
                            df = pd.read_csv(f_obs, header=None)
                        source_observations.append(
                            (obs_file, np.array(df).reshape(-1))
                        )
            else:
                obs_files = [
                    f
                    for f in os.listdir(obs_dir)
                    if f.startswith("observations_") and f.endswith(".csv")
                ]
                obs_files = sorted(obs_files, key=self._observation_index)
                for obs_file in obs_files:
                    f_obs = os.path.join(obs_dir, obs_file)
                    df = pd.read_csv(f_obs, header=None)
                    source_observations.append((obs_file, np.array(df).reshape(-1)))

            for obs_file, d_temp in source_observations:
                obs_idx = self._observation_index(obs_file)
                obs.append(d_temp)

                d_point_temp = default_obs_points[obs_idx]
                obs_point_temp.append(
                    [int(var_idx), int(d_point_temp[0]), int(d_point_temp[1])]
                )

        self.obs_num = len(obs)
        self.set_obs_place(obs_point_temp)
        self.nvar = int(np.max(self.obs_point[:, 0])) + 1

        step_num = min(np.size(obs_i) for obs_i in obs)
        self.obs = np.array([obs_i[:step_num] for obs_i in obs])
