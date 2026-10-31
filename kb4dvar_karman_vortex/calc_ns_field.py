import os
import numpy as np
from ctypes import c_void_p


class NSFWD:
    def __init__(self, libfile: str):

        # Optional directory containing MinGW/gfortran runtime DLLs on Windows.
        # Keep the handle alive while the loaded library is in use.
        self._dll_directory_handle = None
        if os.name == "nt":
            dll_dir = os.environ.get("KB4DVAR_DLL_DIR")
            if dll_dir:
                self._dll_directory_handle = os.add_dll_directory(dll_dir)

        so_abspath = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + libfile
        self.libs = np.ctypeslib.load_library(so_abspath, ".")

        # 計算条件を条件ファイルから読み取る
        self.read_initial_condition = self.libs.read_initial_condition
        self.read_initial_condition.restype = c_void_p

        # メッシュ条件をpython側にも転記する
        self.return_mesh_num = self.libs.return_mesh_num
        self.return_mesh_num.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.int32),  # imax
            np.ctypeslib.ndpointer(dtype=np.int32),  # imx1
            np.ctypeslib.ndpointer(dtype=np.int32),  # jmax
            np.ctypeslib.ndpointer(dtype=np.int32),  # jmx1
            np.ctypeslib.ndpointer(dtype=np.int32),  # t_steps_tmp
        ]
        self.return_mesh_num.restype = c_void_p

        self.imax = None
        self.imx1 = None
        self.jmax = None
        self.jmx1 = None
        self.t_no = None
        self.qin_time_sr = None

        # 前進計算を行う
        self.calc_karman_vortex = self.libs.calc_karman_vortex
        self.calc_karman_vortex.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),
        ]
        # qin_time_sr
        self.calc_karman_vortex.restype = c_void_p

        # 元々の4dVarケースと同じように初期化計算を実行する
        self.set_validation_case = self.libs.set_validation_case
        self.set_validation_case.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),
        ]
        self.set_validation_case.restype = c_void_p

    def initialize_condition(self, t_no: int):

        imax = np.array([0], dtype=np.int32)
        imx1 = np.array([0], dtype=np.int32)
        jmax = np.array([0], dtype=np.int32)
        jmx1 = np.array([0], dtype=np.int32)
        t_steps = np.array([0], dtype=np.int32)

        self.return_mesh_num(imax, imx1, jmax, jmx1, t_steps)

        self.imax = int(imax[0])
        self.imx1 = int(imx1[0])
        self.jmax = int(jmax[0])
        self.jmx1 = int(jmx1[0])
        self.t_no = int(t_no)

        # Fortran と合わせて (imx1, jmx1, 3, t_no)
        self.qin_time_sr = np.zeros(
            (self.imx1, self.jmx1, 3, self.t_no), dtype=np.float64, order="F"
        )
        return self.qin_time_sr

    def set_init_state(self, u0: np.ndarray, v0: np.ndarray, p0: np.ndarray):
        # self._require_state()
        self.qin_time_sr[:, :, :, :] = 0.0

        self.qin_time_sr[0 : self.imx1, 0 : self.jmax, 0, 0] = u0
        self.qin_time_sr[0 : self.imax, 0 : self.jmx1, 1, 0] = v0
        self.qin_time_sr[0 : self.imax, 0 : self.jmax, 2, 0] = p0

        return self.qin_time_sr

    def proceed_forward(self):
        # self._require_state()
        self.calc_karman_vortex(self.qin_time_sr)
        return self.qin_time_sr

    def set_validation_initial_state(self):
        # self._require_state()
        self.set_validation_case(self.qin_time_sr)
        return self.qin_time_sr

    def get_result(self):
        # self._require_state()
        # Python側で扱いやすいように (t, y, x) に寄せる
        u_time_sr = np.transpose(
            self.qin_time_sr[0 : self.imx1, 0 : self.jmax, 0, :], (2, 1, 0)
        )
        v_time_sr = np.transpose(
            self.qin_time_sr[0 : self.imax, 0 : self.jmx1, 1, :], (2, 1, 0)
        )
        p_time_sr = np.transpose(
            self.qin_time_sr[0 : self.imax, 0 : self.jmax, 2, :], (2, 1, 0)
        )
        return u_time_sr, v_time_sr, p_time_sr

    def _require_state(self):
        if self.qin_time_sr is None:
            raise RuntimeError("initialize_condition(t_no) must be called first.")

