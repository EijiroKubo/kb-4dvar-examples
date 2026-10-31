
import numpy as np
import matplotlib.pyplot as plt


class DCT:
    """ 離散コサイン変換のためのプログラム    
    Note: 以下の参考ページのものを二次元ように修正していっている
        参考: https://tony-mooori.blogspot.com/2016/02/dctpythonpython.html
        [1]『画像処理とパターン認識入門』酒井幸市 著        
        [2] scipy.fftpack.dct http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.fftpack.dct.html
    """       

    def __init__(self, H, W):
        self.H = H  # 縦のデータ数
        self.W = W  # 横のデータ数
        
        # 縦方向(H)用の変換行列
        self.phi_H = np.array([self.generate_phi(i, self.H) for i in range(self.H)])
        # 横方向(W)用の変換行列
        self.phi_W = np.array([self.generate_phi(i, self.W) for i in range(self.W)])

    def generate_phi(self, k, N):
        """ 指定されたサイズNに対するDCT基底関数を生成 """
        if k == 0:
            return np.ones(N) / np.sqrt(N)
        else:
            return np.sqrt(2.0 / N) * np.cos((k * np.pi / (2 * N)) * (np.arange(N) * 2 + 1))

    def dct_2d(self, data):
        """ 
        2次元離散コサイン変換 (非正方形対応)
        公式: C = A_H * data * (A_W)^T
        """
        # dataが(H, W)の場合
        # self.phi_H.dot(data) -> (H, H) * (H, W) = (H, W)
        # .dot(self.phi_W.T)   -> (H, W) * (W, W) = (H, W)
        return self.phi_H.dot(data).dot(self.phi_W.T)

    def inv_dct_2d(self, c):
        """ 
        2次元離散コサイン逆変換 (非正方形対応)
        公式: data = (A_H)^T * c * A_W
        """
        return self.phi_H.T.dot(c).dot(self.phi_W)
    
