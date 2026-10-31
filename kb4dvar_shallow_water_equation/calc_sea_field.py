import os 
import numpy as np
import pandas as pd 
from ctypes import *

class SEAFIELD:
    def __init__(self, libfile):

        if os.name == "nt": 
            # windows環境で計算実行する場合、依存DLLのパスを追加            
            os.add_dll_directory(r"C:\TDM-GCC-64\bin")
        else: 
            pass 

        so_abspath = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + libfile       
        self.libs = np.ctypeslib.load_library( so_abspath , "." )

        #計算ループ数、出力間隔（ステップ）数、Δtの共有用
        self.settings = self.libs.settings
        self.settings.argtypes = [ 
            np.ctypeslib.ndpointer(dtype=np.float64),     # dt
            np.ctypeslib.ndpointer(dtype=np.float64),     # time_proceed
            np.ctypeslib.ndpointer(dtype=np.float64),     # t_update_in
            np.ctypeslib.ndpointer(dtype=np.int32) ]      # number of threads
        self.settings.restype = c_void_p       

        #海の条件(mesh.csv)を読み込み初期設定を行う
        self.setcond = self.libs.setcond
        self.setcond.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.int32)        # 初期条件ファイルの読み込みスイッチ          
        ] 
        self.setcond.restype = c_void_p       
        
        self.receive_initial_condition = self.libs.receive_initial_condition
        self.receive_initial_condition.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.int32) ,      # i_no
            np.ctypeslib.ndpointer(dtype=np.int32) ,      # j_no
            np.ctypeslib.ndpointer(dtype=np.float64),     # dy
            np.ctypeslib.ndpointer(dtype=np.float64)      # dx
        ] 
        self.receive_initial_condition.restype = c_void_p 

        self.return_mesh_num = self.libs.return_mesh_num
        self.return_mesh_num.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.int32) ,      # i_no
            np.ctypeslib.ndpointer(dtype=np.int32)        # j_no
        ]
        self.return_mesh_num.restype = c_void_p  

        # 波の計算
        self.calc_sea = self.libs.calc_sea
        self.calc_sea.argtypes = [ 
            np.ctypeslib.ndpointer(dtype=np.float64),     # 浸水深の格納先となる配列
            np.ctypeslib.ndpointer(dtype=np.int32)        # 流速を初期化するか
            ]     
        self.calc_sea.restype = c_void_p       

        # 波の計算 with 流速値の交換
        self.calc_sea_with_velocity = self.libs.calc_sea_with_velocity
        self.calc_sea_with_velocity.argtypes = [ 
            np.ctypeslib.ndpointer(dtype=np.float64),    # 浸水深の格納先となる配列
            np.ctypeslib.ndpointer(dtype=np.float64),    # vel_xiの格納先となる配列
            np.ctypeslib.ndpointer(dtype=np.float64)     # vel_etaの格納先となる配列
            ]     
        self.calc_sea_with_velocity.restype = c_void_p       


    # 計算時間間隔などの各種設定を行う
    def common_settings( 
            self , 
            dt , 
            time_proceed , 
            t_update , 
            n_thread 
        ):

        self.dt                 = np.array( dt           , dtype='float64' )  # 変微分方程式のΔt
        self.time_proceed       = np.array( time_proceed , dtype='float64' )  # 一度の計算で進める秒数
        self.t_update           = np.array( t_update     , dtype='float64' )  # 一コマの秒数
        self.number_of_thread   = np.array( n_thread     , dtype='int32' )    # 一コマの秒数

        # SharedObjectに共有
        self.settings( 
            self.dt , 
            self.time_proceed ,
            self.t_update , 
            self.number_of_thread 
        ) 

    # 初期条件などの設定を行う
    def initialize_condition(self , 
            initial_depth: float , 
            read_initial_file: bool = True ,
            i_no: int = None , 
            j_no: int = None , 
            dy: float = None , 
            dx: float = None 
        ): 

        # 初期ファイルの読み込み有り無し
        if read_initial_file == True:
            read_initial_file_tmp   = np.array( [1] , dtype='int32')
        else: 
            read_initial_file_tmp   = np.array( [0] , dtype='int32')

            if any( [i_no,j_no,dy,dx] )  == None:
                raise KeyError(f"must fill the i,j,dx,dy values")
            i_no_tmp   = np.array( i_no , dtype='int32' )
            j_no_tmp   = np.array( j_no , dtype='int32' )
            dy_tmp     = np.array( dy   , dtype='float64' ) 
            dx_tmp     = np.array( dx   , dtype='float64' ) 


            self.receive_initial_condition( i_no_tmp , j_no_tmp , dy_tmp, dx_tmp )

        # 初期設定を実施
        self.setcond( read_initial_file_tmp ) 

        # メッシュ数の設定値を取得
        i_no = np.array( [1] , dtype='int32') # dmy
        j_no = np.array( [1] , dtype='int32') # dmy
        self.return_mesh_num( i_no, j_no )    

        self.t_no = int( self.time_proceed / self.t_update ) 

        # 浸水データのマトリックス
        self.h_time_sr = np.ones( ( self.t_no , j_no[0] , i_no[0] ) , dtype='float64') * initial_depth 
        
        # 流速関係のマトリックス
        self.u_time_sr = np.zeros( ( self.t_no , j_no[0] , i_no[0] ) , dtype='float64') 
        self.v_time_sr = np.zeros( ( self.t_no , j_no[0] , i_no[0] ) , dtype='float64')
        
        # インスタンス変数に設定
        self.i_no = i_no[0]
        self.j_no = j_no[0]


    # 波の計算を継続する（流速の初期化等は行わない）
    def proceed_sea_field_ordinary( self ):
        self.calc_sea( self.h_time_sr ,np.array( [0] , dtype='int32') )

    # 流速の初期化を行った上で、波の計算を継続する
    def proceed_sea_field__with_initialize_velocity( self ):
        self.calc_sea( self.h_time_sr ,np.array( [1] , dtype='int32')  )
        
    # 流速も相互にやり取りする
    def proceed_sea_field_with_velocity( self ):
        self.calc_sea_with_velocity( self.h_time_sr , self.u_time_sr , self.v_time_sr )

    def allocate_reshaped_res(self):
        # dmdコードでの計算結果と配列を合わせておく
        print( "Reshaped" , np.shape( self.h_time_sr ) )
        self.h_time_sr_reshaped = np.transpose( self.h_time_sr , (2, 1, 0) )
        
