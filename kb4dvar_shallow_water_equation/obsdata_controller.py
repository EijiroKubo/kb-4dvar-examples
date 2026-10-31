import os 
import zipfile
import numpy as np 
import pandas as pd 


class observation:

    def __init__(self): 
        
        x_tmp = np.arange( 0 , 200 , 4 ) + 2
        y_tmp = np.arange( 0 , 200 , 4 ) + 2 
        
        data_tmp = [] 
        for j in x_tmp:
            for i in y_tmp: 
                data_tmp.append( [j,i] )
            
        self.data = np.array( data_tmp )
        
        '''
        self.data = np.array( [
            [ 10 , 10 ] , 
            [ 30 , 10 ] ,
            [ 30 , 30 ]               
        ] ) 
        '''
        




class ObsDataProcessor:    

    def __init__( self ):
        self.obs_pass = os.path.join( os.getcwd() , "observe" ) 


    def set_obs_place(self , obs_point ):
        self.obs_point = obs_point 


    def obs_function( self, X ): 
        ''' X[i,j,t]'''
        x_temp = [] 
        for i in range( np.shape(self.obs_point)[0]):
            # Y , X , T 
            x_temp.append( X[int(self.obs_point[i][1]),int(self.obs_point[i][0]),:] )

        x_temp = np.array( x_temp )
        x_temp = np.where( x_temp < 0 , 0 , x_temp )
        return x_temp 
    

    def obs_function_for_adj( self, X ): 
        ''' 観測行列のマトリクスを生成する 
        Input: X ( t, ny, nx )
        Output: ( ny ,nx ,t )
        '''
        
        t , ny ,nx = np.shape( X )
        x_temp = np.zeros( ( ny , nx , t ) )
    
        
        for i in range( np.shape(self.obs_point)[0]):
            ii = int(self.obs_point[i][1].item())  # Y
            jj = int(self.obs_point[i][0].item())  # X 
            x_temp[ii,jj,:] = X[:,ii,jj].real
                    
        return x_temp 


    def set_obsmatrix_for_adj(self,ny, nx, nt ):
        ''' 実観測データを観測行列の形に整形しておく '''

        self.obs_mat = np.zeros( (ny,nx,nt) ) 
        iter = nt

        for i in range( np.shape(self.obs_point)[0]):
            ii = int(self.obs_point[i][1].item())   # Y 
            jj = int(self.obs_point[i][0].item())   # X 
            self.obs_mat[ii,jj,:] = self.obs[ i ][:iter]


    def create_test_observation( self, X ):
        obs = self.obs_function( X )
        
        t = 0 
        for i in range( np.shape( self.obs_point )[0]): 
            f_name = "observations_" + str( i ) + ".csv"
            f = os.path.join( self.obs_pass, f_name )
            df = pd.DataFrame( obs[i] )
            df.to_csv( f , sep= "," , index = None , header = None )

            f_place_name = "obs_place_" + str( i ) + ".csv" 
            f_place = os.path.join( self.obs_pass , f_place_name )
            d = np.array([int(self.obs_point[i][0]),int(self.obs_point[i][1])])
            df_place = pd.DataFrame( d )
            df_place.to_csv( f_place , sep="," , index = None , header = None )
        
            
    def read_observation(self):
        ''' 観測データを取得 '''
        obs_zip = os.path.join( os.getcwd() , "observe.zip")
        self.obs = [] 
        obs_point_temp = [] 
        with zipfile.ZipFile( obs_zip ) as archive:
            self.obs_num = sum(
                name.startswith("observations_") for name in archive.namelist()
            )
            for i in range( self.obs_num ):
                with archive.open("observations_" + str(i) + ".csv") as f_obs:
                    df = pd.read_csv( f_obs , header=None )
                d_temp = np.array( df )
                self.obs.append( d_temp )

                with archive.open("obs_place_" + str(i) + ".csv") as f_obs_point:
                    df_point = pd.read_csv( f_obs_point , header = None )
                d_point_temp = [ df_point.iloc[0].values ,df_point.iloc[1].values ]
                obs_point_temp.append( d_point_temp )


        self.set_obs_place( obs_point_temp )            

        _step_temp = np.shape( d_temp )[0]
        self.obs = np.reshape( np.array( self.obs ) , ( self.obs_num , _step_temp ) ) 
            

