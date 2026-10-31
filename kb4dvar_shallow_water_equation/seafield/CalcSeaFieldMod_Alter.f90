! Rev 2025/02/24
! 適当に作った非線形長波理論のプログラム
! pythonで動くsharedobjectとして設計

Module globals 
    implicit none 

    !python側に共有する内容
    real(8),allocatable,dimension(:,:,:):: depth_timesr             ! python側に共有する
    !real(8),allocatable,dimension(:,:,:):: h_time_sr                ! python側と共有

    real(8),allocatable,dimension(:,:)::x, y ,z ,zb                 !cartesian coordinate 
    real(8),allocatable,dimension(:,:)::xi, eta                     !global coordinate
    real(8),allocatable,dimension(:,:)::d_xi, d_eta                 !differential 
    real(8),allocatable,dimension(:,:)::dxi_dx, dxi_dy              !global coordinate
    real(8),allocatable,dimension(:,:)::deta_dx, deta_dy            !global coordinate
    real(8),allocatable,dimension(:,:)::dx_dxi, dx_deta             ! coordinate
    real(8),allocatable,dimension(:,:)::dy_dxi, dy_deta             ! coordinate

    real(8),allocatable,dimension(:,:)::a1, a2 ,a3 ,a4 ,a5 ,a6      !cofficient 
    real(8),allocatable,dimension(:,:)::h, nn ,D                    !cartesian coordinate 
    real(8),allocatable,dimension(:,:)::Jac                         !det Jacobian 

    real(8),allocatable,dimension(:,:)::u, v                     !velocity
    real(8),allocatable,dimension(:,:)::xi_u, eta_u              !velocity

    ! integer,parameter::i_no = 200, j_no = 200
    integer:: i_no , j_no 
    integer :: t_loops , t_save , t_no
    real(8),parameter::h_init_val = 10.0  
    real(8),parameter::gg = 9.80665 
    real(8),parameter::n_cof = 2.0
    real(8),parameter::reflect = 0.5
    real(8),parameter::reduction = 0.0001
    real(8)::dt , time_proceed , t_update 

    integer,parameter::i_target = 95 , j_target = 55
    integer::parallel_threads
    integer:: condition_check 

End Module globals

Module SeaField
    !$ use omp_lib
    use globals
    !$ use,intrinsic :: iso_fortran_env    
    implicit none 
    integer:: i , j , k , t     

contains 

    !mesh.csvの情報を設定
    subroutine setcond( read_case_sw ) bind(C)
        use globals
        implicit none
        integer,intent(in):: read_case_sw 
        integer::i,j,k
        integer:: ii , jj 
        real(8)::dx , dy 
        real::rnd_n 

        if( read_case_sw > 0 ) then 
            call read_case
        elseif( condition_check /= 1 ) then  
            write(*,*) "must use receive_initial_condition before..."
        end if 

        ii = i_no
        jj = j_no

        allocate( d_xi(i_no,j_no) , d_eta(i_no,j_no) ) 
        allocate( dxi_dx(i_no,j_no) ,dxi_dy(i_no,j_no) , deta_dx(i_no,j_no) ,deta_dy(i_no,j_no) ) 
        allocate( dx_dxi(i_no,j_no) ,dx_deta(i_no,j_no) , dy_dxi(i_no,j_no) ,dy_deta(i_no,j_no) ) 
        allocate( Jac(i_no,j_no) ) 
        allocate( a1(i_no,j_no) , a2(i_no,j_no) , a3(i_no,j_no) , a4(i_no,j_no) , a5(i_no,j_no) , a6(i_no,j_no) ) 
        allocate( u(i_no,j_no) , v(i_no,j_no) ) 
        allocate( xi_u(i_no,j_no) , eta_u(i_no,j_no) )         

        d_xi(:,:) = 1.0 ; d_eta(:,:) = 1.0 
        dxi_dx(:,:) = 0.0 ; dxi_dy(:,:) = 0.0 
        deta_dx(:,:) = 0.0 ; deta_dy(:,:) = 0.0 
        dx_dxi(:,:) = 0.0 ; dx_deta(:,:) = 0.0 
        dy_dxi(:,:) = 0.0 ; dy_deta(:,:) = 0.0 
        Jac(:,:) = 0.0 

        do i = 2,  ii - 1
        do j = 2 , jj - 1

            dx_dxi(i,j)  = ( x(i,j+1) - x(i,j-1) ) / ( d_xi(i,j)  + d_xi(i,j-1) )
            dx_deta(i,j) = ( x(i+1,j) - x(i-1,j) ) / ( d_eta(i,j) + d_eta(i-1,j) )
            dy_dxi(i,j)  = ( y(i,j+1) - y(i,j-1) ) / ( d_xi(i,j)  + d_xi(i,j-1) )
            dy_deta(i,j) = ( y(i+1,j) - y(i-1,j) ) / ( d_eta(i,j) + d_eta(i-1,j) )
            
            Jac(i,j)     = 1.0 / ( dx_dxi(i,j) * dy_deta(i,j) - dx_deta(i,j) * dy_dxi(i,j) )

            dxi_dx(i,j)  = Jac(i,j) * dy_deta(i,j)
            dxi_dy(i,j)  = -1.0 * Jac(i,j) * dx_deta(i,j)
            deta_dx(i,j) = -1.0 *Jac(i,j) * dy_dxi(i,j)
            deta_dy(i,j) = Jac(i,j) * dx_dxi(i,j)
        end do 
        end do 


        do i = 2 , ii - 1
        do j = 2 , jj - 1            
            a1(i,j) = dxi_dx(i,j)  &
                    & * ( dx_dxi(i,j+1) - dx_dxi(i,j) ) / d_xi(i,j)  &
                    & + dxi_dy(i,j) &
                    & * ( dy_dxi(i,j+1) - dy_dxi(i,j) ) / d_xi(i,j) 

            a2(i,j) = 2.0 * ( dxi_dx(i,j) &
                    & * ( dx_dxi(i+1,j) - dx_dxi(i,j) ) / d_eta(i,j) & 
                    & + dxi_dy(i,j) & 
                    & * ( dy_dxi(i+1,j) - dy_dxi(i,j) ) / d_eta(i,j) )

            a3(i,j) = dxi_dx(i,j) &
                    & * ( dx_deta(i+1,j) - dx_deta(i,j) ) / d_eta(i,j)   & 
                    & + dxi_dy(i,j) &
                    & * ( dy_deta(i+1,j) - dy_deta(i,j) ) / d_eta(i,j)  

            a4(i,j) = deta_dx(i,j) &
                    & * ( dx_dxi(i,j+1) - dx_dxi(i,j) ) / d_xi(i,j)   &
                    & + deta_dy(i,j) &
                    & * ( dy_dxi(i,j+1) - dy_dxi(i,j) ) / d_xi(i,j) 

            a5(i,j) = 2.0 * ( deta_dx(i,j) &
                    & * ( dx_dxi(i+1,j) - dx_dxi(i,j) ) / d_eta(i,j) & 
                    & + deta_dy(i,j) &
                    & * ( dy_dxi(i+1,j) - dy_dxi(i,j) ) / d_eta(i,j) )

            a6(i,j) = deta_dx(i,j) &
                    & * ( dx_deta(i+1,j) - dx_deta(i,j) ) / d_eta(i,j) & 
                    & + deta_dy(i,j) & 
                    & * ( dy_deta(i+1,j) - dy_deta(i,j) ) / d_eta(i,j) 

        end do 
        end do 

        ! 流速関連:u ,v ,xi_u, eta_u の初期化
        call initialize_velocity 

    end subroutine setcond

    subroutine initialize_velocity() 
        use globals 
        implicit none 

        u(:,:) = 0.0 
        v(:,:) = 0.0 
        xi_u(:,:) = 0.0 
        eta_u(:,:) = 0.0 

    end subroutine initialize_velocity

    subroutine return_mesh_num(i_no_tmp , j_no_tmp ) bind(C)
        use globals 
        implicit none 
        integer,intent(inout):: i_no_tmp , j_no_tmp

        i_no_tmp = i_no 
        j_no_tmp = j_no 

    end subroutine return_mesh_num

    subroutine read_case
        ! meshの内容を読み込む
        use globals
        implicit none
        integer::i,j,k , idmy
        integer:: ii , jj , i_no_tmp, j_no_tmp
        real(8)::dx , dy 
        real::rnd_n 

        OPEN(11, FILE='./seafield/mesh.csv' , status='old')
        read(11,*) i_no_tmp , j_no_tmp  !header
        i_no = i_no_tmp
        j_no = j_no_tmp

        allocate( x(i_no,j_no) , y(i_no,j_no) , z(i_no,j_no) , zb( i_no-1 , j_no-1 ) )
        allocate( D(i_no,j_no) , h(i_no,j_no) , nn(i_no,j_no) )

        !set x, y, z 
        do i = 1 , i_no 
            do j = 1 , j_no

                read(11,*) idmy, x(i,j), y(i,j) , z(i,j)                 
                h(i,j) = h_init_val

                nn(:,:) = 0.0
                D(i,j) = h(i,j) + nn(i,j) 

            end do 
        end do 

        do i = 2 , i_no - 1 ; do j = 2 , j_no - 1 
            zb(i,j) = ( z(i,j) + z(i+1,j) + z(i,j+1) + z(i+1, j+1) ) * 0.25 
        end do  ; end do 
        close(11)

        condition_check = 1 

    end subroutine read_case

    subroutine receive_initial_condition( i_no_in , j_no_in, dy_in, dx_in ) bind(C)
        ! python側から初期条件を受け取る場合（read_caseは不使用）
        use globals
        implicit none
        integer,intent(in):: i_no_in , j_no_in 
        real(8),intent(in):: dy_in , dx_in 


        integer::i,j,k , idmy
        integer:: ii , jj 
        real(8)::dx , dy 
        real::rnd_n 

        real(8),parameter:: z_int = 1.0

        i_no = i_no_in 
        j_no = j_no_in 

        allocate( x(i_no,j_no) , y(i_no,j_no) , z(i_no,j_no) , zb( i_no-1 , j_no-1 ) )
        allocate( D(i_no,j_no) , h(i_no,j_no) , nn(i_no,j_no) )

        !set x, y, z 
        x(1,1) = 0.0 
        y(1,1) = 0.0 
        
        do i = 1 , i_no 
            do j = 1 , j_no

                x(i,j) = dble(j) * dx_in + x(1,1)
                y(i,j) = dble(i) * dy_in + y(1,1)
                z(i,j) = z_int       
                h(i,j) = h_init_val

                nn(:,:) = 0.0
                D(i,j) = h(i,j) + nn(i,j) 

            end do 
        end do 

        do i = 2 , i_no - 1 ; do j = 2 , j_no - 1 
            zb(i,j) = ( z(i,j) + z(i+1,j) + z(i,j+1) + z(i+1, j+1) ) * 0.25 
        end do  ; end do 
    
        condition_check = 1 

    end subroutine receive_initial_condition


    ! その他の設定
    ! dt_in:氾濫解析の更新時間、time_proceed_in:一度の計算で進めるβ秒数、t_update_in:アップデート間隔
    subroutine settings( dt_in , time_proceed_in , t_update_in, threads_in ) bind(C)
        use globals 
        implicit none
        real(8),intent(in):: dt_in , time_proceed_in , t_update_in 
        integer,intent(in):: threads_in  

        dt = dt_in 
        time_proceed = time_proceed_in 
        t_update = t_update_in 

        ! 計算用の各種設定値
        t_save = int( t_update / dt )
        t_loops = int( time_proceed / dt ) 
        t_no   = int( time_proceed / t_update )
        
        ! 並列スレッド数
        parallel_threads = threads_in
        if (parallel_threads > 1 ) then 
            !$ call omp_set_num_threads(parallel_threads)        
            !$ write(*,*) "threads;" , parallel_threads
        end if 

        ! 追記＝＞一旦python側から配列型と合わせて明示的に受け取るようにする。
        ! 計算速度に課題があれば、再構築
        ! allocate( h_time_sr(i_no , j_no , t_no ) )

    end subroutine settings 

    ! 水深分布を設定(未実装)
    ! 事前にsetcond, settingsが必要
    subroutine initialize_depth_mat( h_time_sr ) bind(C)
        use globals
    
        implicit none 
        real(8),dimension(i_no,j_no,t_no),intent(inout)::h_time_sr
        
        h_time_sr(:,:,:) = 0.0 
        
    end subroutine initialize_depth_mat

    ! 流れの計算を実行: 初回計算前には、必ずsetcondを回しておくこと(継続計算時は不要)
    subroutine calc_sea(h_time_sr , init_velocity ) bind(C)
        !$ use omp_lib
        use globals
        implicit none 
        real(8),dimension(i_no,j_no,t_no),intent(inout)::h_time_sr
        integer,intent(in):: init_velocity
        integer::t , time_out , i, j  
        
        ! h_time_srの(:,:,1)は初期時刻として与えること
        D(:,:) = h_time_sr(:,:,1)
        nn(:,:) = D(:,:) - h(:,:) ! マニュアル上、D,h,nnを分けて管理

        ! 流速関連:u ,v ,xi_u, eta_u の初期化
        if( init_velocity == 1 ) then             
            call initialize_velocity 
        end if 

        time_out = 0 
        Do t = 1 , t_loops
            call MomentumEquation
            call ContinuousEquation 
            call ReductionSeaWave

        ! update 
            if( mod(t , t_save) == 0 ) then 
                time_out = time_out + 1 
                h_time_sr(:,:,time_out) = D(:,:) !事前にinitialize_depth_matが必要
            end if 

        End Do 
    End subroutine calc_sea


    ! 流れの計算を実行: 初回計算前には、必ずsetcondを回しておくこと(継続計算時は不要)
    subroutine calc_sea_with_velocity(h_time_sr , u_time_sr , v_time_sr ) bind(C)
        !$ use omp_lib
        use globals
        implicit none 
        real(8),dimension(i_no,j_no,t_no),intent(inout)::h_time_sr,u_time_sr,v_time_sr
        integer::t , time_out , i, j  
        
        ! h_time_srの(:,:,1)は初期時刻として与えること
        D(:,:) = h_time_sr(:,:,1)
        u(:,:) = u_time_sr(:,:,1)
        v(:,:) = v_time_sr(:,:,1)

        ! 流速の設定
        call set_xiu_etau_val

        nn(:,:) = D(:,:) - h(:,:) ! マニュアル上、D,h,nnを分けて管理

        time_out = 0 
        Do t = 1 , t_loops
            call MomentumEquation
            call ContinuousEquation 
            call ReductionSeaWave

            ! update 
            if( mod(t , t_save) == 0 ) then 
                time_out = time_out + 1 
                h_time_sr(:,:,time_out) = D(:,:) !事前にinitialize_depth_matが必要

                call set_u_v_value
                u_time_sr(:,:,time_out) = xi_u(:,:)
                v_time_sr(:,:,time_out) = eta_u(:,:)
            end if 
        End Do 

    End subroutine calc_sea_with_velocity


    ! 一般座標系の流速から物理座標の流速を算出
    subroutine set_u_v_value
        use globals 
        implicit none 

        !$omp parallel do private(j)
        do i = 1, i_no 
            do j = 1 , j_no 
                u(i,j)=1.0/Jac(i,j)*(deta_dy(i,j)*xi_u(i,j)-dxi_dy(i,j)*eta_u(i,j)) 
                v(i,j)=1.0/Jac(i,j)*(-1.0*deta_dx(i,j)*xi_u(i,j)+ dxi_dx(i,j)*eta_u(i,j))
            end do 
        end do 
        !$omp end parallel do
    end subroutine 


    ! 物理座標系から一般座標系の流速を設定
    subroutine set_xiu_etau_val
        use globals 
        implicit none 

        !$omp parallel do private(j)
        do i = 1, i_no 
            do j = 1 , j_no 
                xi_u(i,j)=dxi_dx(i,j)*u(i,j)+dxi_dy(i,j)*v(i,j)
                eta_u(i,j)=deta_dx(i,j)*u(i,j)+deta_dy(i,j)*v(i,j)
            end do 
        end do 
        !$omp end parallel do

    end subroutine set_xiu_etau_val


    ! 計算を続けているとノイズが蓄積してNanが発生するので元の状態に段々と戻る処理を入れる
    subroutine ReductionSeaWave
        !$ use omp_lib
        use globals
        implicit none 

        D(:,:) = D(:,:) * ( 1.0 - reduction ) + h(:,:) *( reduction )
        nn(:,:) = nn(:,:) * ( 1.0 - reduction )
        xi_u(:,:) = xi_u(:,:) * (1.0 - reduction)
        eta_u(:,:) = eta_u(:,:) * (1.0 - reduction)

    end subroutine ReductionSeaWave

    !Calc MomentumEQ
    subroutine MomentumEquation
        use globals
        implicit none        

        real(8)::sec, thir , forth , fifth , sixth , seventh ,eighth
        real(8),dimension( i_no, j_no ) ::xi_u_next , eta_u_next 
        

        !xi -direction 
        !$omp parallel do private(j,sec, thir, forth , fifth , sixth , seventh , eighth ) !shared( xi_u_next )        
        do i = 2 , i_no - 1 
            do j = 2 , j_no - 1 
                sec  = xi_u(i,j) * ( xi_u( i , j+1 ) - xi_u(i, j-1) )  / ( d_xi(i,j) + d_xi(i,j-1) ) 
                thir = eta_u(i,j) * ( xi_u( i+1 , j ) - xi_u(i-1, j) ) / ( d_eta(i,j) + d_eta(i-1,j) ) 
                
                forth = a1(i,j) * xi_u(i,j) * xi_u(i,j)
                fifth = a2(i,j) * xi_u(i,j) * eta_u(i,j)             
                sixth = a3(i,j) * eta_u(i,j) * eta_u(i,j)                

                seventh = gg * ( (dxi_dx(i,j) ** 2.0 + dxi_dy(i,j) ** 2.0 ) &
                        & * ( nn(i,j+1) - nn(i,j-1) ) / ( d_xi(i,j) + d_xi(i,j-1) ) & 
                        & + ( dxi_dx(i,j) * deta_dx(i,j) + dxi_dy(i,j) * deta_dy(i,j) ) & 
                        & * ( nn(i+1,j) - nn(i-1,j) )   / ( d_eta(i,j) + d_eta(i-1,j) ) )                         

                eighth = gg * xi_u(i,j) * n_cof ** 2.0 / ( D(i,j) ** ( 10.0 / 3.0 ) )  / Jac(i,j) *  & 
                        & ( (deta_dy(i,j) * xi_u(i,j) - dxi_dy(i,j) * eta_u(i,j)) ** 2.0 + & 
                        &   ( -1.0 * deta_dx(i,j) * xi_u(i,j) + dxi_dx(i,j) * eta_u(i,j) ) ** 2.0 ) ** 0.5 
                
                !4th~6th able to be omitted 
                xi_u_next( i , j ) = -1.0 * ( sec + thir + forth + fifth + sixth + seventh  + eighth ) * dt 

            end do 
        end do 
        !$omp end parallel do

        
        !eta -direction 
        !$omp parallel do private(j,sec, thir, forth , fifth , sixth , seventh , eighth ) !shared( eta_u_next )        
        do i = 2 , i_no - 1 
            do j = 2 , j_no - 1 
                sec  = xi_u(i,j)  * ( eta_u( i , j+1 ) - eta_u(i, j-1) ) / ( d_xi(i,j)  + d_xi(i,j-1)  ) 
                thir = eta_u(i,j) * ( eta_u( i+1 , j ) - eta_u(i-1, j) ) / ( d_eta(i,j) + d_eta(i-1,j) ) 
                forth = a4(i,j) * xi_u(i,j) * xi_u(i,j)
                fifth = a5(i,j) * xi_u(i,j) * eta_u(i,j)
                sixth = a6(i,j) * eta_u(i,j) * eta_u(i,j)

                seventh = gg * ( (deta_dx(i,j) ** 2.0 + deta_dy(i,j) ** 2.0 ) & 
                        & * ( nn(i+1,j) - nn(i-1,j) ) / ( d_eta(i,j) + d_eta(i-1,j) )  & 
                        & + ( deta_dx(i,j) * dxi_dx(i,j) + deta_dy(i,j) * dxi_dy(i,j) ) & 
                        & * ( nn(i,j+1) - nn(i,j-1) ) / ( d_xi(i,j) + d_xi(i,j-1) ) ) 

                eighth = gg * eta_u(i,j) * n_cof ** 2.0 / ( D(i,j) ** (10.0 / 3.0 ) )  / Jac(i,j) *  & 
                        & ( (deta_dy(i,j) * xi_u(i,j) - dxi_dy(i,j) * eta_u(i,j)) ** 2.0 + & 
                        &   ( -1.0 * deta_dx(i,j) * xi_u(i,j) + dxi_dx(i,j) * eta_u(i,j) ) ** 2.0 ) ** 0.5 

                !4th~6th able to be omitted 
                eta_u_next( i , j ) = -1.0 * ( sec + thir + forth + fifth + sixth + seventh + eighth ) * dt 

            end do 
        end do 
        !$omp end parallel do


        !update u value
        do i = 2 , i_no - 1 
            do j = 2 , j_no - 1 
                xi_u(i,j)  = xi_u(i,j)  +  xi_u_next(i,j)
                eta_u(i,j) = eta_u(i,j) +  eta_u_next(i,j) 

            end do 
        end do 

        !boundary 
        do i = 1 , i_no 
            xi_u(i,1) =   ( xi_u(i,2) + xi_u(i,3) ) * (-0.5) * reflect
            eta_u(i,1) =  ( eta_u(i,2) + eta_u(i,3) ) * 0.5 

            xi_u(i,j_no)  = ( xi_u(i,j_no-1) + xi_u(i,j_no-2) ) * (-0.5) * reflect
            eta_u(i,j_no) = ( eta_u(i,j_no-1) + eta_u(i,j_no-2) ) * 0.5 

        end do 

        do j = 1, j_no 
            xi_u(1,j)  = ( xi_u(2,j) + xi_u(3,j) ) * 0.5 
            eta_u(1,j) = ( eta_u(2,j) + eta_u(3,j) ) * (-0.5) * reflect

            xi_u(i_no,j)  = ( xi_u(i_no-1,j) + xi_u(i_no-2,j) ) * 0.5  
            eta_u(i_no,j) = ( eta_u(i_no-1,j) + eta_u(i_no-2,j) ) * (-0.5) * reflect
        end do 

    End Subroutine MomentumEquation

    !Calc ContinuousEQ
    Subroutine ContinuousEquation 
        use globals
        implicit none 
        real(8)::sec , thir ,dh
        real(8),dimension(i_no, j_no)::sec_pre , thir_pre 

        do i = 1, i_no 
        do j = 1, j_no 
            D(i,j) = h(i,j) + nn(i,j)
        end do 
        end do 

        do i = 1 , i_no 
            do j = 1 , j_no 
                sec_pre(i,j) = D(i,j) * xi_u(i,j)  / Jac(i,j) 
                thir_pre(i,j) = D(i,j) * eta_u(i,j)  / Jac(i,j) 

            end do 
        end do 

        !boundary 
        do i = 1 , i_no 
            sec_pre(i,1)  = ( sec_pre(i,2) + sec_pre(i,3)) * (-0.5)
            thir_pre(i,1) = ( thir_pre(i,2) + thir_pre(i,3) ) * (-0.5)

            sec_pre(i,j_no)  = ( sec_pre(i,j_no-1) + sec_pre(i,j_no-2) ) * (-0.5)
            thir_pre(i,j_no) = ( thir_pre(i,j_no-1) + thir_pre(i,j_no-2) ) * (-0.5)

        end do 
        do j = 1, j_no 
            sec_pre(1,j)  = ( sec_pre(2,j) + sec_pre(3,j)) * (-0.5)
            thir_pre(1,j) = ( thir_pre(2,j) + thir_pre(3,j) ) * (-0.5)

            sec_pre(i_no,j)  = ( sec_pre(i_no-1,j) + sec_pre(i_no-2,j) ) * (-0.5)
            thir_pre(i_no,j) = ( thir_pre(i_no-1,j) + thir_pre(i_no-2,j) ) * (-0.5)
        end do 

        !bound_other
            sec_pre(1,1)  = D(2,2) * xi_u(2,2)   / Jac(2,2) * (-1.0)
            thir_pre(1,1) = D(2,2) * eta_u(2,2)  / Jac(2,2) * (-1.0)

            sec_pre(1,j_no)  = D(2,j_no-1) * xi_u(2,j_no-1)   / Jac(2,j_no-1) * (-1.0)
            thir_pre(1,j_no) = D(2,j_no-1) * eta_u(2,j_no-1)  / Jac(2,j_no-1) * (-1.0)

            sec_pre(i_no,1)  = D(i_no-1,2) * xi_u(i_no-1,2)   / Jac(i_no-1,2) * (-1.0)
            thir_pre(i_no,1) = D(i_no-1,2) * eta_u(i_no-1,2)  / Jac(i_no-1,2) * (-1.0)

            sec_pre(i_no,j_no)  = D(i_no-1,j_no-1) * xi_u(i_no-1,j_no-1)   / Jac(i_no-1,j_no-1) * (-1.0)
            thir_pre(i_no,j_no) = D(i_no-1,j_no-1) * eta_u(i_no-1,j_no-1)  / Jac(i_no-1,j_no-1) * (-1.0)


        !Calc ContinuousEQ
        !$omp parallel do private(j, sec , thir , dh) !shared( D, nn ) 
        do i = 2, i_no - 1 
            do j = 2 , j_no -1 
                sec  = ( sec_pre(i,j+1) - sec_pre(i,j-1) ) / ( d_xi(i,j) + d_xi(i,j-1)) 
                thir = ( thir_pre(i+1 , j ) - thir_pre(i-1,j) ) / ( d_eta(i,j) + d_eta(i-1,j) ) 

                dh = ( sec + thir ) * Jac(i,j) * dt  * (-1.0)
                D(i,j) = D(i,j) + dh 
                nn(i,j) = D(i,j) - h(i,j)

            end do 
        end do 
        !$omp end parallel do        

        !boundary 
        do i = 1 , i_no 
            j = 1
                D(i,j)  = ( D(i,j+2) + D(i,j) )  * 0.5
                nn(i,j) = D(i,j) - h(i,j)

            j = j_no
                D(i,j)  = ( D(i,j-2) + D(i,j) )  * 0.5
                nn(i,j) = D(i,j) - h(i,j)
        end do 

        do j = 1, j_no 
            i = 1 
                D(i,j) = ( D(i+2, j ) + D(i,j) ) * 0.5 
                nn(i,j) = D(i,j) - h(i,j)

            i = i_no
                D(i,j)  = ( D(i-2,j) + D(i,j) )  * 0.5
                nn(i,j) = D(i,j) - h(i,j)
        end do         

        !bound_other
        i = 1 ; j = 1
            D(i,j)  = ( D(i+2,j+2) + D(i,j) )  * 0.5
            nn(i,j) = D(i,j) - h(i,j)

        i = 1 ; j = j_no
            D(i,j)  = ( D(i+2,j-2) + D(i,j) )  * 0.5
            nn(i,j) = D(i,j) - h(i,j)            

        i = i_no ; j = 1
            D(i,j)  = ( D(i-2,j+2) + D(i,j) )  * 0.5
            nn(i,j) = D(i,j) - h(i,j)                        

        i = i_no ; j = j_no 
            D(i,j)  = ( D(i-2,j-2) + D(i,j) )  * 0.5
            nn(i,j) = D(i,j) - h(i,j)

    End Subroutine ContinuousEquation 

End Module SeaField