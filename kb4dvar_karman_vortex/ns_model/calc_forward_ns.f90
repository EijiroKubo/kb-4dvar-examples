
! どうでも良いがカルマンフィルタはkalman, カルマン渦はkarmanである。
Module KarmanVortexField
  ! どうやらこれを入れておいた方が良さそう
  use mod_variables
  implicit none

contains 

  subroutine calc_karman_vortex( qin_time_sr ) bind(C)
    ! sub_4dvar.f90から前進計算部分のみ抽出した      
    use mod_variables
    implicit none  
    
    real(8),dimension(imx1, jmx1, 3, nsall),intent(inout)::qin_time_sr

    ustg(1:imx1,1:jmax) = qin_time_sr(1:imx1,1:jmax,1,1)
    vstg(1:imax,1:jmx1) = qin_time_sr(1:imax,1:jmx1,2,1)
    pcnt(1:imax,1:jmax) = qin_time_sr(1:imax,1:jmax,3,1)

      !ustg(1:imx1,1:jmax) = us_time_sr(:,:,1)
      !vstg(1:imax,1:jmx1) = vs_time_sr(:,:,1)
      !pcnt(1:imax,1:jmax) = p_time_sr(:,:,1)

    do istep=2,nsall
      call sub_bc_outer (ustg,vstg,pcnt)
      call sub_bc_wall  (ustg,vstg,pcnt)
      call sub_rhs_fwd
      call sub_HSMAC_fwd(ustg,vstg,pcnt)

      qin_time_sr(1:imx1,1:jmax,1,istep) = ustg(1:imx1,1:jmax)
      qin_time_sr(1:imax,1:jmx1,2,istep) = vstg(1:imax,1:jmx1) 
      qin_time_sr(1:imax,1:jmax,3,istep) = pcnt(1:imax,1:jmax) 

    enddo
  
  end subroutine calc_karman_vortex


  ! ==============================================
  ! 以下, 桑名自作プログラムのIFに合わせてモジュール追加
  ! ==============================================

  subroutine read_initial_condition() bind(C)
    ! 正解データ(通常の4dvar計算)と同じ条件ファイルを読み込む
    use mod_variables 
    implicit none 

    open(28,file='./ns_model/4DVar.inp',form='formatted')
    read(28,*) modef               ! Mode (1:TLM&ADJ Check, 2:4DV, 3:4DV restart)
    read(28,*) n_prb               ! Problem 1:Karman vortex, 2:Vortex advection
    read(28,*) ncycl,ITRMAX        ! Number of 4DV cycles, iteration on each cycle
    read(28,*) vars%ref            ! Error variance of pseudo measurement
    read(28,*) vars%mes            ! Measurement error variance
    read(28,*) vars%bkg            ! Coefficit of background term
    read(28,*) iskip,jskip,tskip   ! Every iskip & jskip & tskip for measuremnt 
    read(28,*) mlft,mrht           ! Left and right of measurement area
    read(28,*) mbot,mtop           ! Bottom and top of measurement area
    read(28,*)
    read(28,*) nstep,nsall         ! DA window steps, total # of time steps
    read(28,*) Re                  ! Reynolds number 
    read(28,*) jmax                ! Number of mesh: jmax (imax x jmax, imax = 3*jmax)
    read(28,*) dt                  ! Time step
    read(28,*) irestart            ! 0:initial, 1:restart
    read(28,*) ioutput             ! 0:No output to screen,1:otherwise
    read(28,*) iskip_plot          ! Output plot3d interval
    read(28,*) iskip_hist          ! Output history interval
    close(28)

    ! Specify the number of points per one side of square cylinder
      jobj = int(0.100*jmax)  ! 8 points per one side when jmax = 80
    ! jobj = int(0.064*jmax)  ! jmax = 320: the original setting 2018/12/26
      iobj = jobj

    ! その他の検討
    imax = 3*jmax
    imx1 = imax+1
    jmx1 = jmax+1


  !-----------------------------------------------------------------------
  !-Initial setup
  !-----------------------------------------------------------------------
    call sub_initial
    if(irestart==1) then
      if(ioutput==1) write(*,*) 'Restart'
      call sub_restart
      tstep_ini = tstep
      ostep_ini = ostep
    else
      if(ioutput==1) write(*,*) 'Impulse start'
      tstep = 0.0
      ostep = 0
    endif

  end subroutine read_initial_condition

  subroutine return_mesh_num( imax_tmp, imx1_tmp, jmax_tmp,jmx1_tmp ,t_steps_tmp ) bind(C)
    use mod_variables
    implicit none 
    integer,intent(inout):: imax_tmp, imx1_tmp, jmax_tmp, jmx1_tmp , t_steps_tmp

    imax_tmp = imax 
    imx1_tmp = imx1
    jmax_tmp = jmax
    jmx1_tmp = jmx1
    t_steps_tmp  = nsall 

  end subroutine return_mesh_num

  ! 検証用：重要コード
  subroutine set_validation_case( qin_time_sr ) bind(C)
    ! 検証用に計算条件をずらす
    use mod_variables
    implicit none

    integer :: iref
    real(8),dimension(imx1, jmx1, 3, nsall),intent(inout)::qin_time_sr

    qin_time_sr(:,:,:,:) = 0.d0

    if(n_prb==2) then
      call sub_addvtx(-4.d0,0.d0)
    else

      ! 注意
      ! sub_4dvarの初期状態10000+180stepのずらしを検証用として設定
      do iref=1, 10000 + int(7.2d0/dt)/2
        call sub_bc_outer (ustg,vstg,pcnt)
        call sub_bc_wall  (ustg,vstg,pcnt)
        call sub_rhs_fwd
        call sub_HSMAC_fwd(ustg,vstg,pcnt)
      enddo
    endif

    qin_time_sr(1:imx1,1:jmax,1,1) = ustg(1:imx1,1:jmax)
    qin_time_sr(1:imax,1:jmx1,2,1) = vstg(1:imax,1:jmx1)
    qin_time_sr(1:imax,1:jmax,3,1) = pcnt(1:imax,1:jmax)

  end subroutine set_validation_case

End Module KarmanVortexField
