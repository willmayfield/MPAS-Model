#!/bin/csh

#-----------------------------------------------------------------------
# Purpose: create namelist templates for MPAS
#-----------------------------------------------------------------------

set DATE        = $DATE       # From driver
set FCST_RANGE  = $FCST_RANGE # From driver

# Need to specify defaults for variables defined in scripts other than driver.csh
if ( ! $?case_number          )  set case_number = 7
if ( ! $?config_stop_time     )  set config_stop_time = `${TOOL_DIR}/da_advance_time.exe ${DATE} 0 -w`
if ( ! $?config_static_interp )  set config_static_interp = .false.
if ( ! $?config_vertical_grid )  set config_vertical_grid = .false.
if ( ! $?config_met_interp    )  set config_met_interp    = .false.
if ( ! $?config_input_sst     )  set config_input_sst     = .false.
if ( ! $?config_frac_seaice   )  set config_frac_seaice   = .false.
if ( ! $?config_input_name    )  set config_input_name    = 'dum'
if ( ! $?config_output_name   )  set config_output_name   = 'dum'
if ( ! $?update_sst_interval )   set update_sst_interval = none
if ( ! $?this_ungrib_vertical_levels ) set this_ungrib_vertical_levels = 32 #atj: changed these from "num_ungrib..." to "this_ungrib..." since diff for ic and lbc
if ( ! $?this_ungrib_soil_levels )     set this_ungrib_soil_levels = 4

#####

set START_DATE_MPAS = `${TOOL_DIR}/da_advance_time.exe $DATE 0 -w`
set END_DATE_MPAS   = `${TOOL_DIR}/da_advance_time.exe $DATE $FCST_RANGE -w`

if ( $update_sst_interval == none ) then
   set local_update_sst = .false.
else
   set local_update_sst = .true.
endif

if ( $MPAS_REGIONAL == true || $MPAS_REGIONAL == .true. ) then
   #set config_blend_bdy_terrain = .false. #turn off for static file generation (?)
   set config_blend_bdy_terrain = .true.
   set config_fg_interval = `expr $LBC_FREQ \* 3600`
   set config_apply_lbcs = .true.
else
   set config_blend_bdy_terrain = .false.
   set config_fg_interval = 86400
   set config_apply_lbcs = .false.
endif

#---------------------------------------
# Set PIO stride, etc
#  This is a moving target that depends on the number of processors used.
# $num_mpas_cells is set in driver.csh
#---------------------------------------
#atj: for now, turn off all options except what worked in my initial test
#if ( $num_mpas_cells == 6488066 ) then
#   set config_pio_num_iotasks = 36 # use this for 15-/3-km grid with 576 processors
#   set config_pio_stride      = 16
#else if ( $num_mpas_cells == 2621442 ) then
#   set config_pio_num_iotasks = 36  # global uniform 15-km mesh with 576 processors
#   set config_pio_stride      = 16
##else if ( $num_mpas_cells == 655362 ) then
##   set config_pio_num_iotasks = 1  # global uniform 30-km mesh
##   set config_pio_stride      = 1
#else
#atj moved below   set config_pio_num_iotasks = 0
#   set config_pio_stride      = 1
#endif

#---------------------------------------

echo "Filling Namelist for $1"

if ( $1 == "mpas" || $1 == "MPAS" ) then
   set config_pio_num_iotasks = $2 #atj
   set config_pio_stride      = $3 #atj
   goto MPAS
else if ( $1 == "mpas_init" || $1 == "MPAS_init" ) then
   set config_pio_num_iotasks = 0
   set config_pio_stride      = 1
   goto MPAS_init
else
   echo "wrong usage"
   echo "use as $0 (mpas,mpas_init)"
   exit
endif

#--------------------------------------------
# MPAS_initialization namelist.input
#--------------------------------------------

MPAS_init:

rm -f ./namelist.init_atmosphere
cat > ./namelist.init_atmosphere << EOF
&nhyd_model
   config_init_case = ${case_number} 
   config_start_time      = '${START_DATE_MPAS}'
   config_stop_time       = '${config_stop_time}'
   config_theta_adv_order = 3
   config_coef_3rd_order = 0.25
    config_interface_projection = 'linear_interpolation' !atj: following Jeff's scripts
/

&dimensions
   config_nvertlevels     = $num_mpas_vert_levels
   config_nsoillevels     = $num_mpas_soil_levels
   config_nfglevels       = $this_ungrib_vertical_levels
   config_nfgsoillevels   = $this_ungrib_soil_levels
   config_nsoilcat        = $num_soilcat 
/

&data_sources
   config_geog_data_path  = '${WPS_GEOG_DIR}/'
   config_met_prefix      = '${ungrib_prefx_model}'
   config_sfc_prefix      = '${ungrib_prefx_sst}'
   config_fg_interval     = $config_fg_interval
   config_landuse_data = 'MODIFIED_IGBP_MODIS_NOAH_15s'
   config_topo_data = 'GMTED2010'
   config_vegfrac_data = 'MODIS'
   config_albedo_data = 'MODIS'
   config_maxsnowalbedo_data = 'MODIS'
   config_supersample_factor = 12
   config_use_spechumd = false
   config_30s_supersample_factor = 3
   config_soilcat_data = "BNU"  
/

&vertical_grid
   config_ztop            = ${z_top_meters} 
   config_nsmterrain      = 1
   config_smooth_surfaces = .true.
   config_dzmin = 0.3
   config_nsm = 30
   config_tc_vertical_grid = .true.
   config_blend_bdy_terrain = ${config_blend_bdy_terrain}
   config_specified_zeta_levels = '/scratch2/BMC/fv3lam/HWT/CONUS_grid_static/L60.txt'
/

&interpolation_control
    config_extrap_airtemp = 'linear'
/


&preproc_stages
   config_static_interp   = $config_static_interp 
   config_native_gwd_static = $config_static_interp
   config_vertical_grid   = $config_vertical_grid
   config_met_interp      = $config_met_interp
   config_input_sst       = $config_input_sst
   config_frac_seaice     = $config_frac_seaice
/

&io
   config_pio_num_iotasks    = $config_pio_num_iotasks
   config_pio_stride         = $config_pio_stride
/

&decomposition
   config_block_decomp_file_prefix = '${graph_info_prefx}'
/
EOF

exit 0

#--------------------------------------------
# MPAS namelist.input
#--------------------------------------------

MPAS:

rm -f ./namelist.atmosphere
cat > ./namelist.atmosphere << EOF2
&nhyd_model
   config_time_integration_order = 2
   config_dt = $time_step
   config_run_duration = '${FCST_RANGE}:00:00'
   config_start_time   = '${START_DATE_MPAS}'
   config_stop_time    = '${END_DATE_MPAS}'
   config_split_dynamics_transport = .true.
   config_number_of_sub_steps = 4
   config_dynamics_split_steps = 3
   config_h_mom_eddy_visc2    = 0.0
   config_h_mom_eddy_visc4    = 0.0
   config_v_mom_eddy_visc2    = 0.0
   config_h_theta_eddy_visc2  = 0.0
   config_h_theta_eddy_visc4  = 0.0
   config_v_theta_eddy_visc2  = 0.0
   config_horiz_mixing        = '2d_smagorinsky'
   config_len_disp            = ${config_len_disp}
   config_visc4_2dsmag        = 0.05
   config_w_adv_order         = 3
   config_theta_adv_order     = 3
   config_scalar_adv_order    = 3
   config_u_vadv_order        = 3
   config_w_vadv_order        = 3
   config_theta_vadv_order    = 3
   config_scalar_vadv_order   = 3
   config_scalar_advection    = .true.
   config_positive_definite   = .false.
   config_monotonic           = .true.
   config_coef_3rd_order      = 0.25
   config_epssm               = 0.1
   config_smdiv               = 0.1
/

&damping
   config_mpas_cam_coef            = 2.0
   config_rayleigh_damp_u          = true
   config_zd                       = 16000.0
   config_xnutr                    = 0.2
   config_number_cam_damping_levels         = 8
   !config_nlevels_cam_damp         = 8
/

&io
   config_pio_num_iotasks    = $config_pio_num_iotasks ! use this for 15-/3-km grid
   config_pio_stride         = $config_pio_stride
/

&decomposition
   config_block_decomp_file_prefix = '${graph_info_prefx}'
/

&restart
   config_do_restart = .false.,    ! False for cold starts
/

&printout
    config_print_global_minmax_vel  = true
    config_print_detailed_minmax_vel = true
    config_print_global_minmax_sca  = true
/

&limited_area
   config_apply_lbcs = $config_apply_lbcs
/

&IAU
    config_IAU_option = 'off'
    config_IAU_window_length_s = 21600.
/

&physics
   config_sst_update          = ${local_update_sst}
   config_sstdiurn_update     = .false.
   config_deepsoiltemp_update = .false.
   config_radtlw_interval     = '00:${radiation_frequency}:00'
   config_radtsw_interval     = '00:${radiation_frequency}:00'
   config_bucket_update       = 'none' !'1_00:00:00'
   config_microp_re           = .true.
   config_lsm_scheme                = 'sf_ruc'
   num_soil_layers                  = 9
   config_physics_suite       = '${physics_suite}'
   config_convection_scheme   = '${deep_conv_param}'
/
EOF2

#&soundings
#    config_sounding_interval = '1:00:00'
#    config_sounding_output_path = './soundings/'
#/

#------------------------------------------

exit 0
