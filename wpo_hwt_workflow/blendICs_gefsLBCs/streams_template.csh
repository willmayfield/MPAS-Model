#!/bin/csh
#
#-----------------------------------------------------------------------
# Script streams_template.csh
#
# Purpose: create streams templates for MPAS
#
#-----------------------------------------------------------------------

# Need to specify defaults for variables defined in scripts other than the driver
if ( ! $?update_sst_interval )     set update_sst_interval    = "none"
if ( ! $?ic_output_interval  )     set ic_output_interval     = "initial_only"
if ( ! $?lbc_output_interval )     set lbc_output_interval    = "none"

#####
if ( $MPAS_REGIONAL == .true. || $MPAS_REGIONAL == true ) then
   set lbc_input_interval = ${LBC_FREQ}:00:00
else
   set lbc_input_interval = none
endif

echo "Filling streams for $1"

if ( $1 == "mpas" || $1 == "MPAS" ) then
   goto MPAS
else if ( $1 == "mpas_init" || $1 == "MPAS_init" ) then
   goto MPAS_init
echo "wrong usage"
   echo "use as $0 (mpas,mpas_init)"
   exit
endif

#--------------------------------------------
# MPAS_initialization namelist.input
#--------------------------------------------

MPAS_init:

rm -f ./streams.init_atmosphere
cat > ./streams.init_atmosphere << EOF
<streams>
<immutable_stream name="input"
                  type="input"
                  filename_template="${config_input_name}"
                  input_interval="initial_only" />

<immutable_stream name="output"
        type="output"
        io_type="pnetcdf,cdf5"
        filename_template="${config_output_name}"
        clobber_mode="replace_files"
	precision="single"
        packages="initial_conds"
        output_interval="${ic_output_interval}" />

<immutable_stream name="surface"
                  type="output"
		  io_type="pnetcdf,cdf5"
                  filename_template="./sfc_update.nc"
                  filename_interval="none"
                  packages="sfc_update"
                  clobber_mode="replace_files"
	          precision="single"
                  output_interval="24:00:00"/>

<immutable_stream name="lbc"
                  type="output"
		  io_type="pnetcdf,cdf5"
                  filename_template="lbc.\$Y-\$M-\$D_\$h.\$m.\$s.nc"
                  filename_interval="output_interval"
                  packages="lbcs"
		  clobber_mode="replace_files"
                  output_interval="${lbc_output_interval}" />

</streams>
EOF

exit 0

#--------------------------------------------
# MPAS namelist.input
#--------------------------------------------

MPAS:

rm -f ./streams.atmosphere
cat > ./streams.atmosphere << EOF2
<streams>
<immutable_stream name="input"
                  type="input"
                  filename_template="./init.nc"
		  io_type="pnetcdf,cdf5"
                  input_interval="initial_only"/>

<immutable_stream name="restart"
        type="output"
        filename_template="restart.\$Y-\$M-\$D_\$h.\$m.\$s.nc"
        input_interval="initial_only"
        io_type="pnetcdf,cdf5"
        clobber_mode="replace_files"
	precision="single"
        output_interval="${restart_output_interval}:00:00" />

<stream name="output"
        type="output"
        filename_template="./history.\$Y-\$M-\$D_\$h.\$m.\$s.nc"
        clobber_mode="replace_files"
	precision="single"
        io_type="netcdf4"
        output_interval="${diag_output_interval}:00:00">

    <file name="${SCRIPT_DIR}/stream_list.atmosphere.output"/>
</stream>

<stream name="diagnostics"
        type="output"
        filename_template="diag.\$Y-\$M-\$D_\$h.\$m.\$s.nc"
        clobber_mode="replace_files"
	precision="single"
	io_type="netcdf4"
        output_interval="${diag_output_interval}:00:00">

    <file name="${SCRIPT_DIR}/stream_list.atmosphere.diagnostics"/>
</stream>

<stream name="surface"
        type="input"
        io_type="pnetcdf,cdf5"
        filename_interval="none"
        filename_template="./sfc_update.nc"
        input_interval="$update_sst_interval">
    <file name="${SCRIPT_DIR}/stream_list.atmosphere.surface"/>
</stream>

<immutable_stream name="lbc_in"
                 type="input"
                 io_type="pnetcdf,cdf5"
                 filename_template="lbc.\$Y-\$M-\$D_\$h.\$m.\$s.nc"
                 filename_interval="input_interval"
                 input_interval="${lbc_input_interval}" />

</streams>
EOF2

##########

exit 0
