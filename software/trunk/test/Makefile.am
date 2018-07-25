BUILT_SOURCES   = .license.stamp

AUTOMAKE_OPTIONS = foreign
SUBDIRS = .

check_PROGRAMS  = # Append below
check_PROGRAMS += input_read_vel_grid
check_PROGRAMS += df_bkw
check_PROGRAMS += df_maxwell
check_PROGRAMS += df_rot_boltz
check_PROGRAMS += df_vib_boltz
check_PROGRAMS += coll_com_vel
check_PROGRAMS += coll_num_coll_pairs
check_PROGRAMS += coll_post_coll_vel
check_PROGRAMS += coll_rand_partner
check_PROGRAMS += coll_rand_vel_on_sphere
check_PROGRAMS += coll_rel_vel
check_PROGRAMS += ecoll_efrac_array
check_PROGRAMS += ecoll_g_prime
check_PROGRAMS += ecoll_pick_level
check_PROGRAMS += ecoll_rot_level
check_PROGRAMS += ecoll_vib_level
check_PROGRAMS += ecoll_rot_vib_level
check_PROGRAMS += ecoll_split_depl
check_PROGRAMS += interp_x
check_PROGRAMS += interp_y
check_PROGRAMS += interp_z
check_PROGRAMS += interp_xy
check_PROGRAMS += interp_xz
check_PROGRAMS += interp_yz
check_PROGRAMS += interp_xyz
check_PROGRAMS += props_dens
check_PROGRAMS += props_x_vel
check_PROGRAMS += props_y_vel
check_PROGRAMS += props_z_vel
check_PROGRAMS += props_tr_temp
check_PROGRAMS += props_rot_temp
check_PROGRAMS += props_vib_temp
check_PROGRAMS += props_temp
check_PROGRAMS += props_tr_energy
check_PROGRAMS += props_rot_energy
check_PROGRAMS += props_vib_energy
check_PROGRAMS += props_energy
check_PROGRAMS += props_rot_dof
check_PROGRAMS += props_vib_dof
check_PROGRAMS += props_Z_rot
check_PROGRAMS += props_Z_vib
check_PROGRAMS += bkw_relax
check_PROGRAMS += maxwell_relax
check_PROGRAMS += in_out

AM_CPPFLAGS = -I$(top_srcdir)/src
AM_FCFLAGS  = -I$(top_srcdir)/src $(GRVY_FCFLAGS)
AM_LDFLAGS  = -L$(top_builddir)/src -ldvm $(GRVY_FLIBS)

AM_CPPFLAGS += $(GSL_CFLAGS)
AM_LDFLAGS += $(GSL_LIBS)

input_read_vel_grid_SOURCES     = input_read_vel_grid_test.f90
df_bkw_SOURCES                  = df_bkw_test.f90
df_maxwell_SOURCES              = df_maxwell_test.f90
df_rot_boltz_SOURCES            = df_rot_boltz_test.f90
df_vib_boltz_SOURCES            = df_vib_boltz_test.f90
coll_com_vel_SOURCES            = coll_com_vel_test.f90
coll_num_coll_pairs_SOURCES     = coll_num_coll_pairs_test.f90
coll_post_coll_vel_SOURCES      = coll_post_coll_vel_test.f90
coll_rand_partner_SOURCES       = coll_rand_partner_test.f90
coll_rand_vel_on_sphere_SOURCES = coll_rand_vel_on_sphere_test.f90
coll_rel_vel_SOURCES            = coll_rel_vel_test.f90
ecoll_efrac_array_SOURCES       = ecoll_efrac_array_test.f90
ecoll_g_prime_SOURCES           = ecoll_g_prime_test.f90
ecoll_pick_level_SOURCES        = ecoll_pick_level_test.f90
ecoll_rot_level_SOURCES         = ecoll_rot_level_test.f90
ecoll_vib_level_SOURCES         = ecoll_vib_level_test.f90
ecoll_rot_vib_level_SOURCES     = ecoll_rot_vib_level_test.f90
ecoll_split_depl_SOURCES        = ecoll_split_depl_test.f90
interp_x_SOURCES                = interp_x_test.f90
interp_y_SOURCES                = interp_y_test.f90
interp_z_SOURCES                = interp_z_test.f90
interp_xy_SOURCES               = interp_xy_test.f90
interp_xz_SOURCES               = interp_xz_test.f90
interp_yz_SOURCES               = interp_yz_test.f90
interp_xyz_SOURCES              = interp_xyz_test.f90
props_dens_SOURCES              = props_dens_test.f90
props_x_vel_SOURCES             = props_x_vel_test.f90
props_y_vel_SOURCES             = props_y_vel_test.f90
props_z_vel_SOURCES             = props_z_vel_test.f90
props_tr_temp_SOURCES           = props_tr_temp_test.f90
props_rot_temp_SOURCES          = props_rot_temp_test.f90
props_vib_temp_SOURCES          = props_vib_temp_test.f90
props_temp_SOURCES              = props_temp_test.f90
props_tr_energy_SOURCES         = props_tr_energy_test.f90
props_rot_energy_SOURCES        = props_rot_energy_test.f90
props_vib_energy_SOURCES        = props_vib_energy_test.f90
props_energy_SOURCES            = props_energy_test.f90
props_rot_dof_SOURCES           = props_rot_dof_test.f90
props_vib_dof_SOURCES           = props_vib_dof_test.f90
props_Z_rot_SOURCES             = props_Z_rot_test.f90
props_Z_vib_SOURCES             = props_Z_vib_test.f90
bkw_relax_SOURCES               = bkw_relax_test.f90
maxwell_relax_SOURCES           = maxwell_relax_test.f90
in_out_SOURCES                  = in_out_test.f90

dvm_test_SRC = $(input_read_vel_grid_SOURCES) \
               $(df_bkw_SOURCES) \
               $(df_maxwell_SOURCES) \
               $(df_rot_boltz_SOURCES) \
               $(df_vib_boltz_SOURCES) \
               $(coll_com_vel_SOURCES) \
               $(coll_num_coll_pairs_SOURCES) \
               $(coll_post_coll_vel_SOURCES) \
               $(coll_rand_partner_SOURCES) \
               $(coll_rand_vel_on_sphere_SOURCES) \
               $(coll_rel_vel_SOURCES) \
               $(ecoll_efrac_array_SOURCES) \
               $(ecoll_g_prime_SOURCES) \
               $(ecoll_pick_level_SOURCES) \
               $(ecoll_rot_level_SOURCES) \
               $(ecoll_vib_level_SOURCES) \
               $(ecoll_rot_vib_level_SOURCES) \
               $(ecoll_split_depl_SOURCES) \
               $(interp_x_SOURCES) \
               $(interp_y_SOURCES) \
               $(interp_z_SOURCES) \
               $(interp_xy_SOURCES) \
               $(interp_xz_SOURCES) \
               $(interp_yz_SOURCES) \
               $(interp_xyz_SOURCES) \
               $(props_dens_SOURCES) \
               $(props_x_vel_SOURCES) \
               $(props_y_vel_SOURCES) \
               $(props_z_vel_SOURCES) \
               $(props_tr_temp_SOURCES) \
               $(props_rot_temp_SOURCES) \
               $(props_vib_temp_SOURCES) \
               $(props_temp_SOURCES) \
               $(props_tr_energy_SOURCES) \
               $(props_rot_energy_SOURCES) \
               $(props_vib_energy_SOURCES) \
               $(props_energy_SOURCES) \
               $(props_rot_dof_SOURCES) \
               $(props_vib_dof_SOURCES) \
               $(props_Z_rot_SOURCES) \
               $(props_Z_vib_SOURCES) \
               $(bkw_relax_SOURCES) \
               $(maxwell_relax_SOURCES) \
               $(in_out_SOURCES)


TESTS = input_prompt.sh \
        input_read_vel_grid \
        dist_func_prompt.sh \
        df_bkw \
        df_maxwell \
        df_rot_boltz \
        df_vib_boltz \
        collision_prompt.sh \
        coll_com_vel \
        coll_num_coll_pairs \
        coll_post_coll_vel \
        coll_rand_partner \
        coll_rand_vel_on_sphere \
        coll_rel_vel \
        energy_collision_prompt.sh \
        ecoll_efrac_array \
        ecoll_g_prime \
        ecoll_pick_level \
        ecoll_rot_level \
        ecoll_vib_level \
        ecoll_rot_vib_level \
        ecoll_split_depl \
        interpolation_prompt.sh \
        interp_x \
        interp_y \
        interp_z \
        interp_xy \
        interp_xz \
        interp_yz \
        interp_xyz \
        properties_prompt.sh \
        props_dens \
        props_x_vel \
        props_y_vel \
        props_z_vel \
        props_tr_temp \
        props_rot_temp \
        props_vib_temp \
        props_temp \
        props_tr_energy \
        props_rot_energy \
        props_vib_energy \
        props_energy \
        props_rot_dof \
        props_vib_dof \
        props_Z_rot \
        props_Z_vib \
        bkw_relax.sh \
        maxwell_relax.sh \
        in_out.sh \
        finalize.sh

CLEANFILES = *.mod

# Embedded license header support
#--------------------------------- 
.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE $(dvm_test_SRC)
	echo 'updated source license headers' >$@