# wrapper for the main FortEPiaNO code
import numpy as np
import six
import sys
import prepareIni as pim
from fortepianoWrapper import fortepianowrapper as fpw


def get_version():
    """Print the version of the fortran code"""
    return fpw.get_version() if six.PY2 else str(fpw.get_version(), "utf-8")


def runConfig(args, values, logname="", verbose=0, no_interpolation=False):
    nnu = values["nnu"]
    dms = [0.0] + [values["dm%d1" % i] for i in range(2, nnu + 1)]
    mixang = np.zeros((nnu, nnu))
    for j in range(2, nnu + 1):
        for i in range(1, j):
            mixang[i - 1, j - 1] = values["th%d%d" % (i, j)]

    fpw.w_initial_operations(logname if logname != "" else (args.inifile + ".log"))
    fpw.pass_basic_params(
        nnu,  # flavorNumber
        True,  # force_replace
        args.outputfolder,  # outputFolder
        0,  # args.num_threads, # num_threads
    )
    fpw.w_check_output_folder()
    fpw.w_check_num_threads()

    fpw.pass_verbosity(
        args.verbose,  # verbose
        args.verbose_deriv_freq,  # Nprintderivs
    )

    fpw.pass_xy(
        args.x_in,  # x_in
        args.x_fin,  # x_fin
        args.Nx,  # Nx
        not args.no_GL,  # use_gauss_laguerre
        args.Ny,  # Ny
        args.y_min,  # y_min
        args.y_max,  # y_max
    )
    fpw.w_set_xy_vectors()

    fpw.pass_precision(
        100,  # maxiter
        1e-7,  # toler_jkyg
        args.dlsoda_atol_z,  # dlsoda_atol_z
        args.dlsoda_atol_d,  # dlsoda_atol_d
        args.dlsoda_atol_o,  # dlsoda_atol_o
        args.dlsoda_rtol,  # dlsoda_rtol
    )
    fpw.w_set_interp()

    fpw.pass_output_config(
        True,  # checkpoint
        args.save_BBN,  # save_BBN
        args.save_fd,  # save_fd
        args.save_energy_entropy,  # save_energy_entropy_evolution
        args.save_Neff,  # save_Neff
        args.save_nuDens,  # save_nuDens_evolution
        args.save_number,  # save_number_evolution
        args.save_w,  # save_w_evolution
        args.save_z,  # save_z_evolution
        args.save_intermediate,  # intermediateSteps%output
    )

    fpw.w_allocate_stuff()

    fpw.pass_collint(
        values["collint_damping_type"],  # collint_damping_type
        values["collint_diagonal_zero"] == "T",  # collint_diagonal_zero
        values["collint_offdiag_damping"] == "T",  # collint_offdiag_damping
        args.collint_d_no_nue,  # collint_d_no_nue
        args.collint_d_no_nunu,  # collint_d_no_nunu
        args.collint_od_no_nue,  # collint_od_no_nue
        args.collint_od_no_nunu,  # collint_od_no_nunu
        [
            [False for j in range(nnu)] for i in range(nnu)
        ],  # damping_read_zero (nf*nf matrix!)
    )
    fpw.pass_ftqed(
        values["ftqed_temperature_corr"] == "T",  # ftqed_temperature_corr
        values["ftqed_log_term"] == "T",  # ftqed_log_term
        values["ftqed_ord3"] == "T",  # ftqed_ord3
        args.ftqed_e_mth_ld,  # ftqed_e_mth_leptondens
    )

    fpw.pass_nu_def(
        values["factors_v"],  # nuFactor (nf vector)
        values["sterile_v"],  # sterile (nf vector)
    )
    fpw.w_set_nu_properties()

    fpw.w_init_fermions()
    fpw.w_zin_solver()

    fpw.pass_nu_mixing(
        args.use_sinsq,  # giveSinSq
        dms,  # massSplittings (nf vector, first element ignored)
        mixang,  # mixingAngles (nf*nf matrix)
    )

    fpw.w_set_mixing_matrix()
    fpw.w_set_mass_matrix()
    fpw.w_init_matrices()
    fpw.w_end_config()


if __name__ == "__main__":
    parser = pim.setParser(prog="fortepiano.py")
    args = parser.parse_args(sys.argv[1:])
    values = pim.getIniValues(args)
    runConfig(args, values)
    # fpw.w_solver()
    # print("Wrapper for FortEPiaNO v%s" % getVersion())
