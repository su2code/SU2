install_data(['continuous_adjoint.py',
             'compute_uncertainty.py',
	     'finite_differences.py',
	     'mesh_deformation.py',
	     'parallel_computation.py',
	     'parallel_computation_fsi.py',
	     'shape_optimization.py',
	     'merge_solution.py',
	     'set_ffd_design_var.py',
	     'compute_polar.py',
	     'discrete_adjoint.py',
	     'direct_differentiation.py',
	     'fsi_computation.py',
	     'SU2_CFD.py'],
	     install_dir: join_paths(get_option('bindir')))

install_data(['SU2/__init__.py'], install_dir: join_paths(get_option('bindir'), 'SU2/'))

install_data(['SU2/eval/design.py',
    	      'SU2/eval/functions.py',
    	      'SU2/eval/gradients.py',
    	      'SU2/eval/__init__.py'],
	      install_dir: join_paths(get_option('bindir'), 'SU2/eval'))

install_data(['SU2/io/config.py',
              'SU2/io/config_options.py',
              'SU2/io/data.py',
              'SU2/io/filelock.py',
              'SU2/io/redirect.py',
              'SU2/io/state.py',
              'SU2/io/tools.py',
              'SU2/io/historyMap.py',
              'SU2/io/__init__.py'],
	      install_dir: join_paths(get_option('bindir'), 'SU2/io'))

install_data(['SU2/opt/project.py',
              'SU2/opt/scipy_tools.py',
              'SU2/opt/__init__.py'],
	      install_dir: join_paths(get_option('bindir'), 'SU2/opt'))

install_data(['SU2/run/adjoint.py',
              'SU2/run/deform.py',
              'SU2/run/direct.py',
              'SU2/run/interface.py',
              'SU2/run/merge.py',
              'SU2/run/geometry.py',
              'SU2/run/projection.py',
              'SU2/run/__init__.py'],
	      install_dir: join_paths(get_option('bindir'), 'SU2/run'))

install_data(['SU2/util/bunch.py',
              'SU2/util/filter_adjoint.py',
              'SU2/util/lhc_unif.py',
              'SU2/util/misc.py',
              'SU2/util/mp_eval.py',
              'SU2/util/ordered_bunch.py',
              'SU2/util/ordered_dict.py',
              'SU2/util/plot.py',
              'SU2/util/polarSweepLib.py',
              'SU2/util/which.py',
              'SU2/util/switch.py',
              'SU2/util/__init__.py'],
	      install_dir: join_paths(get_option('bindir'), 'SU2/util'))

install_data(['FSI_tools/__init__.py',
              'FSI_tools/FSIInterface.py',
              'FSI_tools/FSI_config.py'],
	      install_dir: join_paths(get_option('bindir'), 'FSI_tools'))

install_data(['SU2_Nastran/__init__.py',
              'SU2_Nastran/pysu2_nastran.py'],
        install_dir: join_paths(get_option('bindir'), 'SU2_Nastran'))
