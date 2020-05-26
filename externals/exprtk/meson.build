check_dir = run_command(python, 
                        script_path / 'check_dir.py',
                        'exprtk.hpp')
if check_dir.returncode() != 0
  message('Extracting exprtk ...')
  extract_boost = run_command(python, 
                              script_path / 'extract_file.py',
                              'exprtk.tar.gz', 
                              meson.current_source_dir(), check: true)
else
  message('Exprtk sources found.')
endif

su2_deps += declare_dependency(include_directories: './')
