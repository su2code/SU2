check_dir = run_command('test', '-d', 'boost')
if check_dir.returncode() != 0
  message('Extracting boost ...')
  extract_boost = run_command('tar', '-xzf', 'boost.tar.gz', '-C', './')
else
  message('Boost sources found.')
endif

