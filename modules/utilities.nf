process COLLECTFILEPATHS {

  input:
  path filepath

	output:
  path "data_paths.txt", emit: data_paths

	script:
	"""
  printf "%s\n" ${filepath.join('\n')} > data_paths.txt
	"""
}

