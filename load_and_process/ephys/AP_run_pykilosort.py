from pathlib import Path
from pykilosort import run, add_default_handler, np1_probe, np2_probe, neuropixel_probe_from_metafile

# Print pykilosort starting notice
print('Running pykilosort: ', data_filename)

# (NOT USED) Automatically detect probe type
# (I doubt this works with open ephys, haven't tried yet)
# ProbeType = neuropixel_probe_from_metafile(data_filename)

# Run pykilosort
# NOTE: currently hard-coded np1 probe WITHOUT sync channel
# (open ephys has no sync channel, spikeglx has sync channel)
add_default_handler(level='INFO') # print output as the algorithm runs
run(data_filename, probe=np1_probe(sync_channel=False), dir_path=pykilosort_output_path)

# Print pykilosort ending notice
print('Finished pykilosort: ', pykilosort_output_path)
