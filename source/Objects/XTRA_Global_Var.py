import json

with open('./Directories.json') as fp:
    _param = json.load(fp)

FDIR_HALOS = _param["Halo_dir"]
FDIR_CLUSTERS = _param["Cluster_dir"]
FDIR_SB_MAP_SAVE = _param["SB_Map_dir"]
FDIR_EVENT_MAP = _param["Event_Map_dir"]
del _param
