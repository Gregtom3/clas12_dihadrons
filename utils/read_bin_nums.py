import sys
import yaml

yaml_file = sys.argv[1]
scheme_num = int(sys.argv[2])
with open(yaml_file, 'r') as stream:
    try:
        data = yaml.safe_load(stream)
        for i,binning in enumerate(data['binningStructures']):
            if(i!=scheme_num):
                continue
            nbins=1
            for bin_edges in binning['binEdges']:
                nbins*=(len(bin_edges)-1)
            print(nbins)
    except yaml.YAMLError as exc:
        print(exc)