
import yaml
import argparse
import main

def begin_experiment(params):  
    main.run(params) 

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-config", help = "2-D Waveguide Adaptation")
    args = parser.parse_args()

    if(args.config == None):
        print("\nAttach Configuration File! Run experiment.py -h\n")
        exit()

    params = yaml.load(open(args.config), Loader = yaml.FullLoader)

    print (params['cell_config']['resolution'])
    #input()
    begin_experiment(params)