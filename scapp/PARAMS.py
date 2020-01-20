# advanced parameters for the Recycler2 workflow
# create global variables set at beginning of pipeline and accessible from everywhere

import json
import os


#
MAX_CV = 0.5
MIN_LENGTH = 1000
#

#
CLASSIFICATION_THRESH = 0.5 # threshold for classifying potential plasmid
GENE_MATCH_THRESH = 0.75 # threshold for % identity and fraction of length to match plasmid genes
SELF_LOOP_SCORE_THRESH = 0.9 # threshold for self-loop plasmid score
SELF_LOOP_MATE_THRESH = 0.1 # threshold for self-loop off loop mates
CHROMOSOME_SCORE_THRESH = 0.2 # threshold for high confidence chromosome node score
CHROMOSOME_LEN_THRESH = 10000 # threshold for high confidence chromosome node length
PLASMID_SCORE_THRESH = 0.9 # threshold for high confidence plasmid node score
PLASMID_LEN_THRESH = 10000 # threshold for high confidence plasmid node length
GOOD_CYC_DOMINATED_THRESH = 0.5 # threshold for # of mate-pairs off the cycle in dominated node

def load_params_json():

    global CLASSIFICATION_THRESH
    global GENE_MATCH_THRESH
    global SELF_LOOP_SCORE_THRESH
    global SELF_LOOP_MATE_THRESH
    global CHROMOSOME_SCORE_THRESH
    global CHROMOSOME_LEN_THRESH
    global PLASMID_SCORE_THRESH
    global PLASMID_LEN_THRESH
    global GOOD_CYC_DOMINATED_THRESH


    ''' Load the parameters from the params.json file
    '''
    curr_file_path = os.path.dirname(os.path.abspath(__file__))
    try:
        json_path = os.path.join(curr_file_path,'params.json')
        with open(json_path) as j:
            params = json.load(j)
            MAX_CV = params['MAX_CV']
            MIN_LENGTH = params['MIN_LENGTH']
            CLASSIFICATION_THRESH = params['CLASSIFICATION_THRESH']
            GENE_MATCH_THRESH = params['GENE_MATCH_THRESH']
            SELF_LOOP_SCORE_THRESH = params['SELF_LOOP_SCORE_THRESH']
            SELF_LOOP_MATE_THRESH = params['SELF_LOOP_MATE_THRESH']
            CHROMOSOME_SCORE_THRESH = params['CHROMOSOME_SCORE_THRESH']
            CHROMOSOME_LEN_THRESH = params['CHROMOSOME_LEN_THRESH']
            PLASMID_SCORE_THRESH = params['PLASMID_SCORE_THRESH']
            PLASMID_LEN_THRESH = params['PLASMID_LEN_THRESH']
            GOOD_CYC_DOMINATED_THRESH = params['GOOD_CYC_DOMINATED_THRESH']
    except:
        print("Error loading parameters. Please check params.json file")
        raise
