import argparse
import os
from os import path
import pandas as pd
import numpy as np
import torch
from torch import nn
from utils import FinetunableRibonanzaNet, RNA_Dataset, load_config_from_yaml

USE_GPU = torch.cuda.is_available()

class finetuned_RibonanzaNet(FinetunableRibonanzaNet):
    def __init__(self, config):
        config.dropout=0.3
        super(finetuned_RibonanzaNet, self).__init__(config)

        self.dropout=nn.Dropout(0.0)
        self.ct_predictor=nn.Linear(64,1)

    def forward(self,src):
        _, pairwise_features=super(finetuned_RibonanzaNet, self).forward(src, torch.ones_like(src).long().to(src.device))
        pairwise_features=pairwise_features+pairwise_features.permute(0,2,1,3) #symmetrize
        output=self.ct_predictor(self.dropout(pairwise_features)) #predict

        return output.squeeze(-1)

model = finetuned_RibonanzaNet(load_config_from_yaml(path.join(path.dirname(__file__), '../RibonanzaNet', 'configs/pairwise.yaml')))
if USE_GPU:
    model = model.cuda()
model.load_state_dict(torch.load(path.join(path.dirname(__file__), '../RibonanzaNet-Weights', 'RibonanzaNet-SS.pt'), map_location='cpu'))
model.eval()



def dedupe_lists(list_of_lists):
    # Step 1: Convert each sublist to a sorted tuple
    tuple_set = {tuple(sorted(sublist)) for sublist in list_of_lists}
    
    # Step 2: Convert the set of tuples back to a list of lists
    deduped_list = [list(tup) for tup in tuple_set]
    
    return deduped_list

def detect_crossed_pairs(bp_list):
    """
    Detect crossed base pairs in a list of base pairs in RNA secondary structure.

    Args:
    bp_list (list of tuples): List of base pairs, where each tuple (i, j) represents a base pair.
    
    Returns:
    list of tuples: List of crossed base pairs.
    """
    crossed_pairs = []
    # Iterate through each pair of base pairs
    for i in range(len(bp_list)):
        for j in range(i+1, len(bp_list)):
            bp1 = bp_list[i]
            bp2 = bp_list[j]

            # Check if they are crossed
            if (bp1[0] < bp2[0] < bp1[1] < bp2[1]) or (bp2[0] < bp1[0] < bp2[1] < bp1[1]):
                crossed_pairs.append(bp1)
                crossed_pairs.append(bp2)
    return dedupe_lists(crossed_pairs)

if __name__ == '__main__':
    from tqdm import tqdm

    parser = argparse.ArgumentParser()
    parser.add_argument(dest='sequence', help='Sequence to predict')
    parser.add_argument('--batch-size', type=int, dest='batch_size', help='batch size (number of predictions to run simultaniously - requires more memory, but allows for increased parallelization)', default=1)
    parser.add_argument('--output-confidence', action='store_true', dest='output_pair_confidence', help='inclue pair confidence matrix in output')
    args = parser.parse_args()
    seq = args.sequence
    test_dataset=RNA_Dataset(pd.DataFrame([{'sequence': seq}]))

    test_preds=[]
    for i in tqdm(range(len(test_dataset))):
        example=test_dataset[i]
        sequence=example['sequence']
        if USE_GPU:
            sequence = sequence.cuda()
        sequence = sequence.unsqueeze(0)

        with torch.no_grad():
            pred = model(sequence).sigmoid()
            if USE_GPU:
                pred = pred.cpu()
            test_preds.append(pred.numpy())

    # create dummy arnie config
    os.environ['NUPACKHOME'] = '/tmp'
    from arnie.pk_predictors import _hungarian
    from arnie.utils import get_expected_accuracy, convert_bp_list_to_dotbracket

    def mask_diagonal(matrix, mask_value=0):
        matrix=matrix.copy()
        n = len(matrix)
        for i in range(n):
            for j in range(n):
                if abs(i - j) < 4:
                    matrix[i][j] = mask_value
        return matrix

    hungarian_structures=[]
    confidences=[]
    cp_masked_confidences=[]
    for i in range(len(test_preds)):
        s,bp=_hungarian(mask_diagonal(test_preds[i][0]), theta=0.5,min_len_helix=1) #best theta based on val is 0.5
        hungarian_structures.append(s)
        confidences.append(
            get_expected_accuracy(s, test_preds[i][0], mode='fscore') if len(bp) > 0 else 0
        )
        # With this formulation, we are effectively saying any nucleotides not predicted to be
        # part of a cross pair are "correct", as we specify they should be unpaired
        # and there are no pairing probabilities at that position > 0. This is specifically tuned
        # for F1, in which the "true negative" component is disgarded in its formula.
        # As put by Rhiju:
        # The motivation for this choice is that in cases with inferred pseudoknots,
        # we really don't care about what residues outside the relevant crossed-pairs are doing.
        # We just want an estimate of whether the specific crossed-pairs will show up in the actual
        # structure. There is an analogy to how we set up the 'crossed pair quality' component for OpenKnotScore.
        crossed_pairs = detect_crossed_pairs(bp)
        cross_pair_res = sorted([res for pair in crossed_pairs for res in pair])
        mask_matrix = np.zeros((len(s), len(s)))
        mask_matrix[:,cross_pair_res] = 1
        mask_matrix[cross_pair_res,:] = 1
        masked_pred = test_preds[i][0] * mask_matrix
        cp_masked_confidences.append(
            get_expected_accuracy(convert_bp_list_to_dotbracket(sorted(crossed_pairs), len(s)), masked_pred, mode='fscore') if len(crossed_pairs) > 0 else 0
        )

    print('structure:' + hungarian_structures[0])
    print('confidence_f1:' + str(confidences[0]))
    print('confidence_f1_cross_pair_masked:' + str(cp_masked_confidences[0]))
    if args.output_pair_confidence:
        print('pair_confidence:' + ','.join(str(x) for xs in test_preds[0][0].tolist() for x in xs))
