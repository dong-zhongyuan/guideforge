import sys
from os import path
import pandas as pd
import torch
sys.path.append(path.join(path.dirname(__file__), '../RibonanzaNet'))
from Network import RibonanzaNet
from utils import RNA_Dataset, load_config_from_yaml

USE_GPU = torch.cuda.is_available()

model = RibonanzaNet(load_config_from_yaml(path.join(path.dirname(__file__), '../RibonanzaNet', 'configs/pairwise.yaml')))
if USE_GPU:
    model = model.cuda()
model.load_state_dict(torch.load(path.join(path.dirname(__file__), '../RibonanzaNet-Weights', 'RibonanzaNet.pt'), map_location='cpu'))
model.eval()

if __name__ == '__main__':
    from tqdm import tqdm

    seq = sys.argv[1]
    test_dataset=RNA_Dataset(pd.DataFrame([{'sequence': seq}]))

    test_preds=[]
    for i in tqdm(range(len(test_dataset))):
        example=test_dataset[i]
        sequence=example['sequence']
        if USE_GPU:
            sequence = sequence.cuda()
        sequence = sequence.unsqueeze(0)

        with torch.no_grad():
            pred = model(sequence,torch.ones_like(sequence)).squeeze()
            if USE_GPU:
                pred = pred.cpu()
            test_preds.append(pred.numpy())

    pred_2a3 = test_preds[0][:,0]
    pred_dms = test_preds[0][:,1]

    print('2a3:' + ','.join(str(x) for x in pred_2a3.tolist()))
    print('dms:' +','.join(str(x) for x in pred_dms.tolist()))