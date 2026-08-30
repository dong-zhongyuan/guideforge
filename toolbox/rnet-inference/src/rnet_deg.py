import sys
from os import path
import pandas as pd
import torch
from torch import nn
from utils import FinetunableRibonanzaNet, RNA_Dataset, load_config_from_yaml

USE_GPU = torch.cuda.is_available()

class finetuned_RibonanzaNet(FinetunableRibonanzaNet):
    def __init__(self, config):
        super(finetuned_RibonanzaNet, self).__init__(config)
        self.decoder = nn.Linear(config.ninp,5)
        
    def forward(self,src):
        sequence_features, pairwise_features=super(finetuned_RibonanzaNet, self).forward(src, torch.ones_like(src).long().to(src.device))
        output=self.decoder(sequence_features) #predict
        return output

model = finetuned_RibonanzaNet(load_config_from_yaml(path.join(path.dirname(__file__), '../RibonanzaNet', 'configs/pairwise.yaml')))
if USE_GPU:
    model = model.cuda()
model.load_state_dict(torch.load(path.join(path.dirname(__file__), '../RibonanzaNet-Weights', 'RibonanzaNet-Deg.pt'), map_location='cpu'))
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
            pred = model(sequence)
            if USE_GPU:
                pred = pred.cpu()
            test_preds.append(pred.numpy())

    reactivity = test_preds[0][0,:,0]
    deg_Mg_pH10 = test_preds[0][0,:,1]
    deg_pH10 = test_preds[0][0,:,2]
    deg_Mg_50C = test_preds[0][0,:,3]
    deg_50C = test_preds[0][0,:,4]

    print('reactivity:' + ','.join(str(x) for x in reactivity.tolist()))
    print('deg_Mg_pH10:' +','.join(str(x) for x in deg_Mg_pH10.tolist()))
    print('deg_pH10:' + ','.join(str(x) for x in deg_pH10.tolist()))
    print('deg_Mg_50C:' +','.join(str(x) for x in deg_Mg_50C.tolist()))
    print('deg_50C:' +','.join(str(x) for x in deg_50C.tolist()))