import sys
from os import path
import pandas as pd
import torch
from torch import nn
from utils import FinetunableRibonanzaNet, RNA_Dataset, load_config_from_yaml, Config

USE_GPU = torch.cuda.is_available()

class finetuned_RibonanzaNet(nn.Module):
    """
    Add a final dense layer to RibonanzaNet in order to predict the number of reads 
    during chemical mapping experiments. 
    
    Parameters 
    ----------
    out_size : int
        The number of output features, which corresponds to the number of different 
        chemical mapping experiments. In Ribonanza, this is equal to 2, corresponding 
        to the 2A3 and DMS experiments. 
    """
    def __init__(self, out_size :int, config: Config) -> None:
        super().__init__()
        self.model = FinetunableRibonanzaNet(config)
        self.head = nn.Linear(256, out_size)
        
        
    def forward(self, x : torch.Tensor) -> torch.Tensor:
        """
        A forward pass of the model.
            
        Parameters 
        ----------
        x : torch.Tensor
            The input sequences, embedded using the dictionary 
            {'A' : 0, 'C' : 1, 'G' : 2, 'U' : 3}.
            
        Returns
        -------
        torch.Tensor:
            The number of reads for each sequence in the input and each chemical 
            mapping experiment. 
        """
        # compute the RibonanzaNet embeddings
        embeddings, _ = self.model(x, torch.ones_like(x).long().to(x.device))
        # pass through the dense head and take the mean over the sequence dimension
        y = torch.mean(
            self.head(embeddings),
            dim = -2,
        )
        # pass through a nonlinearity to return positive values
        return torch.exp(y)

model = finetuned_RibonanzaNet(out_size = 2, config = load_config_from_yaml(path.join(path.dirname(__file__), '../RibonanzaNet', 'configs/pairwise.yaml')))
if USE_GPU:
    model = model.cuda()
state_dict = torch.load(path.join(path.dirname(__file__), '../RibonanzaNet-Weights', 'RibonanzaNet-Drop.pt'), map_location='cpu')
state_dict['head.weight'] = state_dict.pop('head.layers.0.weight')
state_dict['head.bias'] = state_dict.pop('head.layers.0.bias')
model.load_state_dict(state_dict)
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

    pred_2a3 = test_preds[0][0,0]
    pred_dms = test_preds[0][0,1]

    print(f'2a3:{pred_2a3}')
    print(f'dms:{pred_dms}')
