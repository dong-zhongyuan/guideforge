import sys
from os import path
import numpy as np
import torch
from torch.utils.data import Dataset
import yaml

sys.path.append(path.join(path.dirname(__file__), '../RibonanzaNet'))
from Network import RibonanzaNet

class RNA_Dataset(Dataset):
    def __init__(self,data):
        self.data=data
        self.tokens={nt:i for i,nt in enumerate('ACGU')}

    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        sequence=[self.tokens[nt] for nt in self.data.loc[idx,'sequence']]
        sequence=np.array(sequence)
        sequence=torch.tensor(sequence)

        return {'sequence':sequence}

class Config:
    def __init__(self, **entries):
        self.__dict__.update(entries)
        self.entries=entries

    def print(self):
        print(self.entries)

def load_config_from_yaml(file_path):
    with open(file_path, 'r') as file:
        config = yaml.safe_load(file)
    return Config(**config)

class FinetunableRibonanzaNet(RibonanzaNet):
    def forward(self, src, src_mask=None, return_aw=False):
        B,L=src.shape

        src = src.long()

        src = self.encoder(src).reshape(B,L,-1)
        pairwise_features=self.outer_product_mean(src)
        pairwise_features=pairwise_features+self.pos_encoder(src)
        for layer in self.transformer_encoder:
            src,pairwise_features=layer(src, pairwise_features, src_mask,return_aw=return_aw)
        return src, pairwise_features
