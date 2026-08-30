import argparse
import pandas as pd
import torch
from torch.utils.data import DataLoader
from rnet_shape import model
from utils import RNA_Dataset

USE_GPU = torch.cuda.is_available()

if __name__ == '__main__':
    from tqdm import tqdm

    parser = argparse.ArgumentParser()
    parser.add_argument(dest='sequence', help='Sequence to predict')
    parser.add_argument('--batch-size', type=int, dest='batch_size', help='batch size (number of predictions to run simultaniously - requires more memory, but allows for increased parallelization)', default=1)
    args = parser.parse_args()
    seq = args.sequence
    batch_size = args.batch_size
    
    rc = {'A':'U','U':'A','C':'G','G':'C'}
    max_len = len(seq)
    m2_sequences = []
    for i in range(max_len):
        m2_sequences.append(seq[:i] + rc[seq[i]] + seq[i+1:])

    data = [ [x] for x in m2_sequences ]
    test_data = pd.DataFrame(data,columns = ['sequence'])
    test_dataset=RNA_Dataset(test_data)
    dataloader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

    test_preds = []
    with torch.no_grad():
        for batch in tqdm(dataloader):
            inputs = batch['sequence']
            if USE_GPU:
                inputs = inputs.cuda()
            outputs = model(inputs,torch.ones_like(inputs))
            test_preds.append(outputs)
    output = torch.cat(test_preds, dim=0)
    if USE_GPU:
        output = output.cpu()
    output = output.numpy()

    print('2a3:' + ','.join(str(x) for xs in output[:,:,0].tolist() for x in xs))
    print('dms:' +','.join(str(x) for xs in output[:,:,1].tolist() for x in xs))
