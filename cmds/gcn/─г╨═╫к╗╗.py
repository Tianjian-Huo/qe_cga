import torch
#from gcn import GCN, GCNF
from f8_10 import GCNF

"""net = GCN(36,28,80,120,72)
net = torch.load('energy.pth', map_location='cpu')
torch.save(net.state_dict(),'energy.pth')"""

netf = GCNF(36,28,80,120)
netf = torch.load('force.pth', map_location='cpu')
torch.save(netf.state_dict(),'force.pth')