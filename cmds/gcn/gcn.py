'''cd对角线元素改为0，改变self.beta*h_prime的位置，使得临近原子和中心原子分来聚合
注意力还乘adj，效果明显好
同gcn_torch4_60''''''cutoff=6,加2个bn层'''
import torch
from torch import nn
import torch.nn.functional as F
from math import pi


class Bond_Feather_Extract(nn.Module): 
   
    def __init__(self, nfd):
        super(Bond_Feather_Extract, self).__init__()
        #Parameter类是Tensor的子类，如果一个Tensor是Parameter，那么它会自动被添加到模型的参数列表里。
        self.fp = nn.Parameter(torch.Tensor(1,1,1,nfd))
        self.coeff = nn.Parameter(torch.Tensor(self.fp.size()))
        self.reset_parameters()

    def reset_parameters(self):
        self.fp.data.uniform_(1, 6)
        self.coeff.data.uniform_(2, 3)
        
    def forward(self, x, cd):
        '''输入b*N*N，输出b*N*nfd
        输入b*N*N*N，输出b*N*N*nfd'''
        xx = torch.unsqueeze(x, -1)
        fd = torch.exp(-100.*self.coeff*(xx-self.fp)*(xx-self.fp))
        N = x.shape[1]
        mask = (torch.ones(N,N,dtype=x.dtype)-torch.eye(N,dtype=x.dtype)).to(x.device)
        mask = torch.unsqueeze(torch.unsqueeze(mask,0), -1)
        fd = fd * mask * cd.unsqueeze(-1)
        return fd.sum(axis=-2) / 12
       
       
class Angle_Feather_Extract(nn.Module):

    def __init__(self, nfa=5):
        '''nfd：距离特征数
        nfa：角度特征数'''
        super(Angle_Feather_Extract, self).__init__()
        self.fp = nn.Parameter(torch.Tensor(1,1,1,1,nfa))
        self.coeff = nn.Parameter(torch.Tensor(self.fp.size()))
        self.reset_parameters()

    def reset_parameters(self):
        self.fp.data.uniform_(-1, 1)
        self.coeff.data.uniform_(2, 3)
        
    def calc_angle(self, d):
        '''余弦定理计算角度
        对角线元素为NaN'''
        d_2 = d*d
        d1 = torch.unsqueeze(d_2,-1)
        d2 = torch.unsqueeze(d_2,-2)
        d3 = torch.unsqueeze(d_2,-3)
        a = (d1+d2-d3) / (2*torch.unsqueeze(d,-1)*torch.unsqueeze(d,-2)) #cos
        a[torch.isnan(a)] = 0
        return a
    
    def forward(self, d, cd):
        '''d:邻接矩阵，b*N*N'''
        a = self.calc_angle(d) #计算角度
        a = a.unsqueeze(-1)
        temp = torch.exp(-100.*self.coeff*(a-self.fp)*(a-self.fp)) #角度特征
        N = d.shape[1]
        cd0 = cd.clone()
        cd0[cd0==1] = 0
        dij = cd0.view(-1,N,N,1,1)
        dik = cd0.view(-1,N,1,N,1)
        fa = (temp*dij*dik).sum(axis=(-2,-3)) /(N*12)
        #fa[torch.isnan(fa)] = 0
        return fa
        
        
class Feather_Extract(nn.Module):
    '''从距离矩阵中提取键长和键角特征'''

    def __init__(self, nfd=3, nfa=5):
        '''nfd：键长特征数
        nfa：键角特征数'''
        super(Feather_Extract, self).__init__()
        self.fsd = Bond_Feather_Extract(nfd)
        self.fsa = Angle_Feather_Extract(nfa)
        self.cutoff = 6
    
    def forward(self, d):
        '''d:邻接矩阵，b*N*N'''
        cd = d.clone()
        temp = cd > self.cutoff
        cd = (1+torch.cos(pi*d/self.cutoff))/2
        cd[temp] = 0
        cd[cd==1] = 0
        fb = self.fsd(d, cd) #键长特征
        fa = self.fsa(d, cd) #键角特征
        return cd, torch.cat([fb,fa], axis=-1) #合并特征
 
 
class GraphAttentionLayer(nn.Module):
    '''图卷积层'''
    def __init__(self, n_in, n_out, use_bias=True):
        super(GraphAttentionLayer, self).__init__()
        self.n_in = n_in
        self.n_out = n_out
        self.use_bias = use_bias
        #定义GCN层的权重矩阵
        self.weight = nn.Parameter(torch.Tensor(n_in, n_out))
        if self.use_bias:
            self.bias = nn.Parameter(torch.Tensor(n_out))
        else:
            self.register_parameter('bias', None)
        self.a = nn.Parameter(torch.Tensor(2*n_out+1))
        self.beta = nn.Parameter(torch.Tensor(1))
        self.reset_parameters() #使用自定义的参数初始化方式
 
    def reset_parameters(self):
        #自定义参数初始化方式
        #权重参数初始化方式
        nn.init.kaiming_uniform_(self.weight)
        self.a.data.uniform_(-0.01,0.01)
        self.beta.data.uniform_(-0.4, -0.1)
        if self.use_bias: #偏置参数初始化为0
            nn.init.zeros_(self.bias)
 
    def forward(self, adj, x):
        B, N = x.shape[:2]
        h = torch.matmul(x, self.weight) #B*N*n_out
        a_input = torch.cat([h.repeat(1,1,N).view(B,N*N,self.n_out),
                             h.repeat(1,N,1),
                             adj.view(B,N*N,1)], dim=2).view(B,N,N,2*self.n_out+1)
        attention = F.softmax(torch.nn.LeakyReLU()(torch.matmul(a_input, self.a)), dim=2) #B*N*N
        h_prime = torch.bmm(attention, h) #B*N*n_out
        if self.use_bias:
            out = self.beta*h_prime + h + self.bias
        else:
            out = self.beta*h_prime + h
        return F.elu(out) 
        
        

        
        
class ForceGraphAttentionLayer(nn.Module):
    '''力的图卷积层'''
    def __init__(self, n_in, use_bias=True):
        '''n_in：原子的特征数'''
        super(ForceGraphAttentionLayer, self).__init__()
        self.use_bias = use_bias
        self.a = nn.Parameter(torch.Tensor(2*n_in+1, 8))
        self.a.data.uniform_(-0.01,0.01)
        self.a2 = nn.Parameter(torch.Tensor(8))
        self.a2.data.uniform_(-0.01,0.01)
 
    def forward(self, adj, c, x):
        '''adj：距离矩阵,B*N*N
        c：坐标矩阵,B*N*3
        x：特征,B*N*n_in'''
        B, N = x.shape[:2]
        a_input = torch.cat([x.repeat(1,1,N).view(B,N*N,-1),
                             x.repeat(1,N,1),
                             adj.view(B,N*N,1)], dim=2).view(B,N,N,-1)
        a_input = a_input.view(B*N*N,-1)
        attention = torch.nn.LeakyReLU()(torch.mm(a_input, self.a)) #B*N*N*8
        attention = torch.matmul(attention, self.a2)
        attention = attention.view(B,N,N,1)
        r = c.unsqueeze(-2)-c.unsqueeze(-3) #坐标差,B*N*N*3
        r = r / (r.norm(2,dim=-1,keepdim=True)+1e-8)
        return (attention * r).sum(axis=-2)
        
        
class GCN(nn.Module):

    def __init__(self, nfb, nfa, nf1, nf2, nh1):
        super(GCN, self).__init__()
        self.f0 = Feather_Extract(nfb, nfa)
        self.bn1 = nn.BatchNorm1d(nfb+nfa)
        self.gcn1 = GraphAttentionLayer(nfb+nfa, nf1)
        self.gcn2 = GraphAttentionLayer(nf1, nf2)
        self.bn2 = nn.BatchNorm1d(nf2)
        self.fc1 = nn.Linear(nf2, nh1)
        self.fc2 = nn.Linear(nh1, 1)
    
    def forward(self, d):
        '''d:距离矩阵，b*N*N'''
        D, x = self.f0(d)
        x = x.transpose(1,2)
        x = self.bn1(x)
        x = x.transpose(1,2)
        x = self.gcn1(D, x)
        x = self.gcn2(D, x)
        x = x.transpose(1,2)
        x = self.bn2(x)
        x = x.transpose(1,2)
        x = self.fc1(x)
        x = torch.relu(x)
        x = self.fc2(x)
        x = torch.tanh(x)
        x = x.squeeze(-1)
        x = x.sum(axis=1)
        return x
  
  
class GCNF(nn.Module):
    '''预测力'''

    def __init__(self, nfb, nfa, nf1, nf2):
        super(GCNF, self).__init__()
        self.f0 = Feather_Extract(nfb, nfa)
        self.bn1 = nn.BatchNorm1d(nfb+nfa)
        self.gc1 = GraphAttentionLayer(nfb+nfa, nf1)
        self.gc2 = GraphAttentionLayer(nf1, nf2)
        self.bn2 = nn.BatchNorm1d(nf2)
        self.gcf = ForceGraphAttentionLayer(nf2)
    
    def forward(self, c, d):
        ''' c:坐标差矩阵，b*N*N*3
        d:距离矩阵，b*N*N'''
        D, x = self.f0(d)
        
        x = x.transpose(1,2)
        x = self.bn1(x)
        x = x.transpose(1,2)
        
        x = self.gc1(D, x)
        x = self.gc2(D, x)      
        
        x = x.transpose(1,2)
        x = self.bn2(x)
        x = x.transpose(1,2)
        
        f = self.gcf(D, c, x)
        return f
        
        
if __name__ == '__main__':
    gcnf = GCNF(36,28,80,120)
    d = torch.rand(10,4,4)*4
    c = torch.rand(10,4,3)
    for i in range(10):
        for j in range(4):
            d[i,j,j] = 0
    out = gcnf(c,d)
    print(out.size())
    raise
    
    model = torch.load('gcn4.28.7_epoch50_data5.pth', map_location='cpu')
    #print(model.f0.fsd.coeff, model.f0.fsa.coeff)
    import numpy as np
    data =  np.load(r'..\Si36_10000.npz')
    print('loaded')
    x = data['arr_0'][:5]
    e = data['arr_1'][:5]
    x = torch.tensor(x)
    y = model(x)
    print(y)
    print(e)