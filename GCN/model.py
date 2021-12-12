import torch
from torch import nn
from torch.nn import functional as F


class args:
    hidden = [256, 256, 256, 256, 256, 256, 256]


def sparse_dropout(x, rate, noise_shape):
    """

    :param x:
    :param rate:
    :param noise_shape: int scalar
    :return:
    """
    random_tensor = 1 - rate
    random_tensor += torch.rand(noise_shape).to(x.device)
    dropout_mask = torch.floor(random_tensor).byte()
    i = x._indices()
    v = x._values()
    i = i[:, dropout_mask]
    v = v[dropout_mask]
    out = torch.sparse.FloatTensor(i, v, x.shape).to(x.device)
    out = out * (1. / (1-rate))
    return out


class GraphConvolution(nn.Module):

    def __init__(self, input_dim, output_dim, num_features_nonzero,
                 dropout=0.,
                 is_sparse_inputs=False,
                 bias=True,
                 activation=F.relu,
                 featureless=False):
        super(GraphConvolution, self).__init__()
        self.dropout = dropout
        self.bias = bias
        self.activation = activation
        self.is_sparse_inputs = is_sparse_inputs
        self.featureless = featureless
        self.num_features_nonzero = num_features_nonzero
        # 双通道 output_dim = (2,*)
        self.weight = nn.Parameter(torch.randn(input_dim, *output_dim))
        if bias:
            self.bias = nn.Parameter(torch.zeros(output_dim))

    def forward(self, inputs):
        # print('inputs:', inputs)
        # support : sparse adjacent matrix
        x, support = inputs
        if self.training and self.is_sparse_inputs:
            x = sparse_dropout(x, self.dropout, self.num_features_nonzero)
        elif self.training:
            x = F.dropout(x, self.dropout)
        # convolve
        if not self.featureless:  # if it has features x
            if self.is_sparse_inputs:
                xw = torch.sparse.mm(x, self.weight)
            else:
                xw = torch.mm(x, self.weight)
        else:
            xw = self.weight
        out = torch.sparse.mm(support, xw)
        if self.bias is not None:
            out += self.bias
        return self.activation(out), support


class GCN(nn.Module):

    def __init__(self, input_dim, output_dim, num_features_nonzero):
        super(GCN, self).__init__()
        self.input_dim = input_dim
        self.output_dim = output_dim
        print('input dim:', input_dim)
        print('output dim:', output_dim)
        print('num_features_nonzero:', num_features_nonzero)
        self.layers = nn.Sequential(GraphConvolution(self.input_dim, args.hidden[0], num_features_nonzero,
                                                     activation=F.relu,
                                                     dropout=0.5,
                                                     is_sparse_inputs=True),
                                    *[GraphConvolution(hidden, hidden, num_features_nonzero,
                                                       activation=F.relu,
                                                       dropout=0.5,
                                                       is_sparse_inputs=True) for hidden in args.hidden],
                                    )
        self.out = nn.Linear(args.hidden[-1], output_dim)

    def forward(self, inputs):
        x, support = inputs
        # GCN
        x = self.layers((x, support))
        x = x[0]
        x = [F.max_pool1d(i, i.size(2)).squeeze(2) for i in x]
        x = torch.cat(x, 1)
        x = F.relu(x[0])
        x = self.out(x)
        return x

    def l2_loss(self):
        layer = self.layers.children()
        layer = next(iter(layer))
        loss = None
        for p in layer.parameters():
            if loss is None:
                loss = p.pow(2).sum()
            else:
                loss += p.pow(2).sum()
        return loss
