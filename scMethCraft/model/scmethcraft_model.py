import math
from typing import List

import torch
import torch.nn as nn
import numpy as np
from .layers import ConvLayer, DenseLayer
from .utils_model import m_round

def sigmoid(x):
    s = 1 / (1 + np.exp(-x))
    return s

def _col_round(x):
    frac = x - math.floor(x)
    if frac <= 0.5:
        return math.floor(x)
    return math.ceil(x)


def m_get_filter_dim(seq_length: int, pooling_sizes: List[int]):
    filter_dim = seq_length
    for ps in pooling_sizes:
        filter_dim = _col_round(filter_dim / ps)
    return int(filter_dim/2)


class model_of_scBasset(nn.Module):
    """
    PytTorch implementation of scBasset model (Yuan and Kelley, 2022)
    Article link: https://www.nature.com/articles/s41592-022-01562-8
    Original implementation in Keras: https://github.com/calico/scBasset

    Args:
        n_cells: number of cells to predict region accessibility
        n_filters_init: nr of filters for the initial conv layer
        n_repeat_blocks_tower: nr of layers in the convolutional tower
        filters_mult: proportion by which the nr of filters should inrease in the
            convolutional tower
        n_bottleneck_layer: size of the bottleneck layer
        batch_norm: whether to apply batch norm across model layers
        dropout: dropout rate across layers, by default we don't do it for
            convolutional layers but we do it for the dense layers
    """

    def __init__(
        self,
        n_cells: int,
        n_filters_init: int = 288,
        n_repeat_blocks_tower: int = 5,
        filters_mult: float = 1.122,
        n_filters_pre_bottleneck: int = 256,
        n_bottleneck_layer: int = 32,
        batch_norm: bool = True,
        dropout: float = 0.0,
        genomic_seq_length: int = 10000,
    ):
        super().__init__()

        self.stem = ConvLayer(
            in_channels=4,
            out_channels=n_filters_init,
            kernel_size=17,
            pool_size=3,
            dropout=dropout,
            batch_norm=batch_norm,
        )

        tower_layers = []
        curr_n_filters = n_filters_init
        for i in range(n_repeat_blocks_tower):
            tower_layers.append(
                ConvLayer(
                    in_channels=curr_n_filters,
                    out_channels=_round(curr_n_filters * filters_mult),
                    kernel_size=5,
                    pool_size=2,
                    dropout=dropout,
                    batch_norm=batch_norm,
                )
            )
            curr_n_filters = _round(curr_n_filters * filters_mult)
        self.tower = nn.Sequential(*tower_layers)

        self.pre_bottleneck = ConvLayer(
            in_channels=curr_n_filters,
            out_channels=n_filters_pre_bottleneck,
            kernel_size=1,
            dropout=dropout,
            batch_norm=batch_norm,
            pool_size=2,
        )

        # get pooling sizes of the upstream conv layers
        pooling_sizes = [3] + [2] * n_repeat_blocks_tower + [1]
        # get filter dimensionality to account for variable sequence length
        filter_dim = _get_filter_dim(
            seq_length=genomic_seq_length, pooling_sizes=pooling_sizes
        )
        self.bottleneck = DenseLayer(
            in_features=n_filters_pre_bottleneck * filter_dim,
            out_features=n_bottleneck_layer,
            use_bias=True,
            batch_norm=True,
            dropout=0.2,
            activation_fn=nn.Identity(),
        )
        self.final = nn.Linear(n_bottleneck_layer, n_cells)

    def forward(
        self,
        x: torch.Tensor,  # input shape: (batch_size, 4, seq_length)
    ):
        # TODO: add random shift to act as a regularizer on the dataset level
        # TODO: add use reverse complement randomly on the dataset level
        x = self.stem(x)
        x = self.tower(x)
        x = self.pre_bottleneck(x)
        # flatten the input
        x = x.view(x.shape[0], -1)
        x = self.bottleneck(x)
        x = self.final(x)
        return x

class EarlyStopping:
    """Early stops the training if validation loss doesn't improve after a given patience."""
    def __init__(self, patience=10, verbose=False, delta=0):
        """
        Args:
            patience (int): How long to wait after last time validation loss improved.
                            Default: 7
            verbose (bool): If True, prints a message for each validation loss improvement. 
                            Default: False
            delta (float): Minimum change in the monitored quantity to qualify as an improvement.
                            Default: 0
        """
        self.patience = patience
        self.verbose = verbose
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.val_loss_min = np.Inf
        self.delta = delta

    def __call__(self, val_loss):

        score = -val_loss

        if self.best_score is None:
            self.best_score = score
        elif score < self.best_score + self.delta:
            self.counter += 1
            print(f'EarlyStopping counter: {self.counter} out of {self.patience}')
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.best_score = score
            self.counter = 0    





