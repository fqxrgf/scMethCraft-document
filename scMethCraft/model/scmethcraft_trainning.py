from .scmethcraft_model import *
from .utils_model import *

import torch
import pandas as pd
import sklearn.metrics as metrics
import scanpy as sc
import anndata as ad
from scipy.special import expit
import sys
import numpy

import random
device = "cuda:2"

def load_seq(input_path,filename,load_range,mode= "onehot"):
    if mode == "onehot":
        with h5py.File(input_path+filename) as file:
            if load_range != False:
                onehot = file["X"][load_range[0]:load_range[1]]
            else:
                onehot = file["X"]
            onehot = torch.nn.functional.one_hot(torch.tensor(onehot).to(torch.int64)).transpose(-2, -1)
            return onehot.to(torch.float)
        
    if mode == "kmer":
        with h5py.File(input_path+filename) as file:
            if load_range != False:
                kmer = file["Kmer"][load_range[0]:load_range[1]]
            else:
                kmer = file["Kmer"]           
            return torch.tensor(kmer).to(torch.float)
        
    if mode == "both":
        with h5py.File(input_path+filename) as file:
            if load_range != False:
                onehot = file["X"][load_range[0]:load_range[1]]
                kmer = file["Kmer"][load_range[0]:load_range[1]]
                pos  = file["Pos"][load_range[0]:load_range[1]]
            else:
                onehot = file["X"]
                kmer = file["Kmer"]  
                pos  = file["Pos"]
                
            onehot = torch.nn.functional.one_hot(torch.tensor(onehot).to(torch.int64)).transpose(-2, -1)
            return onehot.to(torch.float),torch.tensor(kmer).to(torch.float),torch.tensor(pos)

class MethyDataset(torch.utils.data.Dataset):
    def __init__(self, onehot,kmer, state,pos):
        self.onehot = onehot
        self.kmer = kmer
        self.state = state
        self.pos = pos
    def __getitem__(self, index):
        return self.onehot[index], self.state[index],self.kmer[index],self.pos[index]
    def __len__(self):
        return self.state.shape[0]  
    
    
class SimilarityLayer(nn.Module):
    def __init__(
        self,
        n_cells: int,
        dropout_rate: float = 0.1, 
        batch_norm = True
    ):
        
        super().__init__()
        self.similarity_matrix=torch.nn.Parameter(torch.rand(n_cells,n_cells))
        self.n_cells = n_cells
        self.alpha = 0
        self.dropout_rate = dropout_rate
        self.dropout = self.dropout_rate
        self.eyematrix = 1-torch.eye(self.n_cells,self.n_cells, device=device)
        self.batch_norm = (
            nn.BatchNorm1d(n_cells) if batch_norm else nn.Identity()
        )
        
    def my_dropout(self,input_matrix):
        dropout_matrix = (torch.rand(self.n_cells,self.n_cells, device=device)>self.dropout).float()
        dropout_matrix = torch.mul(dropout_matrix,self.eyematrix)
        input_matrix = torch.mul(dropout_matrix, input_matrix)
        input_matrix = input_matrix/(1-self.dropout)
        if self.alpha>0:
            input_matrix = input_matrix+self.alpha*self.n_cells*torch.mean(input_matrix)*torch.eye(self.n_cells,self.n_cells, device=device)
        return input_matrix
    
        
    def forward(
        self,
        input_vector: torch.Tensor,  

    ):

        fixed_similarity_matrix = self.similarity_matrix
        fixed_similarity_matrix = torch.abs(fixed_similarity_matrix+fixed_similarity_matrix.T)/2
        fixed_similarity_matrix = self.my_dropout(fixed_similarity_matrix)        
        output_vector = torch.matmul(input_vector,fixed_similarity_matrix)/(1+self.alpha)
        output_vector = self.batch_norm(output_vector)  
        return output_vector  

class Sequence_extraction(nn.Module):
    def __init__(
        self,
        n_cells: int,
        K:int = 4**8,  
        n_filters_init: int = 256,
        n_repeat_blocks_tower: int =2,
        filters_mult: float = 1.41421,
        n_filters_pre_bottleneck: int = 256,
        n_bottleneck_layer: int = 25,
        dropout_rate_similarity: float = 0.3, 
        batch_norm: bool = True,
        embedding_dim: int = 16,
        dropout: float = 0.0,
        genomic_seq_length: int = 10000,    ):
        super().__init__()

        self.stem = ConvLayer(
            in_channels=4,
            out_channels=n_filters_init,
            kernel_size=12,
            pool_size=4,
            dropout=dropout,
            batch_norm=batch_norm,
        )

        
        
        tower_layers = []
        curr_n_filters = n_filters_init
        
        for i in range(n_repeat_blocks_tower):
            tower_layers.append(
                ConvLayer(
                    in_channels=curr_n_filters,
                    out_channels=m_round(curr_n_filters * filters_mult),
                    kernel_size=5,
                    pool_size=2,
                    dropout=dropout,
                    batch_norm=batch_norm,
                )
            )
            curr_n_filters = m_round(curr_n_filters * filters_mult)
            
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
        pooling_sizes = [4] + [2] * n_repeat_blocks_tower + [1]
        # get filter dimensionality to account for variable sequence length
        filter_dim = m_get_filter_dim(seq_length=genomic_seq_length, pooling_sizes=pooling_sizes)
        
        self.hidden_onehot = DenseLayer(
            in_features=n_filters_pre_bottleneck * filter_dim,
            out_features=n_bottleneck_layer,
            use_bias=True,
            batch_norm=True,
            dropout=0.2,
            activation_fn=nn.Identity(),
        )
        
        self.hidden_pos_1 = KANLinear(
            in_features=64,
            out_features=32,
        )
    
        self.hidden_pos_2 = KANLinear(
            in_features=32,
            out_features=n_bottleneck_layer*2,
        )
        self.Kmer1 = DenseLayer(
            in_features=4 ** K,
            out_features=4 ** int(K/2),
            use_bias=True,
            batch_norm=True,
            dropout=0.2,
            activation_fn=nn.GELU(),
        )
        self.pos_embedding_attention = nn.Embedding(256,embedding_dim)
        self.hidden_Kmer = DenseLayer(
            in_features=4**int(K/2),
            out_features=n_bottleneck_layer,
            use_bias=True,
            batch_norm=True,
            dropout=0.2,
            activation_fn=nn.GELU(),
        )
        
        self.kmer_embedding_linear = DenseLayer(
            in_features=1,
            out_features=embedding_dim,
            use_bias=True,
            batch_norm=False,
            dropout=0,
            activation_fn=nn.GELU(),
        )
        
        self.transformer = torch.nn.TransformerEncoderLayer(embedding_dim,4,batch_first=True)
        
        self.kmer_embedding_linear_inverse = DenseLayer(
            in_features=embedding_dim,
            out_features=1,
            use_bias=True,
            batch_norm=False,
            dropout=0,
            activation_fn=nn.GELU(),
        )
        
        self.final = nn.Linear(n_bottleneck_layer*2, n_cells)
 
    def forward(
        self,
        onehot: torch.Tensor,  
        kmer: torch.Tensor,
        pos:torch.Tensor
    ):

        
        onehot = self.stem(onehot)
        onehot = self.tower(onehot)
        onehot = self.pre_bottleneck(onehot)
        onehot = onehot.view(onehot.shape[0], -1)
        onehot = self.hidden_onehot(onehot)
        
        kmer = self.Kmer1(kmer)
        res_kmer = kmer
        kmer = self.kmer_embedding_linear(kmer.unsqueeze(2))
        kmer_seq_embedding = kmer + self.pos_embedding_attention.weight      
        kmer_seq_embedding = self.transformer(kmer_seq_embedding) 
        kmer = self.kmer_embedding_linear_inverse(kmer_seq_embedding).squeeze(2)
        kmer = kmer + res_kmer
        
        
        kmer = self.hidden_Kmer(kmer)
        
        pos = self.hidden_pos_1(pos)
        pos = self.hidden_pos_2(pos)
        
        latent = torch.cat((onehot,kmer),1)
        
        res_latent = latent
        latent = torch.mul(latent, pos)
        latent = latent + res_latent
        
        latent = self.final(latent)

        
        return latent
    
class Similarity_weighting(nn.Module):
    def __init__(
        self,
        n_cells: int,
        dropout_rate: float = 0.1, 
        batch_norm = True
        ):
        super().__init__()
        self.SimilarityLayer1 = SimilarityLayer(n_cells,dropout_rate,batch_norm)

    def train(self):
        self.SimilarityLayer1.train()
        self.SimilarityLayer1.dropout = self.SimilarityLayer1.dropout_rate
        self.SimilarityLayer1.alpha = 0
            
    def eval(self):
        self.SimilarityLayer1.eval()
        self.SimilarityLayer1.dropout = 0
        self.SimilarityLayer1.alpha = 0.2
        
    def forward(self, methy_level_vector: torch.Tensor):
        methy_level_vector = self.SimilarityLayer1(methy_level_vector)
        return methy_level_vector    

class KANLinear(torch.nn.Module):
    def __init__(
        self,
        in_features,
        out_features,
        grid_size=5,
        spline_order=3,
        scale_noise=0.1,
        scale_base=1.0,
        scale_spline=1.0,
        enable_standalone_scale_spline=True,
        base_activation=torch.nn.SiLU,
        grid_eps=0.02,
        grid_range=[-1, 1],
    ):
        super(KANLinear, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.grid_size = grid_size
        self.spline_order = spline_order

        h = (grid_range[1] - grid_range[0]) / grid_size
        grid = (
            (
                torch.arange(-spline_order, grid_size + spline_order + 1) * h
                + grid_range[0]
            )
            .expand(in_features, -1)
            .contiguous()
        )
        self.register_buffer("grid", grid)

        self.base_weight = torch.nn.Parameter(torch.Tensor(out_features, in_features))
        self.spline_weight = torch.nn.Parameter(
            torch.Tensor(out_features, in_features, grid_size + spline_order)
        )
        if enable_standalone_scale_spline:
            self.spline_scaler = torch.nn.Parameter(
                torch.Tensor(out_features, in_features)
            )

        self.scale_noise = scale_noise
        self.scale_base = scale_base
        self.scale_spline = scale_spline
        self.enable_standalone_scale_spline = enable_standalone_scale_spline
        self.base_activation = base_activation()
        self.grid_eps = grid_eps

        self.reset_parameters()

    def reset_parameters(self):
        torch.nn.init.kaiming_uniform_(self.base_weight, a=math.sqrt(5) * self.scale_base)
        with torch.no_grad():
            noise = (
                (
                    torch.rand(self.grid_size + 1, self.in_features, self.out_features)
                    - 1 / 2
                )
                * self.scale_noise
                / self.grid_size
            )
            self.spline_weight.data.copy_(
                (self.scale_spline if not self.enable_standalone_scale_spline else 1.0)
                * self.curve2coeff(
                    self.grid.T[self.spline_order : -self.spline_order],
                    noise,
                )
            )
            if self.enable_standalone_scale_spline:
                # torch.nn.init.constant_(self.spline_scaler, self.scale_spline)
                torch.nn.init.kaiming_uniform_(self.spline_scaler, a=math.sqrt(5) * self.scale_spline)

    def b_splines(self, x: torch.Tensor):
        """
        Compute the B-spline bases for the given input tensor.

        Args:
            x (torch.Tensor): Input tensor of shape (batch_size, in_features).

        Returns:
            torch.Tensor: B-spline bases tensor of shape (batch_size, in_features, grid_size + spline_order).
        """
        assert x.dim() == 2 and x.size(1) == self.in_features

        grid: torch.Tensor = (
            self.grid
        )  # (in_features, grid_size + 2 * spline_order + 1)
        x = x.unsqueeze(-1)
        bases = ((x >= grid[:, :-1]) & (x < grid[:, 1:])).to(x.dtype)
        for k in range(1, self.spline_order + 1):
            bases = (
                (x - grid[:, : -(k + 1)])
                / (grid[:, k:-1] - grid[:, : -(k + 1)])
                * bases[:, :, :-1]
            ) + (
                (grid[:, k + 1 :] - x)
                / (grid[:, k + 1 :] - grid[:, 1:(-k)])
                * bases[:, :, 1:]
            )

        assert bases.size() == (
            x.size(0),
            self.in_features,
            self.grid_size + self.spline_order,
        )
        return bases.contiguous()

    def curve2coeff(self, x: torch.Tensor, y: torch.Tensor):
        """
        Compute the coefficients of the curve that interpolates the given points.

        Args:
            x (torch.Tensor): Input tensor of shape (batch_size, in_features).
            y (torch.Tensor): Output tensor of shape (batch_size, in_features, out_features).

        Returns:
            torch.Tensor: Coefficients tensor of shape (out_features, in_features, grid_size + spline_order).
        """
        assert x.dim() == 2 and x.size(1) == self.in_features
        assert y.size() == (x.size(0), self.in_features, self.out_features)

        A = self.b_splines(x).transpose(
            0, 1
        )  # (in_features, batch_size, grid_size + spline_order)
        B = y.transpose(0, 1)  # (in_features, batch_size, out_features)
        solution = torch.linalg.lstsq(
            A, B
        ).solution  # (in_features, grid_size + spline_order, out_features)
        result = solution.permute(
            2, 0, 1
        )  # (out_features, in_features, grid_size + spline_order)

        assert result.size() == (
            self.out_features,
            self.in_features,
            self.grid_size + self.spline_order,
        )
        return result.contiguous()

    @property
    def scaled_spline_weight(self):
        return self.spline_weight * (
            self.spline_scaler.unsqueeze(-1)
            if self.enable_standalone_scale_spline
            else 1.0
        )

    def forward(self, x: torch.Tensor):
        assert x.dim() == 2 and x.size(1) == self.in_features

        base_output = torch.nn.functional.linear(self.base_activation(x), self.base_weight)
        spline_output = torch.nn.functional.linear(
            self.b_splines(x).view(x.size(0), -1),
            self.scaled_spline_weight.view(self.out_features, -1),
        )
        return base_output + spline_output

    @torch.no_grad()
    def update_grid(self, x: torch.Tensor, margin=0.01):
        assert x.dim() == 2 and x.size(1) == self.in_features
        batch = x.size(0)

        splines = self.b_splines(x)  # (batch, in, coeff)
        splines = splines.permute(1, 0, 2)  # (in, batch, coeff)
        orig_coeff = self.scaled_spline_weight  # (out, in, coeff)
        orig_coeff = orig_coeff.permute(1, 2, 0)  # (in, coeff, out)
        unreduced_spline_output = torch.bmm(splines, orig_coeff)  # (in, batch, out)
        unreduced_spline_output = unreduced_spline_output.permute(
            1, 0, 2
        )  # (batch, in, out)

        # sort each channel individually to collect data distribution
        x_sorted = torch.sort(x, dim=0)[0]
        grid_adaptive = x_sorted[
            torch.linspace(
                0, batch - 1, self.grid_size + 1, dtype=torch.int64, device=x.device
            )
        ]

        uniform_step = (x_sorted[-1] - x_sorted[0] + 2 * margin) / self.grid_size
        grid_uniform = (
            torch.arange(
                self.grid_size + 1, dtype=torch.float32, device=x.device
            ).unsqueeze(1)
            * uniform_step
            + x_sorted[0]
            - margin
        )

        grid = self.grid_eps * grid_uniform + (1 - self.grid_eps) * grid_adaptive
        grid = torch.concatenate(
            [
                grid[:1]
                - uniform_step
                * torch.arange(self.spline_order, 0, -1, device=x.device).unsqueeze(1),
                grid,
                grid[-1:]
                + uniform_step
                * torch.arange(1, self.spline_order + 1, device=x.device).unsqueeze(1),
            ],
            dim=0,
        )

        self.grid.copy_(grid.T)
        self.spline_weight.data.copy_(self.curve2coeff(x, unreduced_spline_output))

    def regularization_loss(self, regularize_activation=1.0, regularize_entropy=1.0):
        """
        Compute the regularization loss.

        This is a dumb simulation of the original L1 regularization as stated in the
        paper, since the original one requires computing absolutes and entropy from the
        expanded (batch, in_features, out_features) intermediate tensor, which is hidden
        behind the F.linear function if we want an memory efficient implementation.

        The L1 regularization is now computed as mean absolute value of the spline
        weights. The authors implementation also includes this term in addition to the
        sample-based regularization.
        """
        l1_fake = self.spline_weight.abs().mean(-1)
        regularization_loss_activation = l1_fake.sum()
        p = l1_fake / regularization_loss_activation
        regularization_loss_entropy = -torch.sum(p * p.log())
        return (
            regularize_activation * regularization_loss_activation
            + regularize_entropy * regularization_loss_entropy
        )    
    
    
def output_model(MethyBasset_part1,MethyBasset_part2,savepath = f"../sample_data/output/"):
    
    import os
    if not os.path.exists(savepath):
        os.makedirs(savepath)
        print("Folder created")
    else:
        print("Folder already exists")


    torch.save(MethyBasset_part1.state_dict(), f"{savepath}/scMethCraft_part1.pth")
    torch.save(MethyBasset_part2.state_dict(), f"{savepath}/scMethCraft_part2.pth")