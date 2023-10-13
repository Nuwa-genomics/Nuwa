import torch
import torch.nn as nn


class LinBnDrop(nn.Sequential):
    """Module grouping `BatchNorm1d`, `Dropout` and `Linear` layers, adapted from fastai."""
    
    def __init__(self, n_in, n_out, bn=True, p=0., act=None, lin_first=True):
        layers = [nn.BatchNorm1d(n_out if lin_first else n_in)] if bn else []
        if p != 0: layers.append(nn.Dropout(p))
        lin = [nn.Linear(n_in, n_out, bias=not bn)]
        if act is not None: lin.append(act)
        layers = lin+layers if lin_first else layers+lin
        super().__init__(*layers)


class Encoder(nn.Module):
    """Encoder for CITE-seq data"""
    
    def __init__(self, nfeatures_rna, nfeatures_pro, hidden_rna, hidden_pro, z_dim):
        super().__init__()
        self.nfeatures_rna = nfeatures_rna
        self.nfeatures_pro = nfeatures_pro

        if nfeatures_rna > 0:
            self.encoder_rna = LinBnDrop(nfeatures_rna, hidden_rna, p=0.1, act=nn.LeakyReLU())

        if nfeatures_pro > 0:
            self.encoder_protein = LinBnDrop(nfeatures_pro, hidden_pro, p=0.1, act=nn.LeakyReLU())

        # make sure hidden_rna and hidden_pro are set correctly
        hidden_rna = 0 if nfeatures_rna == 0 else hidden_rna
        hidden_pro = 0 if nfeatures_pro == 0 else hidden_pro
        
        self.encoder = LinBnDrop(hidden_rna + hidden_pro, z_dim, act=nn.LeakyReLU())

    def forward(self, x):
        if self.nfeatures_rna > 0 and self.nfeatures_pro > 0:
            x_rna = self.encoder_rna(x[:, :self.nfeatures_rna])
            x_pro = self.encoder_protein(x[:, self.nfeatures_rna:])
            x = torch.cat([x_rna, x_pro], 1)

        elif self.nfeatures_rna > 0 and self.nfeatures_pro == 0:
            x = self.encoder_rna(x)

        elif self.nfeatures_rna == 0 and self.nfeatures_pro > 0:
            x = self.encoder_protein(x)
            
        return self.encoder(x)


class Decoder(nn.Module):
    """Decoder for CITE-seq data"""
    def __init__(self, nfeatures_rna, nfeatures_pro, hidden_rna, hidden_pro, z_dim):
        super().__init__()
        self.nfeatures_rna = nfeatures_rna
        self.nfeatures_pro = nfeatures_pro

        # make sure hidden_rna and hidden_pro are set correctly
        hidden_rna = 0 if nfeatures_rna == 0 else hidden_rna
        hidden_pro = 0 if nfeatures_pro == 0 else hidden_pro

        hidden = hidden_rna + hidden_pro

        self.decoder = nn.Sequential(
            LinBnDrop(z_dim, hidden, act=nn.LeakyReLU()),
            LinBnDrop(hidden, nfeatures_rna + nfeatures_pro, bn=False)
            )

    def forward(self, x):
        x = self.decoder(x)
        return x

class CiteAutoencoder(nn.Module):
    def __init__(self, nfeatures_rna=0, nfeatures_pro=0, hidden_rna=120, hidden_pro=8, z_dim=20):
        """ Autoencoder for citeseq data """
        super().__init__()
 
        self.encoder = Encoder(nfeatures_rna, nfeatures_pro, hidden_rna, hidden_pro, z_dim)
        self.decoder = Decoder(nfeatures_rna, nfeatures_pro, hidden_rna, hidden_pro, z_dim)


    def forward(self, x):
        x = self.encoder(x)
        x = self.decoder(x)
        return x