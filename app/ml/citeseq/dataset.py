from torch.utils.data import Dataset


class TabularDataset(Dataset):
    """Custome dataset for tabular data"""
    def __init__(self, data):
        self.data = data

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        x = self.data[idx]
        return x, x