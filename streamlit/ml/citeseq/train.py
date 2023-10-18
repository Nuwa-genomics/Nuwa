import torch
import torch.nn as nn
from collections import defaultdict
from tqdm import tqdm

def train_model(model, train_dl, valid_dl, lr, epochs, device, callback_on_epoch=None, verbose=True):
    model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    scheduler = torch.optim.lr_scheduler.OneCycleLR(optimizer, max_lr=lr,
                                                    steps_per_epoch=len(train_dl), epochs=epochs)
    criterion = nn.MSELoss()
    history = defaultdict(list)
    
    mean_losses = []
    for epoch in tqdm(range(1, epochs + 1)):
        model.train()
        
        train_loss = 0
        nsamples_train = 0
        for x, y in train_dl:
            optimizer.zero_grad()

            # Forward pass
            x_prime = model(x.to(device))
            
            loss = criterion(x_prime, y.to(device))

            # Backward pass
            loss.backward()
            optimizer.step()
            scheduler.step()
            
            # log losses
            batch_size = x.shape[0]
            nsamples_train += batch_size
            train_loss += batch_size*(loss.item())
            
        valid_loss = 0
        nsamples_valid = 0
        
        model = model.eval()
        with torch.no_grad():
            for x, y in valid_dl:
                x_prime = model(x.to(device))

                loss = criterion(x_prime, x.to(device))
                
                # log losses
                batch_size = x.shape[0]
                nsamples_valid += batch_size
                valid_loss += batch_size*(loss.item())
                
        train_loss = train_loss / nsamples_train
        valid_loss = valid_loss / nsamples_valid

        history['train'].append(train_loss)
        history['valid'].append(valid_loss)
        
        if verbose and epoch%10==0:
            print(f'Epoch {epoch}: train loss {train_loss}; valid loss {valid_loss}')

        if callback_on_epoch is not None:
            callback_on_epoch(epoch, train_loss, valid_loss)

    return model, history

def get_encodings(model, dl, device):
    model.eval()
    with torch.no_grad():
        encodings = [model.encoder(x.to(device)) for x, _ in dl]
    return torch.cat(encodings, dim=0)