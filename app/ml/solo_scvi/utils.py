import pandas as pd


def get_solo_model_history(solo, vae):

    # solo
    solo_valid_loss_df = pd.DataFrame(columns=['epoch', 'validation_loss'])
    solo_train_loss_df = pd.DataFrame(columns=['epoch', 'train_loss_epoch'])
    for i in solo.history:
        if i == "validation_loss":
            solo_valid_loss_df = solo.history[i]
        if i == "train_loss_epoch":
            solo_train_loss_df = solo.history[i]


    solo_df = solo_valid_loss_df.join(solo_train_loss_df, on='epoch')

    # vae
    vae_train_loss_epoch_df = pd.DataFrame(columns=['train_loss_epoch'])
    vae_reconstruction_loss_train_df = pd.DataFrame(columns=['reconstruction_loss_train'])

    for i in vae.history:
        if i == "train_loss_epoch":
            vae_train_loss_epoch_df[i] = vae.history[i]

        if i == "reconstruction_loss_train":
            vae_reconstruction_loss_train_df[i] = vae.history[i]
                
    vae_df = vae_train_loss_epoch_df.join(vae_reconstruction_loss_train_df)

    return solo_df, vae_df