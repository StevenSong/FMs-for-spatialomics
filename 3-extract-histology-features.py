import argparse
import os

import numpy as np
import pandas as pd
import timm
import torch
from anndata import AnnData, read_h5ad
from PIL import Image
from torch.utils.data import DataLoader, Dataset
from torchvision import transforms
from tqdm import tqdm


class TileDataset(Dataset):
    def __init__(self, tiles: np.ndarray):
        super().__init__()
        self.tiles = tiles
        self.transform = transforms.Compose(
            [
                transforms.Resize(224),
                transforms.ToTensor(),
                transforms.Normalize(
                    mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225)
                ),
            ]
        )

    def __getitem__(self, index) -> Image.Image:
        im = self.tiles[index]
        im = Image.fromarray(im)
        im = self.transform(im)
        return im

    def __len__(self) -> int:
        return len(self.tiles)


def extract_histology(
    *,  # enforce kwargs
    input_h5ad: str,
    output_h5ad: str,
    batch_size: int,
):
    device = "cuda" if torch.cuda.is_available() else "cpu"

    repo_root = os.path.dirname(__file__)
    model_checkpoint = os.path.join(repo_root, "models", "UNI", "pytorch_model.bin")

    model = timm.create_model(
        "vit_large_patch16_224",
        img_size=224,
        patch_size=16,
        init_values=1e-5,
        num_classes=0,
        dynamic_img_size=True,
    )
    model.load_state_dict(torch.load(model_checkpoint, map_location="cpu"), strict=True)
    model = model.eval()
    model = model.to(device)

    adata = read_h5ad(input_h5ad)
    ds = TileDataset(adata.obsm["tiles"])

    embeds = []
    dl = DataLoader(ds, batch_size=batch_size)
    for batch in tqdm(dl):
        ims = batch.to(device)
        with torch.inference_mode():
            embed = model(ims)
        embeds.append(embed.to("cpu"))
    embeds = torch.concat(embeds, dim=0).numpy()

    out_data = AnnData(obs=pd.DataFrame(index=adata.obs.index), obsm={"X_uni": embeds})
    out_data.write(output_h5ad)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_h5ad",
        required=True,
        help="Path to spatial AnnData.",
    )
    parser.add_argument(
        "--output_h5ad",
        required=True,
        help="Path to save extracted histology features.",
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=512,
        help="Batch size for inference.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    extract_histology(
        input_h5ad=args.input_h5ad,
        output_h5ad=args.output_h5ad,
        batch_size=args.batch_size,
    )
