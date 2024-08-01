import argparse
import os
import shutil
import subprocess
import sys
import tempfile

import pandas as pd
from anndata import read_h5ad


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
        help="Path to save extracted expression features.",
    )
    parser.add_argument(
        "--model",
        required=True,
        choices=["uce_4", "uce_33"],
        default="uce_4",
        help="Foundation model to use to extract expression features.",
    )
    parser.add_argument(
        "--species",
        required=True,
        help="Species defined by the foundation model. Most commonly `human` or `mouse`.",
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=32,
        help="Batch size for inference.",
    )
    args = parser.parse_args()

    if args.model == "uce_4":
        args.model_loc = "./model_files/4layer_model.torch"
        args.nlayers = 4
    elif args.model == "uce_33":
        args.model_loc = "./model_files/33l_8ep_1024t_1280.torch"
        args.nlayers = 33

    return args


if __name__ == "__main__":
    args = parse_args()
    repo_root = os.path.dirname(__file__)

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = os.path.join(tmp_dir, "")  # trailing slash necessary for UCE pipeline
        adata = read_h5ad(args.input_h5ad)
        adata.obs = pd.DataFrame(index=adata.obs.index)
        del adata.uns
        del adata.obsm
        tmp_h5ad = os.path.join(tmp_dir, "tmp.h5ad")
        adata.write_h5ad(tmp_h5ad)

        subprocess.run(
            args=[
                sys.executable,  # python
                "eval_single_anndata.py",
                "--adata_path",
                tmp_h5ad,
                "--dir",
                tmp_dir,
                "--species",
                args.species,
                "--model_loc",
                args.model_loc,
                "--nlayers",
                str(args.nlayers),
                "--batch_size",
                str(args.batch_size),
            ],
            cwd=os.path.join(repo_root, "models", "UCE"),
        )
        uce_h5ad = os.path.join(tmp_dir, "tmp_uce_adata.h5ad")
        uce_adata = read_h5ad(uce_h5ad)
        uce_adata.obsm["X_expr"] = uce_adata.obsm["X_uce"]
        del uce_adata.obsm["X_uce"]
        uce_adata.write_h5ad(args.output_h5ad)
        print(f"Moved {uce_h5ad} to {args.output_h5ad}")
