import argparse

import numpy as np
import scanpy
from anndata import AnnData
from PIL import Image
from tqdm import tqdm


def convert(
    *,  # enforce kwargs
    spatial_path: str,
    slide_path: str,
    tile_width: int,
) -> AnnData:
    """
    Convert Visium data to AnnData format and generate histology tiles.

    Parameters
    ----------
    spatial_path : str
        Description of parameter `x`.
    slide_path: str
        Description of parameter `y` (with type not specified).
    tile_width: int
        Descrip

    Returns
    -------
    AnnData
        Spatial data in the same format as `scanpy.read_visium` with the following differences:
        - `library_id` is simply `slide`.
        - `obs` has additional columns `pxl_col/res_in_fullres`. The values are the same as those in `obsm["spatial"]`.
        - `obsm["tiles"]` contains histology tiles.

    See Also
    --------
    https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_visium.html
    """
    adata = scanpy.read_visium(
        path=spatial_path,
        count_file="filtered_feature_bc_matrix.h5",
        library_id="slide",
        load_images=True,
    )
    im = np.asarray(Image.open(slide_path))
    assert im.dtype == np.uint8

    adata.obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]] = adata.obsm["spatial"]

    pos_df = adata.obs
    scale_factors = adata.uns["spatial"]["slide"]["scalefactors"]

    # check that the actual spot size is not smaller than tile size - 1
    # ok if actual spot is larger than tile size
    spot_radius = tile_width // 2
    actual_size = scale_factors["spot_diameter_fullres"]
    assert actual_size - spot_radius * 2 > -1

    count = 0
    padded = 0
    tile_shape = (spot_radius * 2, spot_radius * 2, 3)

    tiles = []
    for row, col in zip(
        pos_df["pxl_row_in_fullres"], tqdm(pos_df["pxl_col_in_fullres"], desc="Tile")
    ):
        count += 1
        row_lo, row_hi = row - spot_radius, row + spot_radius
        col_lo, col_hi = col - spot_radius, col + spot_radius
        tile = im[row_lo:row_hi, col_lo:col_hi]
        if tile.shape != tile_shape:
            new_tile = np.full(tile_shape, 255, dtype=tile.dtype)
            new_tile[: tile.shape[0], : tile.shape[1], : tile.shape[2]] = tile
            tile = new_tile
            padded += 1
        tiles.append(tile)
    tiles = np.stack(tiles)

    adata.obsm["tiles"] = tiles
    print(f"Generated {count} tiles ({padded} were padded).")

    return adata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--spatial_path",
        required=True,
        help="Path to `outs` folder of Spaceranger output for a single slide.",
    )
    parser.add_argument(
        "--slide_path",
        required=True,
        help="Path to full resolution histology image.",
    )
    parser.add_argument(
        "--tile_width",
        type=int,
        required=True,
        help=(
            "Tile width in pixels, centered on each spatial spot. "
            "Tile width is typically close to the `spot_diameter`."
        ),
    )
    parser.add_argument(
        "--output_h5ad",
        required=True,
        help="Path to save final AnnData object.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    adata = convert(
        spatial_path=args.spatial_path,
        slide_path=args.slide_path,
        tile_width=args.tile_width,
    )
    adata.write_h5ad(args.output_h5ad)
