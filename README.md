# FMs-for-spatialomics
Tool to apply RNASeq and pathology foundation models to spatial transcriptomics. Currently, this tool only works for extracting expression and histology based features from 10x visium samples and assumes data have been processed by their Spaceranger pipeline.

## Setup
1. Set the environment variable:
    ```bash
    export REPO_ROOT=/path/to/FMs-for-spatialomics
    ```
1. Create conda environment:
    ```bash
    conda env create -f $REPO_ROOT/env.yml
    ```
1. Clone foundation model repos:
    1. UNI:
        * Request access to [UNI on huggingface](https://huggingface.co/MahmoodLab/UNI).
        * Clone the UNI huggingface repo locally (this method uses an [ssh key with huggingface](https://huggingface.co/settings/keys). Alternate methods to clone the repo are also possible e.g. with an access token and/or the HF CLI):
            ```bash
            # make sure git-lfs is already installed
            git clone git@hf.co:MahmoodLab/UNI $REPO_ROOT/models/UNI
            ```
    1. UCE:
        * Clone the UCE github repo locally:
            ```bash
            git clone git@github.com:snap-stanford/UCE.git $REPO_ROOT/models/UCE
            ```
        * Download the model files:
            ```bash
            wget https://figshare.com/ndownloader/articles/24320806/versions/5 -O $REPO_ROOT/models/UCE/model_files/temp.zip
            ```
        * Unzip model files:
            ```bash
            unzip $REPO_ROOT/models/UCE/model_files/temp.zip -d $REPO_ROOT/models/UCE/model_files
            ```
        * Untar additional model files:
            ```bash
            tar -xvf $REPO_ROOT/models/UCE/model_files/protein_embeddings.tar.gz -C $REPO_ROOT/models/UCE/model_files
            ```
## Extract features
Feature extraction of expression and histology data using modality specific foundation models. The following instructions extract features per capture area (which we refer to as a slide) and should be repeated for each slide of interest.
1. Activate conda environemnt:
    ```bash
    conda activate spatialFM
    ```
1. Prepare per-slide data:
    ```bash
    python $REPO_ROOT/1-convert-to-anndata.py \
    --spatial_path /path/to/spaceranger/output/outs \
    --slide_path /path/to/full/resolution/slide.tif \
    --tile_width N \
    --output_h5ad /path/to/converted.h5ad
    ```
    where the `tile_width` should typically be an integer within 1 pixel of `spot_diameter`.
1. Extract expression features:
    ```bash
    python $REPO_ROOT/2-extract-expression-features.py \
    --input_h5ad /path/to/converted.h5ad \
    --output_h5ad /path/to/expr.h5ad \
    --model uce_4 \
    --species human
    ```
    where non-standard `--species` may also need be configured with the UCE model, and `--model` can be one of `uce_4` or `uce_33` to use the 4 or 33 layer UCE models, respectively. 
1. Extract histology features:
    ```bash
    python $REPO_ROOT/3-extract-histology-features.py \
    --input_h5ad /path/to/converted.h5ad \
    --output_h5ad /path/to/hist.h5ad
    ```
1. Unify features:  
Differing inclusion criteria between the foundation models result in minor differences in which barcoded-spots actually get processed. This final step is to take the intersection of those spots for further analysis.
    ```bash
    python $REPO_ROOT/4-combine-data.py \
    --source_h5ad /path/to/converted.h5ad \
    --expr_h5ad /path/to/expr.h5ad \
    --hist_h5ad /path/to/hist.h5ad \
    --output_h5ad /path/to/extracted.h5ad
    ```
1. Clean up (optional):
    ```bash
    rm /path/to/converted.h5ad
    rm /path/to/uce.h5ad
    rm /path/to/uni.h5ad
    ```

### Automating feature extraction
An example script to run all steps of the feature extraction pipeline is located at `run-extract.sh`. It should similarly be run with the conda environment activated.

### Hardware acceleration
The above inference scripts will default to use the first available GPU, if detected. To disable or alter this behavior, the simplest method is to set the environment variable `CUDA_VISIBLE_DEVICES`. If running out of GPU memory, consider tuning the batch size with the `--batch_size` flags to the feature extraction scripts. The default batch sizes were tuned for a single `V100 16GB` GPU.

## Analysis

An example notebook for downstream analysis using extracted features can be found in the `evaluations` subdirectory of this repo.
