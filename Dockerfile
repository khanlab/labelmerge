FROM condaforge/miniforge3:latest

WORKDIR /src/

# Copy your code
COPY . /src/

# Disable user site packages
ENV PYTHONNOUSERSITE=1

# Use bash for the following RUNs
SHELL ["/bin/bash", "-c"]

# ---- ONE SINGLE RUN ----
RUN set -e && \
    conda install -n base -c conda-forge mamba -y && \
    mamba create -y -n snakebids-env -c conda-forge -c bioconda snakebids unzip && \
    source /opt/conda/etc/profile.d/conda.sh && \
    conda activate snakebids-env && \
    ./labelmerge/run.py test_data/bids_wetrun_testing/tpl-MNI152NLin2009cAsym test_out participant \
            --base-desc 100Parcels7Networks --overlay_bids_dir test_data/bids_wetrun_testing/tpl-MNI152NLin2009cAsym \
            --overlay_desc tn -np --use-conda --conda-create-envs-only --conda-prefix /src/conda-envs && \
    conda clean --all -y && \
    rm -rf /opt/conda/pkgs /root/.caches

# Set snakemake profile
ENV SNAKEMAKE_PROFILE=/src/labelmerge/workflow/profiles/docker-conda

# Set entrypoint
ENTRYPOINT ["/src/entrypoint.sh"]
