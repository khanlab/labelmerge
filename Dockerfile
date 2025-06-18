# 1) Start from a minimal base that already includes mamba
FROM condaforge/mambaforge:latest

# 2) Set working directory
WORKDIR /src/

# 3) Copy your entire repository into /src/
COPY . /src/

# 4) Disable user‐site packages
ENV PYTHONNOUSERSITE=1


RUN mamba create -n snakebids-env \
    -c conda-forge \
    -c bioconda \
    snakebids -y \
    && mamba clean --all --yes

RUN echo "source /opt/conda/etc/profile.d/conda.sh && conda activate snakebids-env" >> ~/.bashrc

RUN bash -lc "\
    source /opt/conda/etc/profile.d/conda.sh && \
    conda activate snakebids-env && \
    ./labelmerge/run.py \
    test_data/bids_wetrun_testing/tpl-MNI152NLin2009cAsym \
    test_out participant \
    --base-desc 100Parcels7Networks \
    --overlay_bids_dir test_data/bids_wetrun_testing/tpl-MNI152NLin2009cAsym \
    --overlay_desc tn \
    --use-conda \
    --conda-create-envs-only \
    --cores all \
    --conda-prefix /src/conda-envs \
    && mamba clean --all --yes \
    && rm -rf /opt/conda/pkgs /root/.cache \
    "

# 8) Point Snakemake to the correct profile
ENV SNAKEMAKE_PROFILE=/src/labelmerge/workflow/profiles/docker-conda

# 9) Default entrypoint
ENTRYPOINT ["/src/entrypoint.sh"]
