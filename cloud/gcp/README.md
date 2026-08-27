# Running on GCP

Use this for the `fastq` input mode (quantifying raw reads with simpleaf/alevin-fry), or
when a dataset is too large to cluster on a laptop. Building a human splici index needs
roughly 16-32 GB RAM; large experiments need more memory for the clustering stage than a
laptop has.

If you already have count matrices and a normal-sized experiment, you probably do not need
this - the local config handles that case comfortably.

## Prerequisites

- GCP project with billing enabled and the Compute Engine API on
- [`gcloud` CLI](https://cloud.google.com/sdk/docs/install), authenticated (`gcloud init`)
- Edit `cloud/gcp/config.sh` and set `PROJECT_ID`

## One-time setup

```bash
./cloud/gcp/create_disk.sh    # persistent disk for references, FASTQs, results
./cloud/gcp/start_vm.sh       # spot VM, disk mounted at /mnt/data
gcloud compute ssh scrnaseq-vm --project=<project> --zone=us-central1-a
```

On the VM (interactively, so your GitHub credentials never enter VM metadata):

```bash
git clone https://github.com/<you>/scrnaseq-pipeline.git
cd scrnaseq-pipeline

wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p /mnt/data/miniforge3
source /mnt/data/miniforge3/etc/profile.d/conda.sh
mamba env create -f envs/python.yml -n scrnaseq
mamba env create -f envs/environment.yml -n scrnaseq-quant
mamba env create -f envs/r.yml -n scrnaseq-r

# reference for the FASTQ path (Ensembl human shown)
mkdir -p /mnt/data/reference && cd /mnt/data/reference
REL=112
curl -O https://ftp.ensembl.org/pub/release-$REL/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
curl -O https://ftp.ensembl.org/pub/release-$REL/gtf/homo_sapiens/Homo_sapiens.GRCh38.$REL.gtf.gz
gunzip *.gz
```

## Routine workflow

```bash
./cloud/gcp/start_vm.sh
gcloud compute ssh scrnaseq-vm --project=<project> --zone=us-central1-a

# on the VM
source /mnt/data/miniforge3/etc/profile.d/conda.sh && conda activate scrnaseq
cd scrnaseq-pipeline
python3 scripts/samplesheet.py /mnt/data/run1/samples.csv fastq
snakemake --configfile /mnt/data/run1/config.yaml --use-conda --cores 16

gsutil -m cp -r /mnt/data/run1/results gs://<your-bucket>/run1/

# back on your laptop
./cloud/gcp/stop_vm.sh
```

Set `chemistry` in the config to match the kit (`10xv2`, `10xv3`, `10xv4`, ...). Getting
this wrong produces almost no valid barcodes rather than an error, so check it against
the library prep before a long run.

### Preemption

Spot VMs can be reclaimed at any time. Restart with `start_vm.sh` and re-run the same
`snakemake` command - each stage writes its own `.h5ad`, so it resumes from the last
completed stage rather than starting over. Add `--rerun-incomplete` if a stage was
interrupted mid-write.

## Cost

`us-central1`, verify against the [pricing calculator](https://cloud.google.com/products/calculator):

- **Compute**: `n2-standard-16` spot ~ $0.23/hr, billed only while running.
- **Storage**: 500 GB `pd-balanced` ~ $50/month, billed **whether or not the VM runs**.

For occasional use, storage dominates. Either tear down between projects, or keep FASTQs
and references in a GCS bucket (~4-5x cheaper per GB) and copy them onto a fresh disk per run.

## Teardown

```bash
gcloud compute instances delete scrnaseq-vm --project=<project> --zone=us-central1-a
gcloud compute disks delete scrnaseq-data --project=<project> --zone=us-central1-a
```
