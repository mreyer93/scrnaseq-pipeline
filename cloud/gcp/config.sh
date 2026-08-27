#!/usr/bin/env bash
# Shared settings for the cloud/gcp/*.sh scripts. Edit these, then run create_disk.sh
# and start_vm.sh. The other scripts source this file automatically.

PROJECT_ID="${PROJECT_ID:-your-gcp-project-id}"
ZONE="${ZONE:-us-central1-a}"

VM_NAME="${VM_NAME:-scrnaseq-vm}"
DISK_NAME="${DISK_NAME:-scrnaseq-data}"

# 16 vCPU / 64 GB. Two things drive the requirement: building a splici index for
# simpleaf (~16-32 GB), and holding a large AnnData object in memory during clustering.
# A 500k-cell experiment wants considerably more RAM than a laptop has.
MACHINE_TYPE="${MACHINE_TYPE:-n2-standard-16}"

# FASTQs dominate the disk requirement; size for your largest expected batch.
DISK_SIZE_GB="${DISK_SIZE_GB:-500}"
DISK_TYPE="${DISK_TYPE:-pd-balanced}"

BOOT_DISK_SIZE_GB="${BOOT_DISK_SIZE_GB:-50}"
IMAGE_FAMILY="${IMAGE_FAMILY:-debian-12}"
IMAGE_PROJECT="${IMAGE_PROJECT:-debian-cloud}"
