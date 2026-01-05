# nf-bcftools-par-merge
Parallelized bcftools merge in Nextflow for large datasets (>1,000 samples). Speed up over a naive `bcftools merge` is proportional to the number of samples & size of VCFs. Work in progress, initially developed for merging tumor WGS VCFs from GDC. Edit sample sheet parsing and the PREP process to adapt for other datasets.

**DAG**
```mermaid
flowchart TB
  subgraph " "
    subgraph params
      v12["min_windows"]
      v11["ref_index"]
      v0["n_samples"]
      v24["out_prefix"]
      v4["sample_sheet"]
      v13["primary_assembly"]
    end
    v7([PREP])
    v14([WINDOW_GENOME])
    v19([FIRST_MERGE])
    v21([FINAL_MERGE])
    v25([CONCAT])
    subgraph publish
      v28["bcf"]
      v29["index"]
    end
    v0 --> v7
    v4 --> v7
    v11 --> v14
    v12 --> v14
    v13 --> v14
    v0 --> v19
    v7 --> v19
    v14 --> v19
    v0 --> v21
    v19 --> v21
    v21 --> v25
    v24 --> v25
    v25 --> v28
    v25 --> v29
  end

```
