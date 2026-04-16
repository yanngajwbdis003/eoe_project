import os
import scanpy as sc
import numpy as np
import scipy.sparse as sp
import pickle
from cnmf import cNMF

INPUT_ADATA_PATH = '../data/eoe_sub_raw_approx.h5ad'   
OUTPUT_BASE_DIR  = './cnmf_output'
SEED             = 42
NUM_HIGHVAR_GENES = 2000

CONFIGS = {
    'stromal': {
        'labels': ['Fibroblast', 'Fibroblast (cycling)'],
        'k_values':  [4, 5, 6, 7],
        'n_iter':    30,          
        'final_k':   5,           
        'density_threshold': 0.05,
    },
    'epithelial': {
        'labels': ['Quiescent basal cell', 'Basal cell (cycling)',
                   'Suprabasal', 'Apical cell'],
        'k_values':  [4, 5, 6, 7, 8],
        'n_iter':    30,
        'final_k':   6,
        'density_threshold': 0.05,
    },
    'full': {
        'labels': None,           
        'k_values':  [12,15,18,20,23],
        'n_iter':    30,
        'final_k':   20,
        'density_threshold': 0.10,
    },
}

def get_raw_counts(adata):
    if adata.raw is not None:
        X = adata.raw.X
        var = adata.raw.var
    else:
        X = adata.X
        var = adata.var
        sample = X[:100, :100]
        if sp.issparse(sample):
            sample = sample.toarray()
        if not np.allclose(sample, np.round(sample)):
            raise ValueError(
                "make sure if it's raw count"
            )
    return X, var


def prepare_subset_h5ad(adata_full, labels, run_name, output_dir):
    if labels is None:
        print(f"    use all cells: {adata_full.shape[0]:,}")
        adata_sub = adata_full.copy()
    else:
        adata_sub = adata_full[
            adata_full.obs['cell_type_anno'].isin(labels)
        ].copy()
        print(f"    extract {adata_sub.shape[0]:,} cells")
        print(f"    condition:\n{adata_sub.obs['disease_status'].value_counts().to_string()}")

    X_raw, var_raw = get_raw_counts(adata_sub)

    adata_counts = sc.AnnData(
        X=X_raw,
        obs=adata_sub.obs.copy(),
        var=var_raw.copy()
    )

    # 
    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f'{run_name}_counts.h5ad')
    adata_counts.write_h5ad(save_path)
    print(f"    saving counts: {save_path}")
    return save_path, adata_sub


def run_cnmf(run_name, counts_path, k_values, n_iter, final_k,
             density_threshold, output_dir):

    print(f"\n{'='*60}")
    print(f"start cNMF: {run_name}")
    print(f"K={k_values}, n_iter={n_iter}, final_k={final_k}")
    print(f"{'='*60}")

    cnmf_obj = cNMF(output_dir=output_dir, name=run_name)

    print(">>> Step 1: Prepare...")
    cnmf_obj.prepare(
        counts_fn=counts_path,
        components=k_values,
        n_iter=n_iter,
        seed=SEED,
        num_highvar_genes=NUM_HIGHVAR_GENES
    )

    print(">>> Step 2: Factorize）")
    cnmf_obj.factorize(worker_i=0, total_workers=1)

    print(">>> Step 3: Combine replicates...")
    cnmf_obj.combine()

    print(">>> Step 4: K selection plot...")
    cnmf_obj.k_selection_plot()

    # Step 5: consensus for final_k
    print(f">>> Step 5: Consensus (K={final_k}, density_threshold={density_threshold})...")
    cnmf_obj.consensus(
        k=final_k,
        density_threshold=density_threshold,
        show_clustering=True
    )

    # Step 6: load + save results
    print(f">>> Step 6: save the results (K={final_k})...")
    usage, spectra_scores, spectra_tpm, top_genes = cnmf_obj.load_results(
        K=final_k, density_threshold=density_threshold
    )

    results_dir = os.path.join(output_dir, f'results_k{final_k}')
    os.makedirs(results_dir, exist_ok=True)

    usage.to_csv(os.path.join(results_dir, 'usage.csv'))
    spectra_scores.to_csv(os.path.join(results_dir, 'spectra_scores.csv'))
    spectra_tpm.to_csv(os.path.join(results_dir, 'spectra_tpm.csv'))
    top_genes.to_csv(os.path.join(results_dir, 'top_genes.csv'))

    with open(os.path.join(output_dir, f'{run_name}_cnmf_obj.pkl'), 'wb') as f:
        pickle.dump(cnmf_obj, f)

    print(f"    {results_dir}")
    return cnmf_obj, usage, spectra_scores, top_genes


# ==========================================
# 3. main flow
# ==========================================
print("=" * 60)
adata = sc.read_h5ad(INPUT_ADATA_PATH)
print(f"dataset: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
print(f"Conditions: {adata.obs['disease_status'].value_counts().to_dict()}")

results = {}

for run_name, cfg in CONFIGS.items():
    print(f"\n{'#'*60}")
    print(f"# start: {run_name.upper()}")
    print(f"{'#'*60}")

    run_output_dir = os.path.join(OUTPUT_BASE_DIR, run_name)

    counts_path, adata_sub = prepare_subset_h5ad(
        adata_full=adata,
        labels=cfg['labels'],
        run_name=run_name,
        output_dir=run_output_dir
    )

    adata_sub.write_h5ad(os.path.join(run_output_dir, f'adata_{run_name}.h5ad'))

    cnmf_obj, usage, spectra_scores, top_genes = run_cnmf(
        run_name=run_name,
        counts_path=counts_path,
        k_values=cfg['k_values'],
        n_iter=cfg['n_iter'],
        final_k=cfg['final_k'],
        density_threshold=cfg['density_threshold'],
        output_dir=run_output_dir
    )

    results[run_name] = {
        'cnmf_obj': cnmf_obj,
        'usage': usage,
        'spectra_scores': spectra_scores,
        'top_genes': top_genes,
    }

    print(f"\n>>> {run_name.upper()} Top 10 genes per GEP (K={cfg['final_k']}):")
    for gep in top_genes.columns:
        print(f"  GEP{gep}: {', '.join(top_genes[gep].head(10).tolist())}")

print("\n" + "=" * 60)
for run_name, cfg in CONFIGS.items():
    run_output_dir = os.path.join(OUTPUT_BASE_DIR, run_name)
    print(f"  {run_name:12s} → {run_output_dir}/results_k{cfg['final_k']}/")
