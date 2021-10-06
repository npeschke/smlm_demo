import gzip
import multiprocessing as mp
import pathlib as pl
import pickle
import shutil

import pandas as pd
from smlm import config as smlm_config
from smlm import helpers
from smlm import plot
from smlm import stats
from smlm.orte import Orte
from smlm.styles import prism

prism.prism_style()

recalculate_voronoi = True

stage_method = smlm_config.STAGE_METHOD
method_display = "Eye based, 5 Stages"

result_dir = pl.Path("results/")
result_dir_fig5 = result_dir.joinpath("fig5")

helpers.assure_dir(result_dir_fig5)

data_dir = pl.Path("data/")
orte_dir = data_dir.joinpath("orte/")

zipped_orte_paths = [path for path in orte_dir.glob(f"*.csv.gz")]

for zipped_orte_path in zipped_orte_paths:
    with gzip.open(zipped_orte_path, "rb") as zip_in:
        with orte_dir.joinpath(zipped_orte_path.stem).open("wb") as csv_out:
            shutil.copyfileobj(zip_in, csv_out)

raw_orte_paths = [path for path in orte_dir.glob(f"*.csv")]
labelling_csv_path = data_dir.joinpath("labelling.csv")

# Labelling
stages = pd.read_csv(labelling_csv_path,
                     names=["file", "manual_stage", "confocal_stage", "automatic_stage",
                            "manual_5_stage", "ilastik_stage", "conf_ccp", "labelling", "include"],
                     dtype={"manual_stage": "Int64", "confocal_stage": "Int64",
                            "automatic_stage": "Int64", "manual_5_stage": "Int64",
                            "ilastik_stage": "Int64", "include": bool})

# Find input files
orte_paths = [raw_orte_path for raw_orte_path in raw_orte_paths if raw_orte_path.name in stages.file.to_list()]

voronoi_orte_file_path = result_dir.joinpath("voronoi_orte.p")
if voronoi_orte_file_path.exists() and not recalculate_voronoi:
    raw_all_cells, raw_all_clusters = pickle.load(voronoi_orte_file_path.open("rb"))

else:
    with mp.Pool(processes=min(len(orte_paths), mp.cpu_count())) as p:
        # Calculating Voronoi tesselation for all input files
        orte_list = p.map(Orte, orte_paths, chunksize=1)

    raw_all_cells = pd.DataFrame()
    raw_all_clusters = pd.DataFrame()
    for orte in orte_list:
        raw_all_cells = raw_all_cells.append(orte.orte_df)
        raw_all_clusters = raw_all_clusters.append(orte.get_named_cluster_meta("dbscan"))

    with voronoi_orte_file_path.open("wb") as vof:
        pickle.dump((raw_all_cells, raw_all_clusters), vof)

    del orte
    del orte_list

all_cells = helpers.merge_localizations_stages(raw_all_cells, stages)
all_clusters = helpers.merge_localizations_stages(raw_all_clusters, stages)

del raw_all_cells
del raw_all_clusters

# Density histogram by apoptotic stage
plot_col = smlm_config.LOG_NORM_DENS_COL
plot_label = smlm_config.LOG_NORM_DENS_LABEL

density_lims = smlm_config.DENSITY_LIM

n_stages = plot.get_n_stages(stage_method, stages)

fig, ax = plot.plot_stage_histogram(plot_col=plot_col, method=stage_method, data=all_cells,
                                    plot_label=plot_label, method_label=method_display,
                                    density_lims=density_lims,
                                    n_stages=n_stages)
ax.margins(0)

fig.savefig(result_dir_fig5.joinpath(f"stage_{plot_col}_{stage_method}.{smlm_config.FIG_TYPE}"))

# Radius - Density histogram by apoptotic stage
plot.plot_r_density_dist_stages_separate(
    data=all_cells,
    method=stage_method,
    stat=smlm_config.HIST_STAT,
    result_dir=result_dir_fig5,
    plot_col_density=plot_col,
    plot_label_density=plot_label,
)

# Descriptive stage data
stats.describe_hist_distribution(
    all_cells,
    desc_col=plot_col,
    group=stage_method
).to_csv(result_dir.joinpath(f"voronoi_{plot_col}_{stage_method}_filtered.csv"))

# Wilcoxon rank sums on stage density differences
wilcox_stats = stats.wilcoxon_rank_sums(plot_col, n_stages, stage_method, data=all_cells)
wilcox_stats.to_csv(result_dir.joinpath(f"wilcoxon_rank_sums_{plot_col}_{stage_method}.csv"),
                    index=False, float_format="%1.10f")
