from recombinhunt.core.method import GenomeView
from plotly import graph_objects as go
import numpy as np


# REGION COLORS
color_maps = [
    ["#f28266", "#f17063", "#ef5d60", "#ee4f64", "#ec4067", "#d9376d", "#c62d72", "#b32478", "#a01a7d"],
    ["#0077b6", "#0096c7", "#00b4d8", "#48cae4", "#90e0ef", "#ade8f4", "#caf0f8"],
    ["#DDFFBB", "#C7E9B0", "#B3C99C", "#A4BC92"][::-1]
]


def plot_likelihood(genome: GenomeView, xaxis = "coordinates", changes="target") -> go.Figure:
    # input validation
    if not xaxis in ("coordinates", "changes"):
        raise ValueError("The parameter 'coordinates' must be one of ('coordinates', 'changes').")
    if not changes in ("target", "all"):
        raise ValueError("The parameter 'changes' must be one of ('target', 'all').")


    regions = genome.regions

    unique_candidates = [x for x in dict.fromkeys([r.designated for r in regions])]
    region_idx2unique_cand_idx_map = {r_idx: unique_candidates.index(regions[r_idx].designated) for r_idx in
                                      range(len(regions))}
    region_idx2c_map = {r_idx: color_maps[unique_c_idx] for r_idx, unique_c_idx in
                        region_idx2unique_cand_idx_map.items()}


    # XAXIS
    if xaxis == "coordinates":
        x_axis = genome.merged_df.merg_pos              # genomic coordinates
    else:   # if xaxis == "changes"
        if changes == "target":
            x_axis = np.maximum(0, genome.merged_df.seq_change.cumsum() - 1)    # number of changes of t
        else:   # if changes == "all"
            x_axis = np.arange(genome.num_genome_positions)                     # number of changes
        x_axis += 1                                                             # make 1-based

    # MASK CHANGES
    if changes == "all":
        mask = np.ones(genome.num_genome_positions, dtype=bool)     # no_filter
    else:
        mask = genome.merged_df.seq_change                          # filter target

    fig = go.Figure()

    def trace_region_candidate(x_axis, mask, region, candidate, color, shift_y_values=None):
        window_mask = mask[region.pos_start: region.pos_end]

        x_window = x_axis[region.pos_start: region.pos_end]
        x_window_masked = x_window[window_mask]

        y = region.likelihood_values()
        y_window = y[region.pos_start: region.pos_end]
        y_window_masked = y_window[window_mask]

        # shift_y_values so as the first y value always starts from 0
        if shift_y_values is not None:
            y_window_masked = y_window_masked + (
                        shift_y_values - y_window_masked[0 if region.search_dir == 0 else -1])

        fig.add_trace(go.Scatter(
            x=x_window_masked,
            y=y_window_masked,
            name=candidate + (" >>" if region.search_dir == 0 else " <<"),
            mode="lines",
            line={"color": color}
        ))
        return y_window_masked[0]   # TODO improve clarity of code for shifting y

    next_shift = 0
    for r_idx, r in enumerate(regions):
        color = region_idx2c_map[r_idx][0]
        next_shift = trace_region_candidate(x_axis, mask, r, r.designated, color, next_shift)

    # hack to show_legend when only single trace
    if len(fig['data']) == 1:
        fig['data'][0]['showlegend'] = True
        fig['data'][0]['name'] = regions[0].designated + (" >>" if regions[0].search_dir == 0 else " <<")

    # axis_description
    xaxis_title = xaxis
    if changes == "target":
        xaxis_title = "target " + xaxis_title
    xaxis_title = xaxis_title.title()
    fig.update_layout(yaxis_title="Cumulative Log(P)", xaxis_title=xaxis_title)
    fig.update_xaxes(range=[0, np.sum(mask)])

    return fig


def plot_likelihood_whole_genome(genome: GenomeView, xaxis = "coordinates", changes="target") -> go.Figure:
    # input validation
    if not xaxis in ("coordinates", "changes"):
        raise ValueError("The parameter 'coordinates' must be one of ('coordinates', 'changes').")
    if not changes in ("target", "all"):
        raise ValueError("The parameter 'changes' must be one of ('target', 'all').")


    regions = genome.regions

    unique_candidates = [x for x in dict.fromkeys([r.designated for r in regions])]
    region_idx2unique_cand_idx_map = {r_idx: unique_candidates.index(regions[r_idx].designated) for r_idx in
                                      range(len(regions))}
    region_idx2c_map = {r_idx: color_maps[unique_c_idx] for r_idx, unique_c_idx in
                        region_idx2unique_cand_idx_map.items()}

    # XAXIS
    if xaxis == "coordinates":
        x_axis = genome.merged_df.merg_pos              # genomic coordinates
    else:   # if xaxis == "changes"
        if changes == "target":
            x_axis = np.maximum(0, genome.merged_df.seq_change.cumsum() - 1)    # number of changes of t
        else:   # if changes == "all"
            x_axis = np.arange(genome.num_genome_positions)                     # number of changes
        x_axis += 1                                                             # make 1-based

    # MASK CHANGES
    if changes == "all":
        mask = np.ones(genome.num_genome_positions, dtype=bool)     # no_filter
    else:
        mask = genome.merged_df.seq_change                          # filter target

    fig = go.Figure()

    def trace_region_candidate(x_axis, mask, region, candidate, color, shift_y_values=None):

        fig.add_trace(go.Scatter(
            x=x_axis[mask],
            y=region.likelihood_values()[mask],
            name=candidate + (" >>" if region.search_dir == 0 else " <<"),
            mode="lines",
            line={"color": color}
        ))

    next_shift = 0
    for r_idx, r in enumerate(regions):
        color = region_idx2c_map[r_idx][0]
        next_shift = trace_region_candidate(x_axis, mask, r, r.designated, color, next_shift)

    # hack to show_legend when only single trace
    if len(fig['data']) == 1:
        fig['data'][0]['showlegend'] = True
        fig['data'][0]['name'] = regions[0].designated + (" >>" if regions[0].search_dir == 0 else " <<")

    # axis_description
    xaxis_title = xaxis
    if changes == "target":
        xaxis_title = "target " + xaxis_title
    xaxis_title = xaxis_title.title()
    fig.update_layout(yaxis_title="Cumulative Log(P)", xaxis_title=xaxis_title)
    fig.update_xaxes(range=[0, np.sum(mask)])

    return fig
