import numpy as np
import pandas as pd
import plotly
import plotly.graph_objs as go
from plotly import tools


def diversity(x, index="shannon"):
    """
    Port of vegan's diversity function:

    > x = c(486.0, 123.0, 2.0, 3.0, 1.0, 1.0, 2.0, 2.0, 20.0, 2.0, 2.0, 189.0, 3.0, 4.0, 2.0, 2.0, 3.0, 3.0, 6.0, 2.0, 2.0, 3.0)
    > vegan::diversity(x, index="shannon")
    [1] 1.320956
    > vegan::diversity(x, index="simpson")
    [1] 0.6138655
    > vegan::diversity(x, index="invsimpson")
    [1] 2.589771

    Args:
        x (numpy.array)
        index (Optional[str]): diversity index; accepts string or list of strings;
            choices are shannon, simpson, invsimpson
    Returns:
        float

    >>> import numpy as np
    >>> x = np.array([
            20.0, 3.0, 189.0, 1.0, 2.0, 486.0, 4.0, 2.0, 2.0, 2.0, 3.0, 2.0, 2.0,
            3.0, 2.0, 3.0, 6.0, 123.0, 1.0, 2.0, 2.0, 3.0, 0.0, np.nan
        ])
    >>> diversity(x)
    1.32095...
    >>> diversity(x, "simpson")
    0.61386...
    >>> diversity(x, "invsimpson")
    2.58977...
    """
    valid_indexes = ["shannon", "simpson", "invsimpson"]
    if isinstance(index, list):
        # validate the requested indexes
        for i in index:
            if i.lower() not in valid_indexes:
                raise ValueError(
                    "indexes must be one of 'shannon', 'simpson', or 'invsimpson'"
                )
        H = list()
        # remove nan
        x = x[~np.isnan(x)]
        # remove zeros after removing nans; only applicable to shannon
        x = x[x > 0]
        x = x / sum(x)
        for i in index:
            i = i.lower()
            if i == "shannon":
                ix = -x * np.log(x)
                H.append(sum(ix))
            else:
                ix = x * x
                if i == "simpson":
                    H.append(1 - sum(ix))
                else:
                    H.append(1 / sum(ix))
    else:
        index = index.lower()
        if not index in valid_indexes:
            raise ValueError(
                "index must be one of 'shannon', 'simpson', or 'invsimpson'"
            )
        # remove nan
        x = x[~np.isnan(x)]
        # remove zeros after removing nans; only applicable to shannon
        x = x[x > 0]
        x = x / sum(x)
        if index == "shannon":
            x = -x * np.log(x)
        else:
            x = x * x
        H = sum(x)
        if index == "simpson":
            H = 1 - H
        elif index == "invsimpson":
            H = 1 / H
    return H


def preprocess_table(path, hierarchy, value_cols=None):
    df = pd.read_table(path)
    if not isinstance(hierarchy, list):
        hierarchy = list(hierarchy)
    if not value_cols:
        header = df.columns.tolist()
        header_set = set(header)
        if not len(header) == len(header_set):
            raise ValueError(f"{table} contains a duplicate column header")
        value_cols = [
            i
            for i in header
            if i in header_set - set(hierarchy)
            and i in df.select_dtypes(include=["number"])
        ]
    else:
        for v in value_cols:
            if v not in df.columns:
                raise ValueError(f"{v} is not a column header in {path}")
    for i in hierarchy:
        if i not in df.columns:
            raise KeyError(f"{i} is not a column header in {path}")
    df[hierarchy] = df[hierarchy].fillna("NA")
    df[value_cols] = df[value_cols].fillna(0)
    # drop extraneous columns
    df = df[hierarchy + value_cols]
    return df, value_cols


def get_dfs_across_hierarchy(
    df, hierarchy, value_cols, percent_cutoff=0.01, reorder="shannon", dependent=None
):
    """
    Summarizes counts across each level of the hierarchy.

    Args:
        df (pandas.DataFrame): contains counts in value_cols and
            identifiers in hierarchy columns
        hierarchy (list): column IDs of hierachical levels in the order
            they should be aggregated in
        value_cols (list): list of columns IDs containing data to
            summarize
        percent_cutoff (Optional[float]): omit levels from the
            summaries that fall below this threshold
        reorder (Optional[str]): organize value_cols such that they
            are sorted by specified diversity metric; one of shannon,
            simpson, or invsimpson
        dependent (Optional[str]): levels of the hierarchy should be
            grouped with previous level as their meaning is dependent
            on being coupled with the previous level; specify a string
            separator used when joining the column values
    """
    if percent_cutoff > 99 or percent_cutoff < 0:
        raise ValueError("Percent cutoff should be between 0 and 99.")
    threshold = (percent_cutoff / 100) * df[value_cols].sum().sum()
    dfs = dict()
    sample_order = None
    for i, lvl in enumerate(hierarchy):
        if i == 0:
            # calculates diversity using most specific assignment
            t = df.groupby(hierarchy).sum().reset_index()
            if reorder:
                sample_order = t[value_cols].apply(
                    lambda x: diversity(x, index=reorder)
                )
                sample_order = sample_order.sort_values(ascending=False).keys().tolist()
            else:
                sample_order = value_cols
        # treat each level of the hierarchy as dependent upon the previous
        if dependent and i > 0:
            dep_lvl = dependent.join(hierarchy[: i + 1])
            tdf = df.copy()
            # group the levels
            tdf[dep_lvl] = tdf[hierarchy[: i + 1]].apply(
                lambda x: dependent.join(x.astype(str)), axis=1
            )
            # sum across grouped levels
            t = tdf.groupby(dep_lvl).sum().reset_index()
            # rename the grouped column back into the level ID
            t.columns = [lvl] + value_cols
        # each level of the hierarchy can be quantified independently
        else:
            t = df.groupby(lvl).sum().reset_index()
        t["sum_of_values"] = t[value_cols].sum(axis=1)
        # remove low abundant traces from the plot
        t = t[t["sum_of_values"] > threshold]
        # sort df by abundance which organizes the traces in the plots
        t.sort_values(by="sum_of_values", ascending=False, inplace=True)
        # drop the extra column
        t = t[[lvl] + sample_order]
        # get the values for this particular level as they
        # become column labels after transposing
        col_labels = t[lvl].tolist()
        # calculate relative abundances
        ra_t = t[value_cols].div(t[value_cols].sum())
        ra_t = ra_t * 100
        ra_t = ra_t.transpose()
        ra_t.columns = col_labels

        t = t.transpose()
        # set the column names as the taxonomy levels
        t.columns = col_labels
        # remove the hierarchy column and its values from the first row (after transpose)
        t.drop([lvl], inplace=True)
        dfs[lvl] = [t, ra_t]
    return dfs


def get_abundance_figure_from_dfs(
    dfs,
    hierarchy,
    title=None,
    colors=None,
    height=700,
    showlegend=False,
    hoverinfo="text+y",
):
    if not isinstance(colors, list):
        colors = None
    if not colors:
        colors = [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
        ]
    colors_len = len(colors)
    fig = tools.make_subplots(rows=2, cols=1, shared_xaxes=True, print_grid=False)
    trace_visibility = list()
    for i, lvl in enumerate(hierarchy):
        for j, table in enumerate(dfs[lvl], start=1):
            for k, col in enumerate(table.columns.tolist()):
                trace = go.Bar(
                    x=table.index,
                    y=table[col],
                    name=col,
                    text=col,
                    hoverinfo=hoverinfo,
                    # always start showing the first level
                    visible=True if i == 0 else False,
                    marker=dict(color=colors[k % colors_len]),
                )
                fig.append_trace(trace, j, 1)
                trace_visibility.append(lvl)

    title = title if title else "Abundance Summary Across Hierarchy"
    buttons = []
    for lvl in hierarchy:
        button = dict(
            label=lvl.title(),
            method="update",
            args=[dict(visible=[i == lvl for i in trace_visibility])],
        )
        buttons.append(button)
    updatemenus = list(
        [
            dict(
                active=0,
                buttons=buttons,
                x=0,
                y=1.1,
                pad={"l": 0, "t": 20},
                direction="down",
                xanchor="left",
                yanchor="top",
            )
        ]
    )
    annotations = list(
        [
            dict(
                text="<b>Level:</b>",
                x=0,
                y=1.1,
                xref="paper",
                yref="paper",
                align="left",
                showarrow=False,
            )
        ]
    )
    fig["layout"].update(
        {
            "title": title,
            "yaxis": {"title": "Counts"},
            "yaxis2": {"title": "Relative Abundance (%)"},
            "showlegend": showlegend,
            "height": height,
            "barmode": "stack",
            "hovermode": "closest",
        }
    )
    fig["layout"]["updatemenus"] = updatemenus
    fig["layout"]["annotations"] = annotations
    return fig


def abundance_figure(
    path,
    hierarchy,
    value_cols=None,
    percent_cutoff=0.01,
    reorder="shannon",
    title=None,
    height=700,
    colors=None,
    dependent=None,
    hoverinfo="text+y",
):
    """
    Convenience function to carry out each step and build a figure from an input table.
    """
    df, value_cols = preprocess_table(path, hierarchy, value_cols)
    dfs = get_dfs_across_hierarchy(
        df, hierarchy, value_cols, percent_cutoff, reorder, dependent
    )
    fig = get_abundance_figure_from_dfs(
        dfs,
        hierarchy,
        title=title,
        colors=colors,
        height=height,
        showlegend=False,
        hoverinfo=hoverinfo,
    )
    return fig


def calculate_diversity(
    df, hierarchy=None, value_cols=None, index=["shannon", "simpson", "invsimpson"]
):
    """
    Calculate diversity across a hierachical level using the values of
    value_cols.

    Args:
        df (pandas.DataFrame)
        hierarchy (Optional[list]): the default is every row (df.index)
        value_cols (Optional[list]): the default is df.columns, though this
            assumes all columns have int or float values
        index (Optional[list]): the diversity indexes to calculate

    Returns:
        pandas.DataFrame
    """
    t = df.copy()
    t = t.groupby(hierarchy).sum().reset_index()
    t = t[value_cols].apply(lambda x: diversity(x, index=index))
    t = pd.DataFrame(t.values.tolist(), index=t.index, columns=index)
    return t
