"""Optional plotly figures (install with: pip install "mlcseq[plot]").

Functions return a figure instead of calling .show(), so they work in notebooks,
scripts (fig.write_html) and tests alike.
"""
from __future__ import annotations

import pandas as pd


def _px():
    try:
        import plotly.express as px
        import plotly.graph_objects as go
    except ImportError as e:  # pragma: no cover
        raise ImportError('Plotting needs plotly: pip install "mlcseq[plot]"') from e
    return px, go


def zones(dfs: list[pd.DataFrame], names=None, y: str = "RT", title: str | None = None):
    """Scatter of mass vs RT (or Vol) for one or more compound sets."""
    px, _ = _px()
    names = names or [f"set_{i + 1}" for i in range(len(dfs))]
    data = pd.concat([d.assign(type=n) for d, n in zip(dfs, names)])
    fig = px.scatter(data, x="Mass", y=y, color="type", title=title)
    return fig


def base_calls(calls, y: str = "RT", title: str | None = None):
    """Compounds plus a line + label for every base call (BaseCalls tuple)."""
    _, go = _px()
    df, pairs = calls
    fig = go.Figure(go.Scatter(x=df.Mass, y=df[y], mode="markers", name="compounds"))
    for heavy, light, name in pairs:
        pair = df[df.Mass.isin([heavy, light])]
        fig.add_trace(go.Scatter(x=pair.Mass, y=pair[y], mode="lines+markers",
                                 line=dict(color="pink"), showlegend=False))
        fig.add_annotation(x=pair.Mass.mean(), y=pair[y].mean(), yshift=5, text=name,
                           showarrow=False)
    if title:
        fig.update_layout(title=title)
    return fig


def lasso_selector(df: pd.DataFrame, selected: list):
    """Interactive Jupyter widget: lasso/click points to separate 5' and 3' zones.

    Selected row positions are appended to `selected`; afterwards use
    `df.iloc[selected]` and save it (e.g. to CSV) so the choice is reproducible.
    """
    _, go = _px()
    fig = go.FigureWidget(go.Scatter(x=df.Mass, y=df.RT, mode="markers"))
    scatter = fig.data[0]
    scatter.marker.color = ["#6371f2"] * len(df)

    def on_pick(trace, points, selector):
        selected.extend(points.point_inds)
        colors = list(scatter.marker.color)
        for i in points.point_inds:
            colors[i] = "#87a14a"
        scatter.marker.color = colors

    scatter.on_selection(on_pick)
    scatter.on_click(on_pick)
    return fig
