"""HTML report generation and marker–cell graph export for CellTyper."""
from __future__ import annotations

import base64
import html
import io
import os
import re
import math
import textwrap
import tempfile
import sys
import warnings
from datetime import datetime, timezone
from typing import Any, List, Optional, Tuple

# If matplotlib is missing, warn only once (many clusters per run).
_mpl_import_logged: bool = False


def _graph_stderr(msg: str) -> None:
    pass


# Matplotlib figure size / resolution for marker–cell graph exports.
_MARKER_GRAPH_FIGSIZE_IN: Tuple[float, float] = (15.0, 9.5)
_GRAPH_DPI = 200

# Muted, publication-style palette (avoid bright blue / orange).
_GRAPH_MARKER_FILL = "#5a7d8c"
_GRAPH_MARKER_EDGE = "#2f4a54"
_GRAPH_CELL_FILL = "#b8a48c"
_GRAPH_CELL_EDGE = "#6d5f52"
_GRAPH_EDGE_COLOR = "#b8c0c8"
_GRAPH_BG = "#f9f8f6"


def _prepare_graph_visual_attrs(graph) -> Tuple[List[str], List[str]]:
    """Gene (Marker) vs cell type (Cell): distinct colors and shapes for igraph plot/SVG."""
    n = graph.vcount()
    colors = ["#94a3b8"] * n
    shapes = ["circle"] * n
    if "type" not in graph.vs.attributes():
        return colors, shapes
    for i in range(n):
        t = graph.vs[i]["type"]
        if t == "Marker":
            colors[i] = _GRAPH_MARKER_FILL
            shapes[i] = "diamond"
        elif t == "Cell":
            colors[i] = _GRAPH_CELL_FILL
            shapes[i] = "rect"
    return colors, shapes


def _marker_cell_plot_kw(graph) -> dict:
    """Per-vertex size, frame, label styling so genes read as hubs vs compact cell-type boxes."""
    n = graph.vcount()
    dense = n > 28
    sizes = [22] * n
    frame_w = [1.0] * n
    frame_c = ["#64748b"] * n
    lbl_sz = [7] * n
    lbl_c = ["#334155"] * n
    if "type" not in graph.vs.attributes():
        return {
            "vertex_size": sizes,
            "vertex_frame_width": frame_w,
            "vertex_frame_color": frame_c,
            "vertex_label_size": lbl_sz,
            "vertex_label_color": lbl_c,
        }
    for i in range(n):
        t = graph.vs[i]["type"]
        if t == "Marker":
            sizes[i] = 42 if dense else 48
            frame_w[i] = 2.8
            frame_c[i] = _GRAPH_MARKER_EDGE
            lbl_sz[i] = 11 if dense else 12
            lbl_c[i] = "#2f4a54"
        elif t == "Cell":
            sizes[i] = 16 if dense else 19
            frame_w[i] = 1.0
            frame_c[i] = _GRAPH_CELL_EDGE
            lbl_sz[i] = 7.5 if dense else 8.5
            lbl_c[i] = "#4a4238"
    return {
        "vertex_size": sizes,
        "vertex_frame_width": frame_w,
        "vertex_frame_color": frame_c,
        "vertex_label_size": lbl_sz,
        "vertex_label_color": lbl_c,
    }


def _marker_graph_figsize(n: int) -> Tuple[float, float]:
    if n <= 24:
        return (16.0, 11.0)
    scale = 1.0 + 0.018 * (n - 24)
    return (min(16.0 * scale, 24.0), min(11.0 * scale, 17.0))


def _layout_for_report(graph):
    n = graph.vcount()
    repulserad = max(60.0, n * 14.0 + 40.0)
    niter = min(12000, 500 + n * 100)
    try:
        return graph.layout(
            "fr",
            niter=niter,
            repulserad=repulserad,
            maxdelta=max(n * 2, 12),
        )
    except TypeError:
        try:
            return graph.layout("fr", niter=min(8000, 300 + n * 80))
        except Exception:
            pass
    try:
        return graph.layout("fr")
    except Exception:
        try:
            return graph.layout_fr()
        except Exception:
            return graph.layout("auto")


def _layout_contract_factor(n: int) -> float:
    if n <= 18:
        return 0.78
    if n <= 32:
        return 0.9
    return 1.0


def _label_display_text(graph, index: int, raw: str) -> Optional[str]:
    lab = str(raw).strip()
    if not lab:
        return None
    vt = _vertex_type_at(graph, index)
    n = graph.vcount()
    dense = n > 24
    if vt == "Cell":
        wrap_at = 36 if dense else 44
        if len(lab) > wrap_at:
            return textwrap.fill(lab, width=wrap_at, break_long_words=False)
        return lab
    return lab


def _outward_unit_from_centroid(
    node_x: float, node_y: float, cx: float, cy: float, span: float
) -> Tuple[float, float]:
    dx, dy = node_x - cx, node_y - cy
    nr = math.hypot(dx, dy)
    if nr < span * 0.008:
        return 1.0, 0.0
    return dx / nr, dy / nr


def _label_alignment_for_outward_vector(ux: float, uy: float) -> Tuple[str, str]:
    if abs(ux) >= abs(uy):
        return ("left" if ux >= 0 else "right", "center")
    return ("center", "bottom" if uy >= 0 else "top")


def _plan_graph_labels(
    graph,
    layout,
    labels: List[str],
    span: float,
    pk: dict,
) -> List[dict]:
    """Place each label beside its node (small outward offset, no global repulsion)."""
    n = graph.vcount()
    xs = [float(layout[i][0]) for i in range(n)]
    ys = [float(layout[i][1]) for i in range(n)]
    mx = sum(xs) / max(n, 1)
    my = sum(ys) / max(n, 1)

    specs: List[dict] = []

    for i in range(n):
        display = _label_display_text(graph, i, labels[i])
        if not display:
            continue
        vt = _vertex_type_at(graph, i)
        ux, uy = _outward_unit_from_centroid(xs[i], ys[i], mx, my, span)
        node_r = span * (0.036 if vt == "Marker" else 0.021)
        gap = span * 0.006
        offset = node_r + gap
        lx = xs[i] + ux * offset
        ly = ys[i] + uy * offset
        ha, va = _label_alignment_for_outward_vector(ux, uy)
        try:
            fs = float(pk["vertex_label_size"][i])
        except (TypeError, ValueError, KeyError, IndexError):
            fs = 8.0
        try:
            lc = pk["vertex_label_color"][i]
        except (KeyError, IndexError):
            lc = "#334155"
        specs.append(
            {
                "index": i,
                "text": display,
                "node_x": xs[i],
                "node_y": ys[i],
                "x": lx,
                "y": ly,
                "ha": ha,
                "va": va,
                "fontsize": fs,
                "color": lc,
                "weight": "bold" if vt == "Marker" else "normal",
                "type": vt,
            }
        )

    return specs


def _try_adjusttext_positions(ax, specs: List[dict]) -> None:
    """Disabled: adjustText drifts labels away from their nodes on dense graphs."""
    return


def _contract_layout_toward_centroid(layout, factor: float = 0.8):
    """Shrink coordinates toward the centroid so vertex labels stay away from the clip rect."""
    if layout is None:
        return layout
    try:
        n = len(layout)
    except Exception:
        return layout
    if n == 0:
        return layout
    coords = [layout[i] for i in range(n)]
    mx = sum(c[0] for c in coords) / n
    my = sum(c[1] for c in coords) / n
    try:
        out = layout.copy()
        for i in range(n):
            x, y = out[i]
            out[i] = (factor * (x - mx) + mx, factor * (y - my) + my)
        return out
    except Exception:
        import igraph as ig

        newc = [
            (factor * (coords[i][0] - mx) + mx, factor * (coords[i][1] - my) + my)
            for i in range(n)
        ]
        try:
            return ig.Layout(newc)
        except Exception:
            return layout


def _vertex_labels(graph) -> Optional[list]:
    try:
        if "name" not in graph.vs.attributes():
            return None
        return [str(x) if x is not None else "" for x in graph.vs["name"]]
    except Exception:
        return None


def _matplotlib_include_outside_labels(ax, fig, pad_frac: float = 0.35) -> None:
    """Widen axes using text bboxes (labels sit outside data coords) and add margin."""
    ax.set_aspect("equal")
    for ch in ax.get_children():
        if hasattr(ch, "set_clip_on"):
            try:
                ch.set_clip_on(False)
            except Exception:
                pass
    fig.canvas.draw()
    try:
        renderer = fig.canvas.get_renderer()
    except Exception:
        renderer = None
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    if renderer is not None:
        for t in ax.texts:
            try:
                bb = t.get_window_extent(renderer=renderer)
                bb = bb.expanded(1.12, 1.14)
                bb_d = bb.transformed(ax.transData.inverted())
                xmin = min(xmin, bb_d.x0)
                xmax = max(xmax, bb_d.x1)
                ymin = min(ymin, bb_d.y0)
                ymax = max(ymax, bb_d.y1)
            except Exception:
                continue

    xr = xmax - xmin
    yr = ymax - ymin
    if xr <= 0:
        xr = 1.0
    if yr <= 0:
        yr = 1.0
    pad_x = xr * pad_frac
    pad_y = yr * pad_frac
    ax.set_xlim(xmin - pad_x, xmax + pad_x)
    ax.set_ylim(ymin - pad_y, ymax + pad_y)


def _matplotlib_push_edges_under_vertices(ax) -> None:
    """Ensure edges are drawn below scatter markers (igraph often paints lines on top)."""
    for ln in ax.lines:
        try:
            ln.set_zorder(1)
            ln.set_clip_on(False)
        except Exception:
            pass
    for c in ax.collections:
        try:
            c.set_zorder(2)
            c.set_clip_on(False)
        except Exception:
            pass


def _svg_bytes_expand_viewbox(svg_bytes: bytes, pad_ratio: float = 0.15) -> bytes:
    """Widen SVG viewBox so labels near borders are not clipped by the viewport."""
    try:
        s = svg_bytes.decode("utf-8")
    except Exception:
        return svg_bytes
    m = re.search(r'viewBox="([-\d.eE]+)\s+([-\d.eE]+)\s+([-\d.eE]+)\s+([-\d.eE]+)"', s)
    if not m:
        return svg_bytes
    x, y, w, h = map(float, m.groups())
    dw = w * pad_ratio
    dh = h * pad_ratio
    new = f'viewBox="{x - dw} {y - dh} {w + 2 * dw} {h + 2 * dh}"'
    return s.replace(m.group(0), new, 1).encode("utf-8")


def _matplotlib_savefig_tight_including_markers(ax, fig, buf: io.BytesIO, fmt: str) -> None:
    save_kw: dict = {
        "format": fmt,
        "bbox_inches": "tight",
        "facecolor": "white",
        "pad_inches": 4.2,
    }
    extras: List[Any] = []
    extras.extend(ax.texts)
    extras.extend(ax.collections)
    extras.extend(ax.patches)
    extras.extend(ax.lines)
    if extras:
        save_kw["bbox_extra_artists"] = extras
    fig.savefig(buf, **save_kw)


def _mpl_scatter_marker_for_igraph_shape(shape: str) -> str:
    if shape == "diamond":
        return "D"
    if shape == "rect":
        return "s"
    return "o"


def _vertex_type_at(graph, index: int) -> Optional[str]:
    try:
        if "type" in graph.vs.attributes():
            return graph.vs[index]["type"]
    except Exception:
        pass
    return None


def _matplotlib_overlay_vertex_markers(
    ax,
    graph,
    layout,
    vcolors: List[str],
    vshapes: List[str],
    vertex_sizes: List[int],
) -> None:
    """Draw markers with matplotlib; igraph sometimes omits or hides vertex bodies on Agg."""
    n = graph.vcount()
    if n == 0:
        return
    xs = [layout[i][0] for i in range(n)]
    ys = [layout[i][1] for i in range(n)]
    z = 25
    for i in range(n):
        sz = float(vertex_sizes[i])
        area = max(49.0, (sz * 2.1) ** 2)
        vt = _vertex_type_at(graph, i)
        if vt == "Marker":
            ec, lw = "#2f4a54", 1.1
        elif vt == "Cell":
            ec, lw = "#6d5f52", 0.65
        else:
            ec, lw = "#334155", 0.7
        mk = _mpl_scatter_marker_for_igraph_shape(vshapes[i])
        ax.scatter(
            xs[i],
            ys[i],
            s=area,
            c=vcolors[i],
            marker=mk,
            zorder=z,
            edgecolors=ec,
            linewidths=lw,
            clip_on=False,
        )


def _matplotlib_savefig_graph_bytes(fig, ax, buf: io.BytesIO, fmt: str) -> None:
    """Save figure; include patches/lines/collections/text for tight bbox."""
    fig.canvas.draw()
    extras: List[Any] = []
    extras.extend(ax.texts)
    extras.extend(ax.patches)
    extras.extend(ax.collections)
    extras.extend(ax.lines)
    kw: dict = {
        "format": fmt,
        "dpi": _GRAPH_DPI,
        "facecolor": "white",
        "bbox_inches": "tight",
        "pad_inches": 0.55,
    }
    if extras:
        kw["bbox_extra_artists"] = extras
    fig.savefig(buf, **kw)


def _draw_graph_labels(ax, specs: List[dict], xmin, xmax, ymin, ymax, span: float):
    label_bbox = dict(
        boxstyle="round,pad=0.18,rounding_size=0.6",
        facecolor="white",
        edgecolor="#e8e4df",
        linewidth=0.45,
        alpha=0.92,
    )
    char_w = span * 0.019
    line_h = span * 0.034
    for spec in specs:
        lx, ly = spec["x"], spec["y"]
        lab = spec["text"]
        ha = spec.get("ha", "left")
        va = spec.get("va", "center")
        lines = lab.split("\n")
        max_chars = max((len(line) for line in lines), default=0)
        tw = char_w * max(max_chars, 2)
        th = line_h * max(len(lines), 1)
        if ha == "left":
            x0, x1 = lx, lx + tw
        elif ha == "right":
            x0, x1 = lx - tw, lx
        else:
            x0, x1 = lx - tw / 2, lx + tw / 2
        if va == "bottom":
            y0, y1 = ly, ly + th
        elif va == "top":
            y0, y1 = ly - th, ly
        else:
            y0, y1 = ly - th / 2, ly + th / 2
        xmin = min(xmin, x0)
        xmax = max(xmax, x1)
        ymin = min(ymin, y0)
        ymax = max(ymax, y1)
        ax.text(
            lx,
            ly,
            lab,
            fontsize=spec["fontsize"],
            color=spec["color"],
            fontweight=spec["weight"],
            ha=ha,
            va=va,
            zorder=40,
            clip_on=False,
            bbox=label_bbox,
        )
    return xmin, xmax, ymin, ymax


def _add_graph_legend(ax) -> None:
    try:
        from matplotlib.lines import Line2D

        handles = [
            Line2D(
                [0],
                [0],
                marker="D",
                color="w",
                markerfacecolor=_GRAPH_MARKER_FILL,
                markeredgecolor=_GRAPH_MARKER_EDGE,
                markersize=10,
                label="Marker gene",
            ),
            Line2D(
                [0],
                [0],
                marker="s",
                color="w",
                markerfacecolor=_GRAPH_CELL_FILL,
                markeredgecolor=_GRAPH_CELL_EDGE,
                markersize=8,
                label="Cell type",
            ),
        ]
        ax.legend(
            handles=handles,
            loc="upper right",
            frameon=True,
            framealpha=0.95,
            edgecolor="#e2e8f0",
            fontsize=8,
            borderpad=0.6,
        )
    except Exception:
        pass


def _marker_cell_graph_matplotlib_manual_export(graph, fmt: str) -> Optional[bytes]:
    """Draw edges + node patches + labels purely in matplotlib (no igraph.plot)."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle, Rectangle, RegularPolygon
    except ImportError:
        global _mpl_import_logged
        if not _mpl_import_logged:
            _graph_stderr(
                "matplotlib is not installed or failed to import. "
                "Install it for marker graphs with visible node shapes: "
                "pip install matplotlib   (or: pip install -r requirements.txt). "
                "Until then, graphs use igraph SVG only (limited node styling)."
            )
            _mpl_import_logged = True
        return None

    try:
        if graph.vcount() == 0:
            return None
    except Exception:
        return None

    layout = _contract_layout_toward_centroid(
        _layout_for_report(graph), _layout_contract_factor(graph.vcount())
    )
    n = graph.vcount()
    xs = [float(layout[i][0]) for i in range(n)]
    ys = [float(layout[i][1]) for i in range(n)]
    raw_labels = _vertex_labels(graph)
    labels = raw_labels if raw_labels is not None else [""] * n
    vcolors, vshapes = _prepare_graph_visual_attrs(graph)
    pk = _marker_cell_plot_kw(graph)

    span = max(max(xs) - min(xs), max(ys) - min(ys), 1e-9)
    diamond_rot = math.pi / 4.0
    label_specs = _plan_graph_labels(graph, layout, labels, span, pk)

    try:
        fig, ax = plt.subplots(figsize=_marker_graph_figsize(n), dpi=_GRAPH_DPI)
        fig.patch.set_facecolor("#ffffff")
        ax.set_facecolor(_GRAPH_BG)

        for s, t in graph.get_edgelist():
            ax.plot(
                [xs[s], xs[t]],
                [ys[s], ys[t]],
                color=_GRAPH_EDGE_COLOR,
                linewidth=0.9,
                alpha=0.5,
                solid_capstyle="round",
                zorder=1,
                clip_on=False,
            )

        for i in range(n):
            x, y = xs[i], ys[i]
            vt = _vertex_type_at(graph, i)
            sh = vshapes[i]
            if vt == "Marker":
                rad = span * 0.028
                ec, lw = _GRAPH_MARKER_EDGE, 1.35
            elif vt == "Cell":
                rad = span * 0.014
                ec, lw = _GRAPH_CELL_EDGE, 0.9
            else:
                rad = span * 0.02
                ec, lw = "#475569", 0.8
            fc = vcolors[i]
            if sh == "diamond":
                shadow = RegularPolygon(
                    (x + span * 0.003, y - span * 0.003),
                    4,
                    radius=rad * 1.02,
                    orientation=diamond_rot,
                    facecolor=(0, 0, 0, 0.08),
                    edgecolor="none",
                    zorder=10,
                    clip_on=False,
                )
                ax.add_patch(shadow)
                pat = RegularPolygon(
                    (x, y),
                    4,
                    radius=rad,
                    orientation=diamond_rot,
                    facecolor=fc,
                    edgecolor=ec,
                    linewidth=lw,
                    zorder=12,
                    clip_on=False,
                )
            elif sh == "rect":
                pat = Rectangle(
                    (x - rad, y - rad),
                    2 * rad,
                    2 * rad,
                    facecolor=fc,
                    edgecolor=ec,
                    linewidth=lw,
                    zorder=12,
                    clip_on=False,
                )
            else:
                pat = Circle(
                    (x, y),
                    rad,
                    facecolor=fc,
                    edgecolor=ec,
                    linewidth=lw,
                    zorder=12,
                    clip_on=False,
                )
            ax.add_patch(pat)

        xmin, xmax = min(xs), max(xs)
        ymin, ymax = min(ys), max(ys)
        node_margin = span * 0.05
        xmin -= node_margin
        xmax += node_margin
        ymin -= node_margin
        ymax += node_margin

        _try_adjusttext_positions(ax, label_specs)
        xmin, xmax, ymin, ymax = _draw_graph_labels(
            ax, label_specs, xmin, xmax, ymin, ymax, span
        )

        _add_graph_legend(ax)

        margin = span * 0.12
        ax.set_xlim(xmin - margin, xmax + margin)
        ax.set_ylim(ymin - margin, ymax + margin)
        ax.axis("off")
        fig.subplots_adjust(0, 0, 1, 1)

        buf = io.BytesIO()
        _matplotlib_savefig_graph_bytes(fig, ax, buf, fmt)
        plt.close(fig)
        out = buf.getvalue()
        if not out:
            return None
        return out
    except Exception as exc:
        try:
            import matplotlib.pyplot as plt

            plt.close("all")
        except Exception:
            pass
        warnings.warn(
            "CellTyper marker–cell graph matplotlib export failed (%s). "
            "Using igraph/Cairo/SVG fallback if available. Detail: %s"
            % (exc.__class__.__name__, exc),
            UserWarning,
            stacklevel=2,
        )
        _graph_stderr(
            "Marker graph: matplotlib patch renderer failed (%s): %s"
            % (exc.__class__.__name__, exc)
        )
        return None


def _marker_cell_graph_matplotlib_igraph_plot_export(graph, fmt: str) -> Optional[bytes]:
    """igraph.plot on a matplotlib Axes; used when patch-only renderer did not run."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import igraph as ig
    except ImportError:
        return None

    labels = _vertex_labels(graph)
    n = graph.vcount()
    layout = _contract_layout_toward_centroid(
        _layout_for_report(graph), _layout_contract_factor(n)
    )
    vcolors, vshapes = _prepare_graph_visual_attrs(graph)
    pk = _marker_cell_plot_kw(graph)

    plot_labels = None
    if labels:
        plot_labels = []
        for i, lab in enumerate(labels):
            plot_labels.append(_label_display_text(graph, i, lab) or "")

    base: dict = {
        **pk,
        "vertex_color": vcolors,
        "vertex_shape": vshapes,
        "edge_width": 0.85,
        "edge_color": _GRAPH_EDGE_COLOR,
        "vertex_label_dist": 0.65,
    }

    hide_v = {"vertex_size": [0] * n, "vertex_frame_width": [0] * n}
    invisible = {"vertex_color": ["#00000000"] * n, "vertex_frame_width": [0] * n}
    attempts: List[dict] = []
    if plot_labels:
        attempts.append({**base, "vertex_label": plot_labels, **hide_v})
    attempts.append({**base, **hide_v})
    if plot_labels:
        attempts.append({**base, "vertex_label": plot_labels, **invisible})
    attempts.append({**base, **invisible})

    for extra in attempts:
        try:
            fig, ax = plt.subplots(figsize=_marker_graph_figsize(n), dpi=_GRAPH_DPI)
            ax.set_facecolor(_GRAPH_BG)
            kw = {"target": ax, "layout": layout, **extra}
            ig.plot(graph, **kw)
            _matplotlib_push_edges_under_vertices(ax)
            _matplotlib_overlay_vertex_markers(
                ax, graph, layout, vcolors, vshapes, pk["vertex_size"]
            )
            _matplotlib_include_outside_labels(ax, fig, pad_frac=0.58)
            for t in ax.texts:
                try:
                    t.set_zorder(45)
                except Exception:
                    pass
            buf = io.BytesIO()
            _matplotlib_savefig_tight_including_markers(ax, fig, buf, fmt)
            plt.close(fig)
            data = buf.getvalue()
            if data:
                return data
        except Exception:
            try:
                plt.close("all")
            except Exception:
                pass
    return None


def _marker_cell_graph_matplotlib_export_with_source(
    graph, fmt: str
) -> Tuple[Optional[bytes], str]:
    out_manual = _marker_cell_graph_matplotlib_manual_export(graph, fmt)
    if out_manual:
        return out_manual, "matplotlib-patches"
    ig_mpl = _marker_cell_graph_matplotlib_igraph_plot_export(graph, fmt)
    if ig_mpl:
        return ig_mpl, "matplotlib-igraph-plot"
    return None, ""


def _marker_cell_graph_matplotlib_export(graph, fmt: str) -> Optional[bytes]:
    """Prefer patch matplotlib; fall back to igraph.plot on matplotlib Axes."""
    data, _ = _marker_cell_graph_matplotlib_export_with_source(graph, fmt)
    return data


def _marker_cell_graph_to_png_matplotlib(graph) -> Optional[bytes]:
    return _marker_cell_graph_matplotlib_export(graph, "png")


def _marker_cell_graph_to_svg_matplotlib(graph) -> Optional[bytes]:
    return _marker_cell_graph_matplotlib_export(graph, "svg")


def _marker_cell_graph_cairo_export(graph, suffix: str) -> Optional[bytes]:
    try:
        import igraph as ig
    except ImportError:
        return None

    igplot = getattr(ig, "plot", None)
    if igplot is None:
        return None

    labels = _vertex_labels(graph)
    layout = _contract_layout_toward_centroid(_layout_for_report(graph), 0.8)
    vcolors, vshapes = _prepare_graph_visual_attrs(graph)
    pk = _marker_cell_plot_kw(graph)
    attempts: List[dict] = []
    full: dict = {
        "bbox": (1800, 1800),
        "layout": layout,
        "edge_width": 1.0,
        "margin": 220,
        "vertex_color": vcolors,
        "vertex_shape": vshapes,
        "vertex_size": pk["vertex_size"],
        "vertex_frame_width": pk["vertex_frame_width"],
        "vertex_frame_color": pk["vertex_frame_color"],
        "vertex_label_size": pk["vertex_label_size"],
        "vertex_label_color": pk["vertex_label_color"],
    }
    if labels:
        plot_labels = []
        for i, lab in enumerate(labels):
            plot_labels.append(_label_display_text(graph, i, lab) or "")
        full = {
            **full,
            "vertex_label": plot_labels,
        }
    attempts.append(full)
    mid = {k: v for k, v in full.items() if k != "vertex_shape"}
    attempts.append(mid)
    attempts.append({"layout": layout, "vertex_color": vcolors, "vertex_shape": vshapes})
    attempts.append({"layout": layout, "vertex_color": vcolors})
    attempts.append({"layout": layout})

    for plot_kw in attempts:
        tmp_path = None
        try:
            with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
                tmp_path = tmp.name
            igplot(graph, target=tmp_path, **plot_kw)
            with open(tmp_path, "rb") as f:
                data = f.read()
                if data:
                    return data
        except Exception:
            continue
        finally:
            if tmp_path:
                try:
                    os.unlink(tmp_path)
                except OSError:
                    pass
    return None


def _marker_cell_graph_to_png_cairo(graph) -> Optional[bytes]:
    return _marker_cell_graph_cairo_export(graph, ".png")


def _marker_cell_graph_to_svg_cairo(graph) -> Optional[bytes]:
    return _marker_cell_graph_cairo_export(graph, ".svg")


def _marker_cell_graph_to_svg_bytes(graph) -> Optional[bytes]:
    try:
        import igraph as ig  # noqa: F401

        g = graph.copy()
        colors, shapes = _prepare_graph_visual_attrs(g)
        g.vs["color"] = colors
        g.vs["shape"] = shapes
        layout_svg = _contract_layout_toward_centroid(_layout_for_report(graph), 0.8)
        buf = io.StringIO()
        kw: dict = {
            "layout": layout_svg,
            "width": 1800,
            "height": 1800,
        }
        if "name" in g.vs.attributes():
            kw["labels"] = "name"
        kw["colors"] = "color"
        kw["shapes"] = "shape"
        kw["font_size"] = 14
        g.write_svg(buf, **kw)
        return buf.getvalue().encode("utf-8")
    except Exception:
        try:
            buf = io.StringIO()
            layout_fb = _contract_layout_toward_centroid(_layout_for_report(graph), 0.8)
            kw2: dict = {"width": 1800, "height": 1800, "layout": layout_fb}
            if "name" in graph.vs.attributes():
                kw2["labels"] = "name"
            graph.write_svg(buf, **kw2)
            return buf.getvalue().encode("utf-8")
        except Exception:
            return None


def marker_cell_graph_image_payload(graph) -> Optional[tuple[bytes, str, str]]:
    """Return (image bytes, MIME type, backend label) for HTML embedding and debugging."""
    if graph is None:
        return None
    try:
        if graph.vcount() == 0:
            return None
    except Exception:
        return None

    png, png_src = _marker_cell_graph_matplotlib_export_with_source(graph, "png")
    if png:
        _graph_stderr("Marker graph: PNG via %s." % png_src)
        return png, "image/png", png_src

    png_cairo = _marker_cell_graph_to_png_cairo(graph)
    if png_cairo:
        _graph_stderr("Marker graph: PNG via igraph Cairo (not matplotlib patch renderer).")
        return png_cairo, "image/png", "igraph-cairo"

    svg, svg_src = _marker_cell_graph_matplotlib_export_with_source(graph, "svg")
    if svg:
        svg = _svg_bytes_expand_viewbox(svg, pad_ratio=0.30)
        _graph_stderr("Marker graph: SVG via %s." % svg_src)
        return svg, "image/svg+xml", svg_src

    svg_cairo = _marker_cell_graph_to_svg_cairo(graph)
    if svg_cairo:
        _graph_stderr("Marker graph: SVG via igraph Cairo.")
        return (
            _svg_bytes_expand_viewbox(svg_cairo, 0.30),
            "image/svg+xml",
            "igraph-cairo-svg",
        )

    svg = _marker_cell_graph_to_svg_bytes(graph)
    if svg:
        _graph_stderr("Marker graph: SVG via igraph Graph.write_svg (limited styling).")
        return svg, "image/svg+xml", "igraph-write_svg"

    _graph_stderr("Marker graph: could not generate any image (check matplotlib and igraph).")
    return None


def marker_cell_graph_to_png_bytes(graph) -> Optional[bytes]:
    """Return PNG bytes only (no SVG). Prefer :func:`marker_cell_graph_image_payload` for HTML."""
    payload = marker_cell_graph_image_payload(graph)
    if payload and payload[1] == "image/png":
        return payload[0]
    return None


_PIE_COLORS = (
    "#5a7d8c",
    "#b8a48c",
    "#8b9a7d",
    "#7d6b8a",
    "#a68b7d",
    "#6b8a9a",
    "#9a8b6b",
    "#8a6b7d",
    "#7a8a6b",
    "#968878",
)


def _consistency_pie_chart_b64(
    counts: dict[str, Any],
    scores: dict[str, Any],
    *,
    n_samples: int = 0,
) -> Optional[str]:
    """Render self-consistency vote distribution as a pie chart PNG (base64)."""
    if not counts:
        return None
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    ordered = sorted(counts.keys(), key=lambda k: (-int(counts[k]), str(k)))
    sizes = [int(counts[k]) for k in ordered]
    total = sum(sizes) or max(n_samples, 1)

    legend_labels: List[str] = []
    for label in ordered:
        c = int(counts[label])
        frac = float(scores.get(label, c / total))
        legend_labels.append(f"{label} ({c}/{total}, {frac:.0%})")

    colors = [_PIE_COLORS[i % len(_PIE_COLORS)] for i in range(len(ordered))]
    max_label_len = max((len(str(k)) for k in ordered), default=20)
    fig_w = min(10.0, 5.2 + 0.08 * max_label_len)

    fig, ax = plt.subplots(figsize=(fig_w, 4.5), dpi=120)
    wedges, _texts, autotexts = ax.pie(
        sizes,
        colors=colors,
        startangle=90,
        autopct=lambda pct: f"{pct:.0f}%" if pct >= 8 else "",
        pctdistance=0.72,
        wedgeprops={"linewidth": 0.8, "edgecolor": "white"},
    )
    for t in autotexts:
        t.set_fontsize(9)
        t.set_color("white")
        t.set_fontweight("bold")

    ax.legend(
        wedges,
        legend_labels,
        title="Cell type",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=8.5,
        frameon=False,
    )
    ax.set_aspect("equal")
    fig.subplots_adjust(right=0.52 if len(ordered) > 2 else 0.62)

    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", facecolor="white", dpi=120)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def _cluster_graph_html_section(cluster_id_esc: str, b64: Optional[str], mime_hint: Optional[str]) -> str:
    """Embed PNG via img data-URI; embed SVG inline (some browsers block SVG in img for local files)."""
    if not b64:
        return ""
    try:
        raw = base64.b64decode(b64, validate=False)
    except Exception:
        return ""
    if not raw:
        return ""

    parts: List[str] = []
    is_png = raw.startswith(b"\x89PNG\r\n\x1a\n") or raw.startswith(b"\x89PNG\n")

    if is_png:
        parts.append(
            "<p><img class='graph' "
            f"src='data:image/png;base64,{b64}' "
            f"alt='marker_cell_graph cluster {cluster_id_esc}'></p>"
        )
        return "\n".join(parts)

    try:
        svg = raw.decode("utf-8")
    except Exception:
        svg = None

    if svg and "<svg" in svg.lower():
        parts.append(f"<div class='graph-svg-wrap'>\n{svg}\n</div>")
        return "\n".join(parts)

    m = mime_hint or "image/png"
    parts.append(
        "<p><img class='graph' "
        f"src='data:{html.escape(m)};base64,{b64}' "
        f"alt='marker_cell_graph cluster {cluster_id_esc}'></p>"
    )
    return "\n".join(parts)


_GENE_CONTEXT_SPECS: List[Tuple[str, str, str, Tuple[str, ...]]] = [
    ("Description & summary", "gene_desc_summary_context", "Description", ("Gene's Description & Summary",)),
    ("Function", "gene_func_summary_context", "Function", ("Gene's function",)),
    ("Cell types / tissue", "gene_cell_summary_context", "Related cell type", ("Cell types or Tissue related to Gene",)),
    ("Pathway", "gene_pathway_context", "Pathway", ("Pathway related to Gene",)),
    ("Protein location", "gene_protein_loc_context", "Protein location", ("Gene's Protein Location",)),
    ("Protein–protein interaction", "ppi_context", "PP interaction", ("Gene's Protein-Protein Interaction",)),
]


def _slice_context_for_marker(blob: str, marker: str, section_labels: Tuple[str, ...]) -> str:
    if not blob or not str(marker).strip():
        return ""
    marker_esc = re.escape(str(marker).strip())
    for label in section_labels:
        label_esc = re.escape(label)
        pat = rf"\[{label_esc}\]\s*\n{marker_esc}:\s*\n(.*?)(?=\n\[[^\]]+\]\s*\n|\Z)"
        m = re.search(pat, blob, re.DOTALL)
        if m:
            return m.group(1).strip()
    return ""


def build_per_marker_gene_contexts(
    markers: List[str],
    *,
    species: str,
    gene_desc_summary_context: str = "",
    gene_func_summary_context: str = "",
    gene_cell_summary_context: str = "",
    gene_pathway_context: str = "",
    gene_protein_loc_context: str = "",
    ppi_context: str = "",
) -> List[dict[str, Any]]:
    """Assemble per-marker gene context sections for HTML report (geneDB or parsed retrieval text)."""
    from geneRetriever import get_gene_db

    blobs = {
        "gene_desc_summary_context": gene_desc_summary_context or "",
        "gene_func_summary_context": gene_func_summary_context or "",
        "gene_cell_summary_context": gene_cell_summary_context or "",
        "gene_pathway_context": gene_pathway_context or "",
        "gene_protein_loc_context": gene_protein_loc_context or "",
        "ppi_context": ppi_context or "",
    }

    try:
        gene_db = get_gene_db(species)
    except Exception:
        gene_db = None

    out: List[dict[str, Any]] = []
    for marker in markers:
        marker_s = str(marker).strip()
        if not marker_s:
            continue

        hit = None
        if gene_db is not None:
            rows = gene_db.loc[gene_db["Gene"].astype(str).str.lower() == marker_s.lower()]
            if len(rows) > 0:
                hit = rows.iloc[0]

        sections: List[dict[str, str]] = []
        for title, blob_key, db_col, parse_labels in _GENE_CONTEXT_SPECS:
            text = ""
            if hit is not None:
                val = hit.get(db_col)
                if val is not None and str(val).strip() and str(val).lower() != "nan":
                    text = str(val).strip()
            else:
                text = _slice_context_for_marker(blobs[blob_key], marker_s, parse_labels)
            if text:
                sections.append({"title": title, "body": text})

        out.append({"gene": marker_s, "sections": sections})

    return out


def _report_styles() -> str:
    return """
:root {
  --bg: #f4f6fb;
  --surface: #ffffff;
  --surface-2: #f8fafc;
  --text: #0f172a;
  --muted: #64748b;
  --border: #e2e8f0;
  --accent: #2563eb;
  --accent-2: #7c3aed;
  --accent-soft: #eff6ff;
  --success: #059669;
  --shadow: 0 10px 40px rgba(15, 23, 42, 0.08);
  --radius: 14px;
  --mono: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  font-family: Inter, ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
  line-height: 1.55;
  color: var(--text);
  background: linear-gradient(180deg, #eef2ff 0%, var(--bg) 220px, var(--bg) 100%);
}
.report-wrap { max-width: 1080px; margin: 0 auto; padding: 2rem 1.25rem 3rem; }
.report-hero {
  background: linear-gradient(135deg, #1d4ed8 0%, #6d28d9 100%);
  color: #fff;
  border-radius: calc(var(--radius) + 4px);
  padding: 1.75rem 1.5rem;
  box-shadow: var(--shadow);
  margin-bottom: 1.75rem;
}
.report-hero h1 { margin: 0 0 .35rem; font-size: 1.85rem; letter-spacing: -0.02em; }
.report-hero .subtitle { margin: 0; opacity: .92; font-size: .98rem; }
.meta-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
  gap: .75rem;
  margin-top: 1.25rem;
}
.meta-card {
  background: rgba(255,255,255,.12);
  border: 1px solid rgba(255,255,255,.18);
  border-radius: 10px;
  padding: .65rem .8rem;
}
.meta-card .label { display: block; font-size: .72rem; text-transform: uppercase; letter-spacing: .06em; opacity: .85; }
.meta-card .value { display: block; font-weight: 600; margin-top: .15rem; word-break: break-word; }
.cluster-card {
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: var(--radius);
  box-shadow: var(--shadow);
  padding: 1.35rem 1.35rem 1.1rem;
  margin-bottom: 1.5rem;
}
.cluster-head { display: flex; flex-wrap: wrap; align-items: center; gap: .75rem 1rem; margin-bottom: 1rem; }
.cluster-badge {
  display: inline-flex;
  align-items: center;
  padding: .35rem .7rem;
  border-radius: 999px;
  background: var(--accent-soft);
  color: var(--accent);
  font-weight: 700;
  font-size: .82rem;
  letter-spacing: .03em;
  text-transform: uppercase;
}
.prediction {
  margin: 0;
  font-size: 1.45rem;
  letter-spacing: -0.02em;
  flex: 1 1 auto;
}
.prediction span { color: var(--success); }
.section-title {
  margin: 1.35rem 0 .65rem;
  font-size: 1rem;
  font-weight: 700;
  color: var(--text);
  display: flex;
  align-items: center;
  gap: .5rem;
}
.section-title::before {
  content: "";
  width: 4px;
  height: 1.1em;
  border-radius: 4px;
  background: linear-gradient(180deg, var(--accent), var(--accent-2));
}
.marker-pills { display: flex; flex-wrap: wrap; gap: .45rem; }
.marker-pill {
  display: inline-flex;
  padding: .28rem .62rem;
  border-radius: 999px;
  background: var(--surface-2);
  border: 1px solid var(--border);
  font-size: .84rem;
  font-family: var(--mono);
  color: #334155;
}
.gene-panel {
  border: 1px solid var(--border);
  border-radius: 12px;
  background: var(--surface-2);
  margin: .55rem 0;
  overflow: hidden;
}
.gene-panel > summary {
  cursor: pointer;
  list-style: none;
  padding: .75rem 1rem;
  font-weight: 700;
  font-family: var(--mono);
  background: #fff;
  border-bottom: 1px solid transparent;
}
.gene-panel[open] > summary { border-bottom-color: var(--border); }
.gene-panel > summary::-webkit-details-marker { display: none; }
.gene-panel > summary::after {
  content: "▾";
  float: right;
  color: var(--muted);
  transition: transform .15s ease;
}
.gene-panel[open] > summary::after { transform: rotate(180deg); }
.gene-panel-body { padding: .25rem .85rem .85rem; }
.context-block { margin: .75rem 0 0; }
.context-block h4 {
  margin: 0 0 .35rem;
  font-size: .78rem;
  text-transform: uppercase;
  letter-spacing: .05em;
  color: var(--muted);
}
.context-pre {
  margin: 0;
  white-space: pre-wrap;
  word-break: break-word;
  font-family: var(--mono);
  font-size: .78rem;
  line-height: 1.5;
  background: #fff;
  border: 1px solid var(--border);
  border-radius: 10px;
  padding: .75rem .85rem;
  max-height: 420px;
  overflow: auto;
}
.data-table { width: 100%; border-collapse: collapse; margin: .5rem 0 0; font-size: .92rem; }
.data-table th, .data-table td { border: 1px solid var(--border); padding: .5rem .65rem; text-align: left; }
.data-table th { background: var(--surface-2); font-size: .78rem; text-transform: uppercase; letter-spacing: .04em; color: var(--muted); }
.data-table tr:nth-child(even) td { background: #fbfdff; }
.consistency-chart {
  margin: .5rem 0 0;
  padding: .85rem 1rem;
  background: var(--surface-2);
  border: 1px solid var(--border);
  border-radius: 12px;
  display: inline-block;
  max-width: 100%;
}
.consistency-chart img {
  max-width: 100%;
  height: auto;
  display: block;
}
.llm-panel {
  border: 1px solid var(--border);
  border-radius: 10px;
  background: #fff;
  margin: .5rem 0;
}
.llm-panel > summary { cursor: pointer; padding: .65rem .85rem; font-weight: 600; }
.llm-panel > summary::-webkit-details-marker { display: none; }
.llm-panel pre {
  margin: 0;
  border-top: 1px solid var(--border);
  padding: .75rem .85rem;
  white-space: pre-wrap;
  word-break: break-word;
  font-size: .82rem;
  line-height: 1.5;
  max-height: 360px;
  overflow: auto;
}
.graph-card {
  border: 1px solid var(--border);
  border-radius: 14px;
  padding: 1rem;
  background: linear-gradient(180deg, #f8fafc 0%, #ffffff 100%);
  box-shadow: inset 0 1px 0 rgba(255,255,255,.8);
}
img.graph {
  max-width: 100%;
  height: auto;
  border-radius: 10px;
  display: block;
  background: #f8fafc;
  border: 1px solid #e2e8f0;
}
.graph-svg-wrap { max-width: 100%; overflow: auto; border-radius: 8px; padding: 8px; background: #fff; }
.graph-svg-wrap svg { max-width: 100%; height: auto; display: block; overflow: visible; }
.graph-legend { font-size: .88rem; color: var(--muted); margin: .5rem 0 0; }
.lg-gene { color: #5a7d8c; font-weight: 700; }
.lg-cell { color: #8a7560; font-weight: 700; }
.note { font-size: .85rem; color: var(--muted); margin-top: .5rem; }
.empty-note { color: var(--muted); font-size: .9rem; font-style: italic; }
"""


def _gene_context_html(marker_contexts: List[dict[str, Any]]) -> str:
    if not marker_contexts:
        return "<p class='empty-note'>No per-gene context available.</p>"

    parts: List[str] = ["<div class='gene-context-list'>"]
    for item in marker_contexts:
        gene = html.escape(str(item.get("gene", "")))
        sections = item.get("sections") or []
        parts.append(f"<details class='gene-panel'><summary>{gene}</summary><div class='gene-panel-body'>")
        if not sections:
            parts.append("<p class='empty-note'>No context retrieved for this marker.</p>")
        else:
            for sec in sections:
                title = html.escape(str(sec.get("title", "Context")))
                body = html.escape(str(sec.get("body", "")))
                parts.append(
                    f"<div class='context-block'><h4>{title}</h4>"
                    f"<pre class='context-pre'>{body}</pre></div>"
                )
        parts.append("</div></details>")
    parts.append("</div>")
    return "\n".join(parts)


def build_html_report(
    *,
    markers_path: str,
    species: Optional[str],
    tissue: str,
    condition: str,
    clusters: List[dict[str, Any]],
    title: str = "CellTyper report",
) -> str:
    """Return a single self-contained HTML document (marker graphs as PNG or SVG data URIs)."""
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    parts: List[str] = []
    parts.append("<!DOCTYPE html>")
    parts.append("<html lang='en'><head><meta charset='utf-8'>")
    parts.append("<meta name='viewport' content='width=device-width, initial-scale=1'>")
    parts.append(f"<title>{html.escape(title)}</title>")
    parts.append(f"<style>{_report_styles()}</style></head><body>")
    parts.append("<div class='report-wrap'>")
    parts.append("<header class='report-hero'>")
    parts.append(f"<h1>{html.escape(title)}</h1>")
    parts.append(f"<p class='subtitle'>Generated {html.escape(now)}</p>")
    parts.append("<div class='meta-grid'>")
    parts.append(
        f"<div class='meta-card'><span class='label'>Markers file</span>"
        f"<span class='value'>{html.escape(markers_path)}</span></div>"
    )
    parts.append(
        f"<div class='meta-card'><span class='label'>Species</span>"
        f"<span class='value'>{html.escape(str(species))}</span></div>"
    )
    parts.append(
        f"<div class='meta-card'><span class='label'>Tissue</span>"
        f"<span class='value'>{html.escape(tissue)}</span></div>"
    )
    parts.append(
        f"<div class='meta-card'><span class='label'>Condition</span>"
        f"<span class='value'>{html.escape(condition)}</span></div>"
    )
    parts.append(
        f"<div class='meta-card'><span class='label'>Clusters</span>"
        f"<span class='value'>{len(clusters)}</span></div>"
    )
    parts.append("</div></header><main>")

    for block in clusters:
        cid = html.escape(str(block.get("cluster_id", "")))
        pred = block.get("predicted_display") or block.get("predicted") or ""
        pred_h = html.escape(str(pred))
        parts.append("<article class='cluster-card'>")
        parts.append("<div class='cluster-head'>")
        parts.append(f"<span class='cluster-badge'>Cluster {cid}</span>")
        parts.append(f"<h2 class='prediction'>Predicted: <span>{pred_h}</span></h2>")
        parts.append("</div>")

        mg = block.get("marker_genes") or []
        if mg:
            parts.append("<h3 class='section-title'>Input marker genes</h3>")
            parts.append("<div class='marker-pills'>")
            for g in mg:
                parts.append(f"<span class='marker-pill'>{html.escape(str(g))}</span>")
            parts.append("</div>")

        marker_contexts = block.get("marker_gene_contexts")
        if marker_contexts is None and mg:
            marker_contexts = build_per_marker_gene_contexts(
                mg,
                species=str(species or ""),
                gene_desc_summary_context=block.get("gene_desc_summary_context") or "",
                gene_func_summary_context=block.get("gene_func_summary_context") or "",
                gene_cell_summary_context=block.get("gene_cell_summary_context") or "",
                gene_pathway_context=block.get("gene_pathway_context") or "",
                gene_protein_loc_context=block.get("gene_protein_loc_context") or "",
                ppi_context=block.get("ppi_context") or "",
            )
        parts.append("<h3 class='section-title'>Gene context (full text per marker)</h3>")
        parts.append(_gene_context_html(marker_contexts or []))

        parsed = block.get("llm_parsed_celltypes") or []
        if parsed:
            parts.append("<h3 class='section-title'>LLM sampled labels</h3>")
            parts.append("<p>")
            parts.append(
                ", ".join(html.escape(str(x)) if x else "(parse failed)" for x in parsed)
            )
            parts.append("</p>")

        counts = block.get("llm_consistency_counts") or {}
        scores = block.get("llm_consistency_scores") or {}
        if counts:
            parts.append("<h3 class='section-title'>Self-consistency</h3>")
            pie_b64 = _consistency_pie_chart_b64(
                counts, scores, n_samples=len(parsed) if parsed else sum(counts.values())
            )
            if pie_b64:
                parts.append(
                    "<div class='consistency-chart'>"
                    f"<img src='data:image/png;base64,{pie_b64}' "
                    "alt='Self-consistency pie chart'></div>"
                )
            else:
                parts.append(
                    "<table class='data-table'><thead><tr>"
                    "<th>Cell type</th><th>Count</th><th>Fraction</th>"
                    "</tr></thead><tbody>"
                )
                for label in sorted(counts.keys(), key=lambda k: (-counts[k], k)):
                    c = counts[label]
                    frac = scores.get(label, c / max(len(parsed), 1) if parsed else 0)
                    parts.append(
                        "<tr><td>{0}</td><td>{1}</td><td>{2:.3f}</td></tr>".format(
                            html.escape(str(label)), c, float(frac)
                        )
                    )
                parts.append("</tbody></table>")

        samples = block.get("llm_sample_texts") or []
        if samples:
            parts.append("<h3 class='section-title'>LLM responses</h3>")
            for i, txt in enumerate(samples, start=1):
                parts.append(f"<details class='llm-panel'><summary>Sample {i}</summary>")
                parts.append(f"<pre>{html.escape(str(txt) if txt is not None else '')}</pre></details>")

        sec = _cluster_graph_html_section(cid, block.get("graph_png_b64"), block.get("graph_image_mime"))
        if sec:
            parts.append("<h3 class='section-title'>Marker–cell graph</h3>")
            parts.append("<div class='graph-card'>")
            parts.append(sec)
            gbe = block.get("graph_export_backend")
            if gbe and not str(gbe).startswith("matplotlib-"):
                parts.append(
                    "<p class='note'><strong>Note:</strong> Matplotlib is not being used for this graph. "
                    "Install it in the same Python you use for CellTyper so marker nodes render as filled ◆/■: "
                    "<code>pip install matplotlib</code></p>"
                )
            parts.append(
                "<p class='graph-legend'>"
                "<span class='lg-gene'>◆</span> Marker gene &nbsp; "
                "<span class='lg-cell'>■</span> Cell type"
                "</p>"
            )
            parts.append("</div>")

        parts.append("</article>")

    parts.append("</main></div></body></html>")
    return "\n".join(parts)


def write_html_report(path: str, html_content: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(html_content)
