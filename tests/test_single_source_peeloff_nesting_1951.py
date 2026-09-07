"""A single-source fan-out bundle nests concentrically into a shared entry port.

Regression for #1951: three lines fan out from one junction, ride a shared
bypass trunk right and down, and peel off into a common LEFT entry port through
a right->down->right approach.  The turn is a half-turn that transposes the
bundle, so approach-X must run *opposite* to trunk-depth order.  The staggering
pass placed the descents in trunk order instead, tripping both the bundle-order
and peel-off-concentric render guards.  The layout must seat the descents on the
slots ``peeloff_target_slots`` (the guards' oracle) earns.
"""

from __future__ import annotations

from nf_metro.layout.engine import compute_layout
from nf_metro.layout.routing import compute_station_offsets, route_edges
from nf_metro.layout.routing.invariants import (
    assert_render_curve_invariants,
    check_bundle_order_preserved,
    check_peeloff_concentric,
)
from nf_metro.parser.mermaid import parse_metro_mermaid

REPRO = """\
%%metro title: Repro
%%metro line: riboseq | Ribo-seq | #e6007e
%%metro line: rnaseq | Matched RNA-seq | #2db572
%%metro line: tiseq | TI-seq | #2b6cb0
%%metro line: annotation | Hybrid annotation | #94a3b8 | dashed

%%metro grid: preprocessing, alignment, novel_transcripts | 0,0
%%metro grid: orf_calling, psite, te, reporting | 0,1
%%metro center_ports: true
%%metro x_spacing: 70

graph LR
    subgraph preprocessing [Read pre-processing]
        equalise[Equalise read lengths]
    end
    subgraph alignment [Alignment]
        star[STAR]
        umi_dedup[UMI-tools dedup]
        salmon[Salmon]
    end
    subgraph novel_transcripts [Transcript discovery]
        hybrid_merge[Merge GTF]
    end
    subgraph orf_calling [ORF calling]
        star_hybrid[STAR hybrid]
    end
    subgraph psite [P-site]
        ribowaltz[riboWaltz]
        plastid_psite[plastid P-site]
        quantify_orf_psite[ORF P-sites]
    end
    subgraph te [Translational efficiency]
        te_prep_orf[ORF count matrix]
    end
    subgraph reporting [Reporting]
        multiqc_final[MultiQC]
    end

    equalise -->|riboseq,rnaseq,tiseq| star
    umi_dedup -->|riboseq| ribowaltz
    umi_dedup -->|riboseq| plastid_psite
    hybrid_merge -->|annotation| star_hybrid
    quantify_orf_psite -->|riboseq| te_prep_orf
    salmon -->|rnaseq| te_prep_orf
    salmon -->|riboseq,rnaseq,tiseq| multiqc_final
"""


def _route():
    graph = parse_metro_mermaid(REPRO)
    compute_layout(graph)
    offsets = compute_station_offsets(graph)
    routes = route_edges(graph, station_offsets=offsets)
    return graph, offsets, routes


def test_single_source_peeloff_nests_concentrically() -> None:
    _graph, _offsets, routes = _route()
    assert check_peeloff_concentric(_graph, routes) == []
    assert check_bundle_order_preserved(routes) == []


def test_single_source_peeloff_passes_render_curve_invariants() -> None:
    graph, offsets, routes = _route()
    assert_render_curve_invariants(graph, routes, offsets)


def test_single_source_peeloff_reverses_x_against_trunk_depth() -> None:
    """The shallowest-trunk line takes the port-nearest (inner) descent column.

    riboseq rides the shallowest trunk and rnaseq the middle one; the half-turn
    transposes X against trunk depth, so riboseq's descent must sit *inboard*
    (larger X) of rnaseq's, not outboard.
    """
    _graph, _offsets, routes = _route()
    descents = {
        r.line_id: r.points[-3][0]
        for r in routes
        if r.edge.source == "__junction_9"
        and r.edge.target.startswith("reporting__entry_left")
    }
    assert {"riboseq", "rnaseq", "tiseq"} <= set(descents)
    assert descents["riboseq"] > descents["rnaseq"] > descents["tiseq"]
