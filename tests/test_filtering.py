import json
import tempfile
from pathlib import Path
import importlib.util


def load_pipeline_module():
    # pipeline code lives under unified_pipeline_clean
    module_path = (
        Path(__file__).resolve().parent.parent
        / "unified_pipeline_clean"
        / "nif_downstream_code"
        / "pangenome_pipeline_consolidated2.py"
    )
    spec = importlib.util.spec_from_file_location("pipeline_module", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_filter_keeps_only_gcf_complete():
    mod = load_pipeline_module()

    # enforce filters
    mod.config.require_refseq_gcf = True
    mod.config.require_complete_genome = True

    with tempfile.TemporaryDirectory() as tmpdir:
        tdir = Path(tmpdir)
        table = tdir / "assembly_quality.tsv"
        table.write_text(
            "assembly_accession\tassembly_level\n"
            "GCF_000001\tComplete Genome\n"
            "GCA_000002\tComplete Genome\n"
            "GCF_000003\tScaffold\n"
        )

        kept = mod.load_and_filter_assemblies(table, tdir)

        assert kept == ["GCF_000001"]
        assert (tdir / "filtered_accessions.txt").exists()
        assert (tdir / "assembly_quality.filtered.tsv").exists()
        assert (tdir / "filter_summary.json").exists()

        summary = json.loads((tdir / "filter_summary.json").read_text())
        assert summary["kept"] == 1


def test_filter_raises_when_none_pass():
    mod = load_pipeline_module()

    mod.config.require_refseq_gcf = True
    mod.config.require_complete_genome = True

    with tempfile.TemporaryDirectory() as tmpdir:
        tdir = Path(tmpdir)
        table = tdir / "assembly_quality.tsv"
        table.write_text(
            "assembly_accession\tassembly_level\n"
            "GCA_000003\tScaffold\n"
        )

        try:
            mod.load_and_filter_assemblies(table, tdir)
        except RuntimeError as exc:
            assert "No assemblies passed" in str(exc)
        else:
            raise AssertionError("Expected RuntimeError when nothing passes filter")
