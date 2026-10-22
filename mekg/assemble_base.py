#!/usr/bin/env python3
"""Build a new self-contained MEKG v2 release with leakage-safe public support."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import shutil
import sqlite3
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


TABLE_IDS = {"nodes": "node_id", "edges": "edge_id", "evidence": "evidence_id", "observations": "observation_id"}
PUBLIC = "PUBLIC_INDEPENDENT_SUPPORT"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gz_rows(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def merge_table(baseline: Path, public_files: list[Path], output: Path, id_field: str) -> tuple[int, int]:
    sources = [baseline] + public_files
    header = None
    seen = set()
    duplicate_count = 0
    count = 0
    with gzip.open(output, "wt", encoding="utf-8", newline="") as out:
        writer = None
        for path in sources:
            with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if header is None:
                    header = reader.fieldnames
                    writer = csv.DictWriter(out, fieldnames=header, delimiter="\t", lineterminator="\n")
                    writer.writeheader()
                if reader.fieldnames != header:
                    raise RuntimeError(f"Header mismatch for {path}: {reader.fieldnames} != {header}")
                for row in reader:
                    key = row[id_field]
                    if key in seen:
                        duplicate_count += 1
                        continue
                    seen.add(key)
                    writer.writerow(row)
                    count += 1
    return count, duplicate_count


def classify(predicate: str) -> tuple[str, float]:
    mapping = {
        "HAS_CLASSICAL_VISIBLE_MUTANT_PHENOTYPE_EVIDENCE": ("maize_direct_genetics", 3.0),
        "PUBLIC_TF_MUTANT_ALTERS_TARGET_EXPRESSION": ("maize_in_vivo_perturbation", 2.4),
        "HAS_REPORTED_GWAS_ASSOCIATION_WITH": ("maize_trait_association", 1.8),
        "SUPPORTED_BY_ARABIDOPSIS_ORTHOLOG_PHENOTYPE": ("cross_species_genetics", 1.3),
        "HAS_FUNCTIONALLY_CHARACTERIZED_RICE_ORTHOLOG": ("cross_species_genetics", 1.3),
        "PUBLIC_ENDOSPERM_REGULON_TARGETS": ("endosperm_TF_regulation", 1.5),
        "HAS_PUBLIC_H3K27AC_PROMOTER_SUPPORT": ("endosperm_chromatin", 1.0),
        "HAS_PUBLIC_H3K4ME1_PROMOTER_SUPPORT": ("endosperm_chromatin", 0.8),
        "HAS_PUBLIC_H3K4ME3_PROMOTER_SUPPORT": ("endosperm_chromatin", 1.0),
        "PUBLIC_RNA_EXPRESSION_AT": ("endosperm_expression", 1.0),
        "HAS_PUBLIC_EARLY_SEED_PROMOTER_ATAC_AT": ("early_seed_ATAC", 0.8),
        "MEMBER_OF_KEGG_PATHWAY": ("pathway_ontology", 0.6),
        "MEMBER_OF_PLANT_REACTOME_PATHWAY": ("pathway_ontology", 0.6),
        "PARTICIPATES_IN_PLANT_REACTOME_REACTION": ("pathway_ontology", 0.6),
        "ANNOTATED_WITH_GO": ("pathway_ontology", 0.5),
        "HAS_PUBLIC_FUNCTION_DESCRIPTION": ("database_function", 0.5),
        "HAS_UNIFORMMU_INSERTION": ("mutant_resource_availability", 0.3),
    }
    if predicate in mapping:
        return mapping[predicate]
    if predicate == "HIGH_CONFIDENCE_ORTHOLOG_OF":
        return "orthology", 0.4
    return "other_public_support", 0.25


def build_public_summary(targets: set[str], edge_files: list[Path], output: Path) -> tuple[list[dict], Counter]:
    records: dict[str, dict] = {}
    totals = Counter()
    for path in edge_files:
        for row in gz_rows(path):
            if row.get("independence_class") != PUBLIC or row.get("used_for_training") != "0":
                raise RuntimeError(f"Leakage flag violation in {path}: {row}")
            endpoints = []
            for role, node in (("outgoing", row["subject_node_id"]), ("incoming", row["object_node_id"])):
                if node.startswith("gene:"):
                    gene = node.removeprefix("gene:")
                    if gene in targets:
                        endpoints.append((gene, role))
            for gene, role in endpoints:
                rec = records.setdefault(gene, {"predicates": Counter(), "categories": Counter(),
                    "outgoing": 0, "incoming": 0, "edges": set(), "weights": {}})
                rec[role] += 1
                rec["edges"].add(row["edge_id"])
                predicate = row["predicate"]
                rec["predicates"][predicate] += 1
                category, weight = classify(predicate)
                rec["categories"][category] += 1
                rec["weights"][category] = max(rec["weights"].get(category, 0.0), weight)
                totals[predicate] += 1
    rows = []
    fields = ["gene_id", "public_evidence_edge_count", "public_evidence_predicate_count",
              "public_evidence_category_count", "public_support_tier", "public_support_score_explanatory",
              "outgoing_edge_count", "incoming_edge_count", "predicate_counts_json",
              "category_counts_json", "independence_class", "used_for_training",
              "score_semantics"]
    for gene in sorted(targets):
        rec = records.get(gene, {"predicates": Counter(), "categories": Counter(), "outgoing": 0,
                                 "incoming": 0, "edges": set(), "weights": {}})
        predicates = rec["predicates"]
        if predicates.get("HAS_CLASSICAL_VISIBLE_MUTANT_PHENOTYPE_EVIDENCE", 0):
            tier = "A_DIRECT_MAIZE_GENETICS"
        elif predicates.get("PUBLIC_TF_MUTANT_ALTERS_TARGET_EXPRESSION", 0) or predicates.get("HAS_REPORTED_GWAS_ASSOCIATION_WITH", 0):
            tier = "B_IN_VIVO_OR_ASSOCIATION"
        elif rec["edges"]:
            tier = "C_FUNCTIONAL_OR_MOLECULAR_SUPPORT"
        else:
            tier = "D_NO_PUBLIC_SUPPORT"
        raw_score = sum(rec["weights"].values())
        score = round(min(1.0, raw_score / 8.0), 6)
        rows.append({
            "gene_id": gene, "public_evidence_edge_count": len(rec["edges"]),
            "public_evidence_predicate_count": len(predicates),
            "public_evidence_category_count": len(rec["categories"]), "public_support_tier": tier,
            "public_support_score_explanatory": score, "outgoing_edge_count": rec["outgoing"],
            "incoming_edge_count": rec["incoming"],
            "predicate_counts_json": json.dumps(dict(sorted(predicates.items())), separators=(",", ":")),
            "category_counts_json": json.dumps(dict(sorted(rec["categories"].items())), separators=(",", ":")),
            "independence_class": PUBLIC, "used_for_training": 0,
            "score_semantics": "explanatory_breadth_only_not_model_probability_not_used_for_ranking",
        })
    with gzip.open(output, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    return rows, totals


def db_columns(connection: sqlite3.Connection, table: str) -> list[str]:
    return [row[1] for row in connection.execute(f"PRAGMA table_info({table})")]


def insert_gz(connection: sqlite3.Connection, table: str, path: Path) -> int:
    columns = db_columns(connection, table)
    if not columns:
        return 0
    count = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        use = [column for column in columns if column in (reader.fieldnames or [])]
        placeholders = ",".join("?" for _ in use)
        sql = f"INSERT OR IGNORE INTO {table} ({','.join(use)}) VALUES ({placeholders})"
        batch = []
        for row in reader:
            batch.append(tuple(row[column] for column in use))
            if len(batch) >= 5000:
                before = connection.total_changes
                connection.executemany(sql, batch)
                count += connection.total_changes - before
                batch.clear()
        if batch:
            before = connection.total_changes
            connection.executemany(sql, batch)
            count += connection.total_changes - before
    return count


def build_sqlite(baseline_db: Path, output_db: Path, public: dict[str, list[Path]], aliases: list[Path], summary_rows: list[dict]) -> dict:
    source = sqlite3.connect(f"file:{baseline_db}?mode=ro", uri=True)
    destination = sqlite3.connect(output_db)
    source.backup(destination)
    source.close()
    destination.execute("PRAGMA journal_mode=DELETE")
    inserted = Counter()
    with destination:
        for table in ("nodes", "edges", "evidence", "observations"):
            for path in public.get(table, []):
                inserted[table] += insert_gz(destination, table, path)
        # Keep public aliases in a dedicated table because the baseline alias
        # table uses a different schema. This also makes provenance isolation
        # explicit and avoids altering the original alias semantics.
        destination.execute("DROP TABLE IF EXISTS public_aliases")
        destination.execute("""CREATE TABLE public_aliases (
            gene_id TEXT, alias_type TEXT, alias TEXT, source_node_id TEXT,
            mapping_status TEXT, used_for_training INTEGER, independence_class TEXT)""")
        alias_sql = "INSERT INTO public_aliases VALUES (?,?,?,?,?,?,?)"
        for path in aliases:
            with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                batch = []
                for row in reader:
                    batch.append(tuple(row.get(field, "") for field in (
                        "gene_id", "alias_type", "alias", "source_node_id", "mapping_status",
                        "used_for_training", "independence_class")))
                    if len(batch) >= 5000:
                        destination.executemany(alias_sql, batch); inserted["public_aliases"] += len(batch); batch.clear()
                if batch:
                    destination.executemany(alias_sql, batch); inserted["public_aliases"] += len(batch)
        destination.execute("DROP TABLE IF EXISTS public_gene_support")
        destination.execute("""CREATE TABLE public_gene_support (
            gene_id TEXT PRIMARY KEY, public_evidence_edge_count INTEGER,
            public_evidence_predicate_count INTEGER, public_evidence_category_count INTEGER,
            public_support_tier TEXT, public_support_score_explanatory REAL,
            outgoing_edge_count INTEGER, incoming_edge_count INTEGER,
            predicate_counts_json TEXT, category_counts_json TEXT,
            independence_class TEXT, used_for_training INTEGER, score_semantics TEXT)""")
        fields = list(summary_rows[0])
        destination.executemany(
            f"INSERT INTO public_gene_support ({','.join(fields)}) VALUES ({','.join('?' for _ in fields)})",
            [tuple(row[field] for field in fields) for row in summary_rows],
        )
        destination.execute("CREATE INDEX IF NOT EXISTS idx_edges_subject_public ON edges(subject_node_id,independence_class)")
        destination.execute("CREATE INDEX IF NOT EXISTS idx_edges_object_public ON edges(object_node_id,independence_class)")
        destination.execute("CREATE INDEX IF NOT EXISTS idx_evidence_edge_public ON evidence(edge_id,independence_class)")
        destination.execute("CREATE INDEX IF NOT EXISTS idx_public_aliases_alias ON public_aliases(alias)")
    destination.close()
    return dict(inserted)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-release", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    baseline = Path(args.baseline_release).resolve()
    work = Path(args.work_dir).resolve()
    output = Path(args.output_dir).resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite existing release: {output}")
    output.mkdir(parents=True)

    components = sorted(path for path in (work / "normalized").iterdir() if path.is_dir())
    public = {table: [component / f"{table}.tsv.gz" for component in components if (component / f"{table}.tsv.gz").exists()]
              for table in TABLE_IDS}
    aliases = [component / "aliases.tsv.gz" for component in components if (component / "aliases.tsv.gz").exists()]
    merged = {}
    duplicates = {}
    for table, id_field in TABLE_IDS.items():
        merged[table], duplicates[table] = merge_table(baseline / f"{table}.tsv.gz", public[table], output / f"{table}.tsv.gz", id_field)

    targets = set()
    for row in gz_rows(baseline / "nodes.tsv.gz"):
        if row["node_type"] == "Gene" and row["canonical_id"].startswith("Zm00001eb"):
            targets.add(row["canonical_id"])
    summary_rows, public_predicates = build_public_summary(targets, public["edges"], output / "public_gene_support.tsv.gz")

    for filename in ("gene_cards.tsv.gz", "mechanism_paths.tsv.gz", "top_unlabeled_candidates.tsv"):
        shutil.copy2(baseline / filename, output / filename)
    shutil.copytree(baseline / "gene_cards", output / "gene_cards")
    shutil.copy2(work / "reports" / "source_registry.json", output / "public_source_registry.json")
    report_dir = output / "public_reports"
    shutil.copytree(work / "reports", report_dir)
    code_dir = output / "code"
    code_dir.mkdir()
    for path in sorted((work / "code").glob("*.py")):
        shutil.copy2(path, code_dir / path.name)

    sqlite_inserted = build_sqlite(baseline / "mekg.sqlite", output / "mekg.sqlite", public, aliases, summary_rows)
    tier_counts = Counter(row["public_support_tier"] for row in summary_rows)
    manifest = {
        "release_id": "MEKG-public-v2-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
        "created_at": datetime.now(timezone.utc).isoformat(), "baseline_release": str(baseline),
        "evidence_layer": PUBLIC, "used_for_training": 0,
        "model_predictions_or_performance_changed": False,
        "counts": merged, "deduplicated_rows": duplicates,
        "public_normalized_components": [component.name for component in components],
        "public_sqlite_rows_inserted": sqlite_inserted,
        "public_support_tiers": dict(sorted(tier_counts.items())),
        "public_predicate_gene_endpoint_counts": dict(sorted(public_predicates.items())),
        "source_registry_sha256": sha256(output / "public_source_registry.json"),
        "known_source_gaps": {
            "MapMan": "The official X4.4 Zea mays download returned a license/login HTML page, not a mapping table; no MapMan gene-bin edges are released.",
            "MaizeMine": "The legacy MaizeMine v15 service closed the connection; MaizeGDB bulk B73 v5 annotation and mutant resources are used instead.",
        },
        "files": {},
        "semantic_guards": [
            "Every public record is PUBLIC_INDEPENDENT_SUPPORT and used_for_training=0.",
            "Public support never changes model probability, model rank, training split, or benchmark metrics.",
            "public_support_score_explanatory is evidence breadth only and is not a calibrated probability.",
            "Cross-species, GWAS, QTL, pathway, expression and chromatin evidence retain their non-causal scope.",
        ],
    }
    for path in sorted(output.iterdir()):
        if path.is_file():
            manifest["files"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    (output / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    (output / "README.md").write_text(
        "# MEKG public evidence v2\n\n"
        "This release extends the unchanged MEKG-small v1.1 model/KG outputs with a separate "
        "`PUBLIC_INDEPENDENT_SUPPORT` layer. Public data never participate in model training, "
        "model selection or performance reporting.\n\n"
        "Use `code/query_mekg_public_v2.py --db mekg.sqlite --gene Zm00001eb...` to combine the "
        "original model card with external evidence. `public_support_score_explanatory` summarizes "
        "evidence breadth only; it is not a predicted probability and must not be used to claim causality.\n",
        encoding="utf-8",
    )
    print(json.dumps({"output": str(output), "counts": merged, "tiers": dict(tier_counts),
                      "sqlite_inserted": sqlite_inserted}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
