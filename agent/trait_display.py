"""Public trait fields for the browser only; no inference or coordinate conversion."""
import math
import re


SCALAR_FIELDS = {
    'trait', 'title', 'study_title', 'study_id', 'PMID', 'population',
    'p_value', 'pvalue', 'p_value_adjusted', 'R2_percent', 'effect',
    'reported_gene_text', 'reported_location', 'chromosome', 'position',
    'start_position', 'stop_position', 'start', 'end', 'marker', 'marker_id',
    'model', 'year', 'window_bp', 'window_policy', 'unique_locus_count',
    'reference_genome', 'coordinate_system', 'coordinate_units', 'mapping_status',
    'directness', 'no_combined_pvalue', 'method_repeats_are_not_biological_replicates',
    'independent_study_count', 'causal_gene_proven', 'biological_causality_claim',
}
LOCUS_FIELDS = {'chromosome', 'position', 'pvalue', 'p_value', 'model', 'year', 'locus_id'}
PRIVATE_TEXT = re.compile(r'(/home/|/tmp/|[A-Z]:\\|sk-[a-zA-Z0-9])')


def is_trait_evidence(fact):
    typ = str(fact.get('evidence_type') or '').upper()
    return 'GWAS' in typ or 'QTL' in typ


def safe_scalar(value):
    if value is None or isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return value if math.isfinite(value) else None
    if isinstance(value, str) and len(value) < 2000 and not PRIVATE_TEXT.search(value):
        return value
    return None


def display_context(fact, attributes):
    """Do not infer a gene assignment from an unassigned marker/QTL record."""
    if not is_trait_evidence(fact) or not isinstance(attributes, dict):
        return {}
    result = {}
    for key in SCALAR_FIELDS:
        if key in attributes:
            value = safe_scalar(attributes[key])
            if value is not None:
                result[key] = value
    loci = attributes.get('locus_results')
    if isinstance(loci, list):
        result['locus_results'] = []
        for locus in loci:
            if not isinstance(locus, dict):
                continue
            clean = {key: safe_scalar(value) for key, value in locus.items() if key in LOCUS_FIELDS}
            clean = {key: value for key, value in clean.items() if value is not None}
            if clean:
                result['locus_results'].append(clean)
    return result


def display_order(facts):
    """Keep stronger records first; show trait associations before generic context."""
    if not any(is_trait_evidence(fact) for fact in facts):
        return facts
    strong_types = {'BINDING_AND_FUNCTIONAL_REPORTER_SUPPORT', 'BINDING_CHROMATIN_RNA_INTERSECTION'}
    def key(fact):
        if fact.get('category') in {'genetic_or_functional', 'molecular_response'} or fact.get('evidence_type') in strong_types:
            return 0
        return 1 if is_trait_evidence(fact) else 2
    return sorted(facts, key=key)
