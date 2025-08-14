#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_astral_map.py
------------------
Vasculha uma pasta com múltiplos arquivos de alinhamento (NEXUS e FASTA),
extrai todos os rótulos de táxons e gera um map file para o ASTRAL.

Formatos suportados (detecção por extensão e/ou conteúdo):
- NEXUS: .nex, .nexus  (TAXLABELS ou DATA/MATRIX, incluindo interleaved)
- FASTA: .fasta, .fa, .fas  (rótulo após '>')

Uso básico:
    python make_astral_map.py --input ./alinhamentos --out-map astral.map

Com uma tabela de grupos (CSV/TSV com colunas "taxon,group"):
    python make_astral_map.py --input ./alinhamentos --groups grupos.csv --out-map astral.map

Opções úteis:
    --pattern "*.nex,*.nexus,*.fasta,*.fa,*.fas"
    --default-group species|NA|none  (padrão: species)
    --out-taxa taxa_list.txt         (salva lista única de táxons)
    --strict                         (falha se algum arquivo não puder ser lido)
    --verbose                        (mensagens extras no stderr)

Não requer bibliotecas externas.
"""

import sys
from pathlib import Path
import re
import csv
import argparse
from typing import Iterable, Set, Dict, List

COMMENT_RE = re.compile(r"\[.*?\]", flags=re.DOTALL)  # remove comentários NEXUS [ ... ]
TOKEN_RE = re.compile(r"'([^']*)'|\"([^\"]*)\"|(\S+)")  # tokens respeitando aspas simples/duplas


def strip_nexus_comments(text: str) -> str:
    return re.sub(COMMENT_RE, "", text)


def tokenize(line: str) -> List[str]:
    # Retorna tokens preservando conteúdo dentro de aspas
    toks: List[str] = []
    for m in TOKEN_RE.finditer(line.strip()):
        tok = next(g for g in m.groups() if g is not None)
        toks.append(tok)
    return toks


def parse_taxlabels_block(text: str) -> List[str]:
    """
    Procura por bloco TAXA/TAXLABELS e retorna rótulos até o ';'
    """
    labels: List[str] = []
    lowered = text.lower()
    idx = lowered.find("taxlabels")
    if idx == -1:
        return labels
    # recorta desde taxlabels até o primeiro ';' subsequente
    sub = text[idx:]
    if ";" not in sub:
        return labels
    sub = sub[: sub.find(";")]
    # remove a palavra taxlabels
    sub = re.sub(r"(?i)^taxlabels", "", sub, count=1).strip()
    # tokeniza respeitando aspas
    labels = tokenize(sub)
    return labels


def parse_matrix_block(text: str) -> List[str]:
    """
    Procura bloco MATRIX dentro de BEGIN DATA; ... MATRIX ... ; e extrai o primeiro token de cada linha
    até o ponto-e-vírgula terminador do bloco MATRIX.
    Funciona para matrizes intercaladas (interleaved) também.
    """
    labels: List[str] = []
    lowered = text.lower()
    midx = lowered.find("matrix")
    if midx == -1:
        return labels

    sub = text[midx + len("matrix"):]
    # corta até o primeiro ';' (fim do bloco MATRIX)
    if ";" in sub:
        sub = sub[: sub.find(";")]

    for raw in sub.splitlines():
        line = raw.strip()
        if not line:
            continue
        # ignora linhas de cabeçalho/format etc (em alguns arquivos MATRIX pode repetir)
        if line.lower().startswith(("matrix", "format", "dimensions", "end", "begin")):
            continue
        toks = tokenize(line)
        if not toks:
            continue
        labels.append(toks[0])
    return labels


def parse_nexus_taxa(path: Path) -> Set[str]:
    """
    Tenta extrair rótulos de táxons de um arquivo NEXUS,
    priorizando TAXLABELS; se não houver, tenta MATRIX.
    """
    try:
        raw = path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        raw = path.read_text(encoding="latin-1", errors="ignore")
    no_comments = strip_nexus_comments(raw)

    taxa: Set[str] = set()

    # 1) Tenta TAXLABELS (mais fiel, independe de interleaving)
    labels = parse_taxlabels_block(no_comments)
    if labels:
        taxa.update(labels)

    # 2) Reforço: extrai de MATRIX também (para casos sem TAXLABELS)
    if not taxa:
        taxa.update(parse_matrix_block(no_comments))

    # limpeza final: remove ; perdidos ou espaços
    cleaned = set(l.strip().rstrip(";") for l in taxa if l.strip())
    return cleaned


def parse_fasta_taxa(path: Path) -> Set[str]:
    taxa: Set[str] = set()
    with path.open(encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                label = line[1:].strip()
                if not label:
                    continue
                # pega o primeiro token do cabeçalho para nome do táxon
                label = label.split()[0]
                taxa.add(label)
    return taxa


def detect_format_from_content(path: Path) -> str:
    # Mini detecção: olha primeiras linhas não vazias
    try:
        with path.open(encoding="utf-8", errors="ignore") as f:
            for _ in range(50):
                line = f.readline()
                if not line:
                    break
                s = line.strip()
                if not s:
                    continue
                if s.startswith(">"):
                    return "fasta"
                if s.upper().startswith("#NEXUS"):
                    return "nexus"
    except Exception:
        pass
    # fallback pela extensão
    ext = path.suffix.lower()
    if ext in (".fasta", ".fa", ".fas"):
        return "fasta"
    if ext in (".nex", ".nexus"):
        return "nexus"
    return "unknown"


def load_groups_csv(csv_path: Path) -> Dict[str, str]:
    """
    Lê CSV/TSV com mapeamento de 'taxon' -> 'group'.
    - Aceita cabeçalho com colunas 'taxon' e 'group' (case-insensitive).
    - Aceita arquivos sem cabeçalho com 2 colunas (taxon, group).
    - Delimitador autodetectado para vírgula ou tab.
    """
    mapping: Dict[str, str] = {}

    sample = csv_path.read_text(encoding="utf-8", errors="ignore")
    # Detecta delimitador básico
    delim = "\t" if sample.count("\t") > sample.count(",") else ","

    rows = list(csv.reader(sample.splitlines(), delimiter=delim))
    if not rows:
        return mapping

    # tenta cabeçalho
    header = [c.strip().lower() for c in rows[0]]
    start_idx = 1
    if "taxon" in header and "group" in header:
        tax_idx = header.index("taxon")
        grp_idx = header.index("group")
    else:
        # sem cabeçalho: assume 2 colunas
        tax_idx, grp_idx = 0, 1
        start_idx = 0

    for r in rows[start_idx:]:
        if len(r) <= max(tax_idx, grp_idx):
            continue
        taxon = r[tax_idx].strip()
        group = r[grp_idx].strip()
        if taxon:
            mapping[taxon] = group
    return mapping


def gather_taxa(paths: Iterable[Path], strict: bool=False, verbose: bool=False) -> Set[str]:
    taxa: Set[str] = set()
    for p in paths:
        try:
            fmt = detect_format_from_content(p)
            if fmt == "fasta":
                ts = parse_fasta_taxa(p)
            elif fmt == "nexus":
                ts = parse_nexus_taxa(p)
            else:
                if verbose:
                    print(f"[ignorado] {p} (formato desconhecido)", file=sys.stderr)
                if strict:
                    raise ValueError(f"Formato desconhecido para {p}")
                continue
            if verbose:
                print(f"[ok] {p} -> {len(ts)} táxon(s)", file=sys.stderr)
            taxa.update(ts)
        except Exception as e:
            msg = f"[aviso] Falha ao ler {p}: {e}"
            if strict:
                raise
            else:
                print(msg, file=sys.stderr)
    return taxa


def main():
    ap = argparse.ArgumentParser(description="Gera map file para ASTRAL a partir de múltiplos alinhamentos (NEXUS/FASTA)")
    ap.add_argument("--input", required=True, help="Diretório com alinhamentos")
    ap.add_argument("--pattern", default="*.nex,*.nexus,*.fasta,*.fa,*.fas",
                    help="Padrões de glob separados por vírgula (padrão: *.nex,*.nexus,*.fasta,*.fa,*.fas)")
    ap.add_argument("--groups", default=None, help="CSV/TSV com colunas 'taxon,group' (opcional)")
    ap.add_argument("--default-group", default="species", choices=["species", "NA", "none"],
                    help="Se não houver grupos, usa: 'species' (taxon mapeia pra si mesmo), 'NA' ou 'none' (string vazia).")
    ap.add_argument("--out-map", required=True, help="Caminho do map file de saída (TSV: taxon\\tgroup)")
    ap.add_argument("--out-taxa", default=None, help="(Opcional) Salva lista única de táxons detectados")
    ap.add_argument("--strict", action="store_true", help="Falha se algum arquivo não puder ser lido ou tiver formato desconhecido")
    ap.add_argument("--verbose", action="store_true", help="Mensagens detalhadas durante a varredura")
    args = ap.parse_args()

    root = Path(args.input).expanduser().resolve()
    if not root.is_dir():
        ap.error(f"--input não é um diretório: {root}")

    patterns = [p.strip() for p in args.pattern.split(",") if p.strip()]
    files: List[Path] = []
    for pat in patterns:
        files.extend(root.rglob(pat))

    if not files:
        print("Nenhum arquivo encontrado com os padrões fornecidos.", file=sys.stderr)
        sys.exit(2)

    if args.verbose:
        print(f"Varrendo {len(files)} arquivo(s)...", file=sys.stderr)
    taxa = gather_taxa(files, strict=args.strict, verbose=args.verbose)
    if not taxa:
        print("Nenhum táxon encontrado. Verifique seus arquivos de alinhamento.", file=sys.stderr)
        sys.exit(3)

    if args.verbose:
        print(f"Detectados {len(taxa)} táxon(s) únicos.", file=sys.stderr)

    mapping: Dict[str, str] = {}
    if args.groups:
        mapping = load_groups_csv(Path(args.groups))
        if args.verbose:
            print(f"Carregados {len(mapping)} mapeamentos de grupos.", file=sys.stderr)

    default_group_mode = args.default_group
    def pick_group(t: str) -> str:
        if t in mapping:
            return mapping[t]
        if default_group_mode == "species":
            return t
        if default_group_mode == "NA":
            return "NA"
        if default_group_mode == "none":
            return ""
        return t  # fallback

    out_map = Path(args.out_map).expanduser().resolve()
    out_map.parent.mkdir(parents=True, exist_ok=True)
    with out_map.open("w", encoding="utf-8") as f:
        for t in sorted(taxa):
            f.write(f"{t}\t{pick_group(t)}\n")

    if args.out_taxa:
        out_taxa = Path(args.out_taxa).expanduser().resolve()
        out_taxa.parent.mkdir(parents=True, exist_ok=True)
        with out_taxa.open("w", encoding="utf-8") as tf:
            for t in sorted(taxa):
                tf.write(f"{t}\n")

    if args.verbose:
        print(f"OK: mapa salvo em {out_map}", file=sys.stderr)
        if args.out_taxa:
            print(f"Lista de táxons salva em {out_taxa}", file=sys.stderr)


if __name__ == "__main__":
    main()
