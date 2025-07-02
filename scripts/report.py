import pandas as pd
import re
import argparse
import sys

def short_to_hgvs(var):
    aa1to3 = {
        "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
        "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
        "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
        "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
        "*": "Ter"
    }
    m = re.match(r"([A-Z*])([0-9]+)([A-Z*])", var)
    if m:
        ref, pos, alt = m.groups()
        if ref in aa1to3 and alt in aa1to3:
            return f"p.{aa1to3[ref]}{pos}{aa1to3[alt]}"
    return var

def clinvar_lookup(variant, clinvar38):
    hgvs = short_to_hgvs(variant)
    matches = clinvar38[clinvar38['Name'].str.contains(hgvs, na=False) & clinvar38['Name'].str.contains('GBA1', na=False)]
    if not matches.empty:
        row = matches.iloc[0]
        return (
            str(row['Name']) if 'Name' in row else '',
            ";".join(matches['ClinicalSignificance'].unique()),
            str(row['Chromosome']) if 'Chromosome' in row else '',
            str(row['PositionVCF']) if 'PositionVCF' in row else '',
            str(row['ReferenceAlleleVCF']) if 'ReferenceAlleleVCF' in row else '',
            str(row['AlternateAlleleVCF']) if 'AlternateAlleleVCF' in row else ''
        )
    else:
        return ("Not found", "Not found", "", "", "", "")

def main():
    parser = argparse.ArgumentParser(
        description="Annotate Gauchian TSV with ClinVar info.",
        epilog="Usage: python report.py -i test.gauchian.tsv -c variant_summary.txt.gz -o gauchian_tsv_with_clinvar.tsv"
    )
    parser.add_argument("-i", "--tsv", required=True, help="Gauchian output tsv file")
    parser.add_argument("-c", "--clinvar", required=True, help="Clinvar variant summary table: variant_summary.txt.gz")
    parser.add_argument("-o", "--output", required=True, help="Output file name.")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    GAUCHIAN_TSV = args.tsv
    CLINVAR_SUMMARY = args.clinvar
    OUTPUT_TSV = args.output

    clinvar = pd.read_csv(CLINVAR_SUMMARY, sep="\t", low_memory=False)
    clinvar38 = clinvar[clinvar['Assembly'] == 'GRCh38']

    df = pd.read_csv(GAUCHIAN_TSV, sep="\t")

    clinvar_names = []
    clinvar_pathogenicities = []
    chr_list = []
    pos_list = []
    ref_list = []
    alt_list = []
    for idx, row in df.iterrows():
        variants = set()
        v1 = row.get('GBAP1-like_variant_exon9-11', '')
        if pd.notnull(v1) and v1 and v1 != 'None':
            for v in str(v1).strip("/").split("/"):
                if v: variants.add(v)
        v2 = row.get('other_unphased_variants', '')
        if pd.notnull(v2) and v2 and v2 != 'None':
            for v in str(v2).strip("/").split("/"):
                if v: variants.add(v)
        if variants:
            names = []
            pathos = []
            chrs = []
            poss = []
            refs = []
            alts = []
            for v in variants:
                name, pathogenicity, chrom, pos, ref, alt = clinvar_lookup(v, clinvar38)
                names.append(name)
                pathos.append(pathogenicity)
                chrs.append(chrom)
                poss.append(pos)
                refs.append(ref)
                alts.append(alt)
            clinvar_names.append(";".join(names))
            clinvar_pathogenicities.append(";".join(pathos))
            chr_list.append(";".join(chrs))
            pos_list.append(";".join(poss))
            ref_list.append(";".join(refs))
            alt_list.append(";".join(alts))
        else:
            clinvar_names.append("")
            clinvar_pathogenicities.append("")
            chr_list.append("")
            pos_list.append("")
            ref_list.append("")
            alt_list.append("")

    df['ClinVar_Name'] = clinvar_names
    df['ClinVar_Pathogenicity'] = clinvar_pathogenicities
    df['ClinVar_Chr'] = chr_list
    df['ClinVar_Pos'] = pos_list
    df['ClinVar_Ref'] = ref_list
    df['ClinVar_Alt'] = alt_list

    df.to_csv(OUTPUT_TSV, sep="\t", index=False)

    print(f"Annotation complete. Results written to {OUTPUT_TSV}")

if __name__ == "__main__":
    main()