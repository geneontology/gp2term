# from ontobio.io.gpadparser import GpadParser
# from ontobio.io.assocwriter import GpadWriter
from ontobio.io.gafparser import GafParser
from ontobio.io.assocparser import AssocParserConfig
from ontobio.io.assocwriter import GafWriter
from ontobio.ontol_factory import OntologyFactory
import os
import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filename')
parser.add_argument('-d', '--directory', help="Process a whole directory of files")
parser.add_argument('-m', '--match_output_suffix')
parser.add_argument('-l', '--leftover_output')
parser.add_argument('-o', '--ontology', help="Filename of GO ontology file to use")

aspect_lookup = {}

def get_ancestors(go_term):
    ### BFO:0000050 = part of
    ### I assume "subClassOf" = is a?
    all_ancestors = onto.ancestors(go_term)
    all_ancestors.append(go_term)
    subont = onto.subontology(all_ancestors)
    return subont.ancestors(go_term, relations=["subClassOf", "BFO:0000050"])

def is_biological_process(go_term):
    bp_root = "GO:0008150"
    if go_term == bp_root:
        return True
    ancestors = get_ancestors(go_term)
    if bp_root in ancestors:
        return True
    else:
        return False

def is_molecular_function(go_term):
    mf_root = "GO:0003674"
    if go_term == mf_root:
        return True
    ancestors = get_ancestors(go_term)
    if mf_root in ancestors:
        return True
    else:
        return False

def is_cellular_component(go_term):
    cc_root = "GO:0005575"
    if go_term == cc_root:
        return True
    ancestors = get_ancestors(go_term)
    if cc_root in ancestors:
        return True
    else:
        return False

def get_go_aspect(term):
    if is_molecular_function(term):
        return 'MF'
    elif is_biological_process(term):
        return 'BP'
    elif is_cellular_component(term):
        return 'CC'
    return None

def lookup_aspect(go_term):
    mf_root = "GO:0003674"
    if go_term in aspect_lookup["MF"] or go_term == mf_root:
        return "MF"
    bp_root = "GO:0008150"
    if go_term in aspect_lookup["BP"] or go_term == bp_root:
        return "BP"
    cc_root = "GO:0005575"
    if go_term in aspect_lookup["CC"] or go_term == cc_root:
        return "CC"

def get_aspect(term):
    if term not in aspect_lookup:
        aspect_lookup[term] = get_go_aspect(term)
    return aspect_lookup[term]

def match_in_with_col(annot, prefix):
    with_support_from = annot["evidence"]["with_support_from"]
    matches = []
    for w in with_support_from:
        if w.startswith(prefix):
            matches.append(w)
    return matches

def get_iprs(annot):
    return match_in_with_col(annot, "InterPro:IPR")

def get_ptns(annot):
    return match_in_with_col(annot, "PANTHER:PTN")

def annots_by_subject(annots, subject_id):
    found_annots = []
    for a in annots:
        if a["subject"]["id"] == subject_id:
            found_annots.append(a)
    return found_annots

def first_qualifier(annot):
    qualifiers = annot.get('qualifiers')
    if qualifiers and len(qualifiers) > 0:
        return qualifiers[0]

if __name__ == "__main__":
    args = parser.parse_args()

    if args.ontology:
        onto = OntologyFactory().create(args.ontology)
    else:
        onto = OntologyFactory().create("go")

    # Needed to retain IBA annotations
    config = AssocParserConfig()
    config.paint = True

    dir = 'resources/'
    if args.filename:
        filenames = [args.filename]
    else:
        filenames = []
        for f in os.listdir(dir):
            filenames.append(dir + f)
    annots = []
    all_annots = []
    leftover_annots = []
    # grouped_annots = {}
    # annots = {
    #     "IEA": [],
    #     "IBA": [],
    # }
    for f in filenames:
        iba_count = 0
        grouped_annots = {}
        print(f)
        mod = f.split(".")[-1]
        # annots[mod] = {
        #     "IEA": [],
        #     "IBA": []
        # }
        mod_annots = GafParser(config).parse(f, skipheader=True)
        all_annots = all_annots + mod_annots
        for a in mod_annots:
            term = a["object"]["id"]
            aspect = get_aspect(term)
            # Cache aspect lookup in ontology and reuse - probably small list
            # aspect = lookup_aspect(term)
            using_annot = False
            if aspect == "BP":
                with_matches = []
                evidence_code = a["evidence"]["type"]
                if evidence_code == "IEA":
                    with_matches = get_iprs(a)
                elif evidence_code == "IBA":
                    with_matches = get_ptns(a)
                    iba_count += len(with_matches)
                if len(with_matches) > 0:
                    annots.append(a)
                    for m in with_matches:
                        using_annot = True
                        if m not in grouped_annots:
                            grouped_annots[m] = {
                                "BP": [a],
                                "MF": []
                            }
                        else:
                            grouped_annots[m]["BP"].append(a)
                if not using_annot:
                    leftover_annots.append(a)
        # print("Total:", len(all_annots))
        # print("Leftovers:", len(leftover_annots))
        print("BP IBA withs count:", iba_count)
        iba_count = 0

        # Find MF annots for BP-derived IPRs
        for a in mod_annots:
            term = a["object"]["id"]
            aspect = get_aspect(term)
            if aspect == "MF":
                with_matches = []
                evidence_code = a["evidence"]["type"]
                if evidence_code == "IEA":
                    with_matches = get_iprs(a)
                elif evidence_code == "IBA":
                    with_matches = get_ptns(a)
                    iba_count += len(with_matches)
                for m in with_matches:
                    if m in grouped_annots:
                        grouped_annots[m]["MF"].append(a)
        print("MF IBA withs count:", iba_count)
        iba_count = 0

        # Filter out IPRs with no MF annots from last step
        filtered_annots = {}
        match_rows = []
        base_f = os.path.basename(f)
        match_outfile = base_f + "_matches.tsv"
        if args.match_output_suffix:
            match_outfile = "{}_{}".format(base_f, args.match_output_suffix)
        with open(match_outfile, 'w') as mof:
            writer = csv.writer(mof, delimiter="\t")
            for k in grouped_annots:
                if len(grouped_annots[k]["MF"]) > 0:
                    # Keep it
                    filtered_annots[k] = grouped_annots[k]
                    for a in grouped_annots[k]["BP"]:
                        gene_id = a["subject"]["id"]
                        gene_id_bits = gene_id.split(":")
                        mf_annots = annots_by_subject(grouped_annots[k]["MF"], gene_id)
                        if len(mf_annots) == 0:
                            # Should probably add this BP annot back to unused list
                            if a not in leftover_annots:
                                leftover_annots.append(a)
                            continue
                        for mfa in mf_annots:
                            id_ns = gene_id_bits[0]
                            id = gene_id_bits[-1]
                            gene_symbol = a["subject"]["label"]
                            relation = first_qualifier(a)
                            bp_term = a["object"]["id"]
                            bp_term_label = onto.label(bp_term)
                            bp_evidence_code = a["evidence"]["type"]
                            bp_reference = ",".join(a["evidence"]["has_supporting_reference"])
                            bp_assigned_by = a["provided_by"]
                            mf_term = mfa["object"]["id"]
                            mf_term_label = onto.label(mf_term)
                            mf_evidence_code = mfa["evidence"]["type"]
                            mf_reference = ",".join(mfa["evidence"]["has_supporting_reference"])
                            mf_assigned_by = mfa["provided_by"]
                            out_fields = [k, id_ns, id, gene_symbol, relation,
                                          bp_term, bp_term_label, bp_evidence_code, bp_reference, bp_assigned_by,
                                          mf_term, mf_term_label, mf_evidence_code, mf_reference, mf_assigned_by]
                            match_rows.append(out_fields)

            # L.sort(key=lambda k: (k[0], -k[1]), reverse=True)
            match_rows.sort(key=lambda k: k[2])
            for mr in match_rows:
                writer.writerow(mr)

    print("Total:", len(all_annots))
    print("Leftovers:", len(leftover_annots))
    outfile = "leftovers.gaf"
    if args.leftover_output:
        outfile = args.leftover_output
    with open(outfile, "w") as lf:
        gaf_writer = GafWriter(lf)
        for a in leftover_annots:
            gaf_writer.write_assoc(a)