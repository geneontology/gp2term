from ontobio.io.gpadparser import GpadParser
from ontobio.io.assocwriter import GpadWriter
from ontobio.ontol_factory import OntologyFactory
import os
import csv

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

def get_iprs(annot):
    with_support_from = annot["evidence"]["with_support_from"]
    iprs = []
    for w in with_support_from:
        if w.startswith("InterPro:IPR"):
            iprs.append(w)
    return iprs

def annot_by_subject(annots, subject_id):
    for a in annots:
        if a["subject"]["id"] == subject_id:
            return a

def first_qualifier(annot):
    qualifiers = annot.get('qualifiers')
    if qualifiers and len(qualifiers) > 0:
        return qualifiers[0]

onto = OntologyFactory().create("go")

dir = 'resources/'
annots = []
all_annots = []
leftover_annots = []
grouped_annots = {}
# annots = {
#     "IEA": [],
#     "IBA": [],
# }
for f in os.listdir(dir):
    print(f)
    mod = f.split(".")[0]
    # annots[mod] = {
    #     "IEA": [],
    #     "IBA": []
    # }
    mod_annots = GpadParser().parse(dir + f, skipheader=True)
    all_annots = all_annots + mod_annots
    for a in mod_annots:
        term = a["object"]["id"]
        aspect = get_aspect(term)
        # Cache aspect lookup in ontology and reuse - probably small list
        # aspect = lookup_aspect(term)
        using_annot = False
        if aspect == "BP":
            evidence_code = a["evidence"]["type"]
            if evidence_code == "IEA" or evidence_code == "IBA":
                iprs = get_iprs(a)
                if len(iprs) > 0:
                    annots.append(a)
                    for i in iprs:
                        using_annot = True
                        if i not in grouped_annots:
                            grouped_annots[i] = {
                                "BP": [a],
                                "MF": []
                            }
                        else:
                            grouped_annots[i]["BP"].append(a)
        if not using_annot:
            leftover_annots.append(a)
print("Total:", len(all_annots))
print("Leftovers:", len(leftover_annots))

# Find MF annots for BP-derived IPRs
for a in all_annots:
    term = a["object"]["id"]
    aspect = get_aspect(term)
    if aspect == "MF":
        iprs = get_iprs(a)
        for i in iprs:
            if i in grouped_annots:
                grouped_annots[i]["MF"].append(a)

# Filter out IPRs with no MF annots from last step
filtered_annots = {}
with open("matching_IPRs.tsv", 'w') as f:
    writer = csv.writer(f, delimiter="\t")
    for k in grouped_annots:
        if len(grouped_annots[k]["MF"]) > 0:
            # Keep it
            filtered_annots[k] = grouped_annots[k]
            for a in grouped_annots[k]["BP"]:
                gene_id = a["subject"]["id"]
                mf_annot = annot_by_subject(grouped_annots[k]["MF"], gene_id)
                if mf_annot is None:
                    # Should probably add this BP annot back to unused list
                    if a not in leftover_annots:
                        leftover_annots.append(a)
                    continue
                gene_id_bits = gene_id.split(":")
                id_ns = gene_id_bits[0]
                id = gene_id_bits[-1]
                relation = first_qualifier(a)
                bp_term = a["object"]["id"]
                bp_term_label = onto.label(bp_term)
                bp_evidence_code = a["evidence"]["type"]
                bp_reference = ",".join(a["evidence"]["has_supporting_reference"])
                bp_assigned_by = a["provided_by"]
                mf_term = mf_annot["object"]["id"]
                mf_term_label = onto.label(mf_term)
                mf_evidence_code = mf_annot["evidence"]["type"]
                mf_reference = ",".join(mf_annot["evidence"]["has_supporting_reference"])
                mf_assigned_by = mf_annot["provided_by"]
                out_fields = [k, id_ns, id, relation,
                              bp_term, bp_term_label, bp_evidence_code, bp_reference, bp_assigned_by,
                              mf_term, mf_term_label, mf_evidence_code, mf_reference, mf_assigned_by]
                writer.writerow(out_fields)

print("Leftovers:", len(leftover_annots))
with open("leftovers.gpad", "w") as lf:
    gpad_writer = GpadWriter(lf)
    for a in leftover_annots:
        gpad_writer.write_assoc(a)