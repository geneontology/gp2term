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
parser.add_argument('-e', '--log_filepath', help="destination of log file")

aspect_lookup = {}
log_file = None
EVI_CODE_TO_WITH_PREFIX = dict(IEA="InterPro:IPR", IBA="PANTHER:PTN")

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

def file_away(annot_dict, annot):
    ### Data structure:
    # annot_dict = {
    #   "IEA": {
    #       "InterPro:IPR58694": [annot1, annot2],
    #       "InterPro:IPR37482": [annot3, annot4],
    #   },
    #   "IBA": {
    #       "PANTHER:PTN58694": [annot5, annot6],
    #       "PANTHER:PTN95847": [annot7, annot8],
    #   }
    # }
    ### TODO: Need method for uniq-ing annot lists due to multiple with matches
    ### TODO: Need method for pulling BP vs MF annots
    evi_code = annot["evidence"]["type"]
    using_annot = False  # don't know yet
    if evi_code not in EVI_CODE_TO_WITH_PREFIX:
        return annot_dict, using_annot
    matches = match_in_with_col(annot, EVI_CODE_TO_WITH_PREFIX[evi_code])
    for m in matches:
        using_annot = True
        if evi_code not in annot_dict:
            # Really should be subkeying with value. Which withs? Create mapping of evi_code to with prefixes
            annot_dict[evi_code] = {m: [annot]}
        elif m not in annot_dict[evi_code]:
            annot_dict[evi_code][m] = [annot]
        else:
            annot_dict[evi_code][m].append(annot)
    return annot_dict, using_annot

def flatten_with_dict(with_dict, uniqify=False):
    all_annots_list = []
    for with_value in with_dict:
        for a in with_dict[with_value]:
            all_annots_list.append(a)
    if uniqify:
        temp_list = []
        for a in all_annots_list:
            if a not in temp_list:
                temp_list.append(a)
        all_annots_list = temp_list
            # if uniqify and a not in all_annots_list:
            #     all_annots_list.append(a)
            # elif not uniqify:
            #     all_annots_list.append(a)
    return all_annots_list

def match_aspect(annots, aspect):
    matching_annots = []
    for annot in annots:
        if annot["aspect"] == aspect:
            matching_annots.append(annot)
    return matching_annots

def log(message):
    print(message)
    log_file.write(message + "\n")

class Gp2TermGafParser(GafParser):
    def __init__(self):
        # Needed to retain IBA annotations
        config = AssocParserConfig()
        config.paint = True
        GafParser.__init__(self, config=config)

class MatchMaker():
    def __init__(self, ontology):
        self.all_annots = []
        self.parser = Gp2TermGafParser()
        self.ontology = ontology
        self.all_promoted_annots = []

    def find_iea_iba_matches(self, filename):
        all_bp_ev_counts = {}
        grouped_annots = {}
        leftover_annots = []
        log(filename)
        mod_annots = self.parser.parse(filename, skipheader=True)
        self.all_annots = self.all_annots + mod_annots
        for a in mod_annots:
            term = a["object"]["id"]
            aspect = a["aspect"]
            if aspect == "P" or aspect == "F":
                grouped_annots, using_annot = file_away(grouped_annots, a)
                if aspect == "P":
                    if not using_annot:
                        leftover_annots.append(a)
                    evidence_code = a["evidence"]["type"]
                    if evidence_code not in all_bp_ev_counts:
                        all_bp_ev_counts[evidence_code] = 1
                    else:
                        all_bp_ev_counts[evidence_code] += 1

        dismissed_annots = []
        # match_rows = []
        base_f = os.path.basename(f)
        match_outfile = base_f + "_matches.tsv"
        if args.match_output_suffix:
            match_outfile = "{}.{}.tsv".format(base_f, args.match_output_suffix)
        with open(match_outfile, 'w') as mof:
            writer = csv.writer(mof, delimiter="\t")
            for ec in grouped_annots:
                match_rows = []
                ### For each evi_code, count unique annots that have with_matches (flatten dict)
                log("BP {} withs count: {}".format(ec, len(match_aspect(flatten_with_dict(grouped_annots[ec], uniqify=True), 'P'))))
                ### Loop through with_value annots, segregate BPs from MFs, if len(BPs) > 0 and len(MFs) > 0 this with_value set gets written out
                for with_value in grouped_annots[ec]:
                    bp_annots = match_aspect(grouped_annots[ec][with_value], 'P')
                    mf_annots = match_aspect(grouped_annots[ec][with_value], 'F')
                    if len(bp_annots) < 1:
                        grouped_annots[ec][with_value] = []  # Delete this key
                    elif len(mf_annots) < 1:
                        dismissed_annots = dismissed_annots + bp_annots  # Cleanup (uniqify, remove annots promoted elsewhere) later
                        grouped_annots[ec][with_value] = []  # Delete this key
                    else:  # Continue on promoting
                        for a in bp_annots:
                            gene_id = a["subject"]["id"]
                            gene_id_bits = gene_id.split(":")
                            id_ns = gene_id_bits[0]
                            id = gene_id_bits[-1]
                            # Find 'with-matched' MF annotations to same gene product
                            gene_mf_annots = annots_by_subject(mf_annots, gene_id)
                            if len(gene_mf_annots) == 0:
                                # Should probably add this BP annot back to unused list
                                if a not in leftover_annots:
                                    leftover_annots.append(a)
                                continue

                            gene_symbol = a["subject"]["label"]
                            relation = first_qualifier(a)
                            bp_term = a["object"]["id"]
                            bp_term_label = self.ontology.label(bp_term)
                            bp_evidence_code = a["evidence"]["type"]
                            bp_reference = ",".join(a["evidence"]["has_supporting_reference"])
                            bp_assigned_by = a["provided_by"]
                            for mfa in gene_mf_annots:
                                mf_term = mfa["object"]["id"]
                                mf_term_label = self.ontology.label(mf_term)
                                mf_evidence_code = mfa["evidence"]["type"]
                                mf_reference = ",".join(mfa["evidence"]["has_supporting_reference"])
                                mf_assigned_by = mfa["provided_by"]
                                out_fields = [with_value, id_ns, id, gene_symbol, relation,
                                              bp_term, bp_term_label, bp_evidence_code, bp_reference, bp_assigned_by,
                                              mf_term, mf_term_label, mf_evidence_code, mf_reference, mf_assigned_by]
                                match_rows.append(out_fields)
                match_rows.sort(key=lambda k: k[2])
                for mr in match_rows:
                    writer.writerow(mr)
        # print("Total:", len(all_annots))
        # print("Leftovers:", len(leftover_annots))

        all_promoted_annots = []
        for ev in grouped_annots:
            promoted_bp_annots = match_aspect(flatten_with_dict(grouped_annots[ev], uniqify=True), 'P')
            all_promoted_annots = all_promoted_annots + promoted_bp_annots
            log("{} {} BP annotations inputted".format(all_bp_ev_counts[ev], ev))
            # 5000 IEA BP annotations ‘involved in’
            log("{} {} BP annotations ‘involved in’".format(len(promoted_bp_annots), ev))
        # self.all_promoted_annots[filename] = all_promoted_annots
        self.all_promoted_annots = self.all_promoted_annots + all_promoted_annots

        ## Something going on below is super slow!
        ### Cleanup leftovers
        for da in dismissed_annots:
            if da not in leftover_annots and da not in all_promoted_annots:
                leftover_annots.append(da)

        # print("Leftovers:", len(leftover_annots))
        log("Leftovers: {}".format(len(leftover_annots)))
        outfile = base_f + "_leftovers.gaf"
        if args.leftover_output:
            outfile = "{}.{}_leftovers.gaf".format(base_f, args.leftover_output)
        with open(outfile, "w") as lf:
            gaf_writer = GafWriter(lf)
            for a in leftover_annots:
                gaf_writer.write_assoc(a)

    def match_experimentals(self, filename):
        ### Take in leftovers file
        annots = self.parser.parse(filename, skipheader=True)
        ### Filter for exp evidence codes (IMP, IDA, IGI, IPI, IEP, ISS)

        ### check GP, term, against self.all_promoted_annots. Do I need to keep these separated?
        if filename in self.all_promoted_annots:
            self.all_promoted_annots[filename]
        else:
            log("No promoted annotations entry for '{}'".format(filename))

# Group all annots by evi code
# Find matches between with values (Sames evi_code structure)
if __name__ == "__main__":
    args = parser.parse_args()

    if args.ontology:
        onto = OntologyFactory().create(args.ontology)
    else:
        onto = OntologyFactory().create("go")
    parser = Gp2TermGafParser()
    mmaker = MatchMaker(onto)

    log_filepath = "gp2t_log.txt"
    if args.log_filepath:
        log_filepath = args.log_filepath
    log_file = open(log_filepath, 'w+')

    dir = 'resources/'
    if args.directory:
        dir = args.directory
    if args.filename:
        filenames = [args.filename]
    else:
        filenames = []
        for f in os.listdir(dir):
            filenames.append(dir + f)

    for f in filenames:
        ### TODO: Methodize this
        mmaker.find_iea_iba_matches(f)

        #

    log("Total: {}".format(len(mmaker.all_annots)))

    log_file.close()
