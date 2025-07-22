import os,sys
import string
from optparse import OptionParser
import glob
import json

import libgly
import csvutil
import section_stats
import biomarker_util

__version__="1.0"
__status__ = "Dev"



def load_protein_names():

    file_list = glob.glob("reviewed/*_protein_*names.csv")
    file_list += glob.glob("reviewed/*_protein_*names_uniprotkb.csv")
    load_obj = {"recommended":{}, "synonyms":{}}
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            for f in f_list:
                name = row[f_list.index(f)]
                if name.strip() == "":
                    continue
                if f in ["recommended_name_full"]:
                    load_obj["recommended"][canon] = name
                else:
                    if canon not in load_obj["synonyms"]:
                        load_obj["synonyms"][canon] = []
                    if name not in load_obj["synonyms"][canon]:
                        load_obj["synonyms"][canon].append(name)


    return load_obj



def get_disease_info(do_id, mondo_id, mim_id):
    disease_info = {}
    json_file = ""
    if do_id != "":
        json_file = "jsondb/diseasedb/doid.%s.json" % (do_id.replace("DOID:", ""))
    elif mondo_id != "":
        json_file = "jsondb/diseasedb/mondo.%s.json" % (mondo_id.replace("MONDO:", ""))
    elif mim_id != "":
        json_file = "jsondb/diseasedb/mim.%s.json" % (mim_id.replace("OMIM:", ""))
    if os.path.isfile(json_file):
        disease_info = json.loads(open(json_file, "r").read())

    return disease_info












###############################
def main():

    global wrk_dir

    
    generated_dir = "/data/projects/glygen/generated/"
    wrk_dir = "/data/shared/repos/glygen-backend-integration/object-maker/"
    jsondb_dir = wrk_dir + "/jsondb/"

    
    global is_cited
    is_cited = libgly.get_is_cited()


    protein_name_dict = load_protein_names()

    map_dict = {}
    csvutil.load_dictionaries(map_dict, "generated/misc/")



    # The canonical ID is a combination of the "biomarker" and "assessed_biomarker_entity" fields. 
    # The isoform ID is a combination of the canonical ID and the condition ID fields. 


    evdn_dict, tag_dict = {}, {}

    seen_tissue = {}
    seen_evdn = {}
    biomarkerid2mid = {}
    #biomarkerid2doid = {}
    load_obj = {"seen":{}}
    file_list = ["reviewed/glycan_biomarkers.csv","reviewed/human_protein_biomarkers.csv"]
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet_as_dict(data_frame, in_file, ",", "biomarker_id")
        tmp_fl = data_frame["fields"]
        #f_list = ["assessed_biomarker_entity", "assessed_entity_type"]
        f_list = ["biomarker_canonical_id"]
        for main_id in data_frame["data"]:
            for tmp_row in data_frame["data"][main_id]:
                if main_id == "":
                    continue
                obj = {"biomarker_id":main_id}
                for f in f_list:
                    obj[f] = tmp_row[tmp_fl.index(f)]
                biomarker  = tmp_row[tmp_fl.index("biomarker")]
                best_biomarker_role = tmp_row[tmp_fl.index("best_biomarker_role")]
                loinc_code = tmp_row[tmp_fl.index("loinc_code")]
                do_desc = tmp_row[tmp_fl.index("do_desc")]
                do_syn = tmp_row[tmp_fl.index("do_syn")].split("|")
                evdn_type = tmp_row[tmp_fl.index("evidence_type")]
                src_xref_key = tmp_row[tmp_fl.index("src_xref_key")]
                src_xref_id = tmp_row[tmp_fl.index("src_xref_id")]
                src_xref_badge = libgly.get_xref_badge(map_dict, src_xref_key)
                src_xref_url =  libgly.get_xref_url(map_dict, src_xref_key, src_xref_id,is_cited)
                xref_key = tmp_row[tmp_fl.index("xref_key")]
                xref_id = tmp_row[tmp_fl.index("xref_id")]
                xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                m_id = tmp_row[tmp_fl.index("glytoucan_ac")] if "glytoucan_ac" in tmp_fl else ""
                m_id = tmp_row[tmp_fl.index("uniprotkb_canonical_ac")] if "uniprotkb_canonical_ac" in tmp_fl else m_id
                component_category = "glycan" if "glytoucan_ac" in tmp_fl else ""
                component_category = "protein" if "uniprotkb_canonical_ac" in tmp_fl else component_category
                src_xref_url = libgly.get_xref_url(map_dict, src_xref_key, src_xref_id,is_cited)
                src_xref_badge = libgly.get_xref_badge(map_dict, src_xref_key)
                rec_name = tmp_row[tmp_fl.index("assessed_biomarker_entity")]
                synonyms = []
                for syn in tmp_row[tmp_fl.index("component_synonyms")].split("|"):
                    if syn.strip() != "":
                        synonyms.append({"synonym":syn})
                cmp_tags = tmp_row[tmp_fl.index("component_level_tags")]
                top_tags = tmp_row[tmp_fl.index("top_level_tags")]
                tag_list_one = cmp_tags.split(";")
                tag_list_two = top_tags.split(";")
                note = tmp_row[tmp_fl.index("evidence")]
                note_obj = {"evidence":note}
                
    
                if main_id not in biomarkerid2mid:
                    biomarkerid2mid[main_id] = []
                biomarkerid2mid[main_id].append(m_id)
                do_id = tmp_row[tmp_fl.index("do_id")]
                condition = tmp_row[tmp_fl.index("condition")]
                tissue_id = tmp_row[tmp_fl.index("uberon_id")]
                tissue_name = tmp_row[tmp_fl.index("anatomical_entity")]
                tissue_obj = {}
                if tissue_id != "":
                    t_ns, t_nm = "UBERON", tissue_name
                    t_xref_key = "tissue_xref_uberon"
                    t_xref_url =  libgly.get_xref_url(map_dict, t_xref_key, tissue_id,is_cited)
                    tissue_obj = {"name":t_nm,"namespace":t_ns,"id":tissue_id,"url":t_xref_url,
                        "loinc_code":loinc_code} 
                disease_obj = {}    
                if do_id != "":
                    disease_obj = get_disease_info(do_id, "", "")
                    disease_obj["id"] = disease_obj["disease_id"]
                    disease_obj.pop("disease_id")
                    #disease_obj["recommended_name"]["name"] = condition
                    #disease_obj["recommended_name"]["description"] = do_desc
                if main_id not in load_obj:
                    load_obj[main_id] = obj
                    load_obj[main_id]["condition"] = disease_obj
                    load_obj[main_id]["biomarker_component"] = {}
                    load_obj[main_id]["best_biomarker_role"] = {}
                    load_obj[main_id]["evidence_source"] = []
                    load_obj[main_id]["crossref"] = {}
 
                load_obj[main_id]["best_biomarker_role"][best_biomarker_role] = True
                
                bco_id = src_xref_url.split("/")[-1]
                cmp_combo_id = "%s|%s" % (main_id,m_id)
                if cmp_combo_id not in load_obj[main_id]["biomarker_component"]:
                    load_obj[main_id]["biomarker_component"][cmp_combo_id]= {
                        "biomarker":biomarker,
                        "assessed_biomarker_entity_id":m_id, 
                        "assessed_entity_type":component_category, 
                        "assessed_biomarker_entity":{
                            "recommended_name": rec_name,
                            "synonyms": synonyms
                        },
                        "specimen": [],
                        "evidence_source":[]
                    }
                tissue_combo_id = "%s|%s|%s" % (main_id, m_id, tissue_id)
                if tissue_obj != {} and tissue_combo_id not in seen_tissue:
                    load_obj[main_id]["biomarker_component"][cmp_combo_id]["specimen"].append(tissue_obj)
                    seen_tissue[tissue_combo_id] = True

                evdn_obj_list = [{"id":bco_id, "database":src_xref_badge, "url":src_xref_url}]
                if xref_key.find("_xref_pubmed") != -1 or xref_key.find("_xref_doi") != -1:
                    evdn_obj_list.append({"id":xref_id, "database":xref_badge, "url":xref_url,"tags":[], 
                        "evidence_list":[]})
                
  

                for evdn_obj in evdn_obj_list:
                    cmp_evdn_comobo_id = "%s|%s|%s|%s" % (main_id,m_id,evdn_obj["id"],evdn_obj["database"])
                    top_evdn_comobo_id = "%s|%s|%s|%s" % (main_id,"top",evdn_obj["id"],evdn_obj["database"])
                    if cmp_evdn_comobo_id not in evdn_dict:
                        evdn_dict[cmp_evdn_comobo_id] = {}
                    if cmp_evdn_comobo_id not in tag_dict:
                        tag_dict[cmp_evdn_comobo_id] = {}
                    evdn_dict[cmp_evdn_comobo_id][note] = True
                    tag_dict[cmp_evdn_comobo_id][cmp_tags] = True
                    
                    if top_evdn_comobo_id not in evdn_dict:
                        evdn_dict[top_evdn_comobo_id] = {}
                    if top_evdn_comobo_id not in tag_dict:
                        tag_dict[top_evdn_comobo_id] = {}
                    evdn_dict[top_evdn_comobo_id][note] = True
                    tag_dict[top_evdn_comobo_id][top_tags] = True

                    if cmp_evdn_comobo_id not in seen_evdn and evdn_type in ["both_levels", "component_level"]:
                        load_obj[main_id]["biomarker_component"][cmp_combo_id]["evidence_source"].append(evdn_obj)
                        oo = load_obj[main_id]["biomarker_component"][cmp_combo_id]["evidence_source"][-1]
                        ooo = json.loads(json.dumps(oo))
                        load_obj[main_id]["biomarker_component"][cmp_combo_id]["evidence_source"][-1] = ooo
                        seen_evdn[cmp_evdn_comobo_id] = True
                    if top_evdn_comobo_id not in seen_evdn and evdn_type in ["both_levels", "top_level"]:
                        load_obj[main_id]["evidence_source"].append(evdn_obj)
                        oo = load_obj[main_id]["evidence_source"][-1]
                        ooo = json.loads(json.dumps(oo))
                        load_obj[main_id]["evidence_source"][-1] = ooo
                        seen_evdn[top_evdn_comobo_id] = True


    #--> publication
    pub_dict = {}
    file_list = ["reviewed/glycan_citations_biomarkers.csv","reviewed/human_protein_citations_biomarkers.csv"]
    seen_pub = {}
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet_as_dict(data_frame, in_file, ",", "biomarker_id")
        tmp_fl = data_frame["fields"]
        f_list =["title","journal_name","publication_date","authors","xref_key","xref_id","src_xref_key"]
        for main_id in data_frame["data"]:
            for tmp_row in data_frame["data"][main_id]:
                if main_id == "":
                    continue
                obj = {}
                for f in f_list:
                    obj[f] = tmp_row[tmp_fl.index(f)]
                pub_obj = {}
                xref_key, xref_id = obj["xref_key"], obj["xref_id"]
                xref_url =  libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                bco_id = src_xref_url.split("/")[-1]
                obj["publication_date"] = obj["publication_date"].split(" ")[0]
                xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                src_xref_url =  libgly.get_xref_url(map_dict, obj["src_xref_key"], main_id,is_cited)
                src_xref_badge = libgly.get_xref_badge(map_dict, src_xref_key)
                ev_obj = {"id":bco_id, "database":src_xref_badge, "url":src_xref_url}
                ref_obj = {"type":xref_badge,"id":xref_id,"url":xref_url}
                pub_obj = {
                    "title": obj["title"]
                    ,"journal": obj["journal_name"]
                    ,"authors": obj["authors"]
                    ,"date": obj["publication_date"]
                    ,"reference": [ref_obj]
                    ,"evidence":{bco_id:ev_obj}
                }
                ref_str = main_id + "|" + json.dumps(ref_obj)
                if ref_str not in seen_pub:
                    seen_pub[ref_str] = True
                    if main_id not in pub_dict:
                        pub_dict[main_id] = {}
                    pub_dict[main_id][xref_id] = pub_obj
                if main_id in pub_dict:
                    if xref_id in pub_dict[main_id]:
                        pub_dict[main_id][xref_id]["evidence"][bco_id] = ev_obj


    #--> crossref
    crossref_dict = {}
    file_list = ["reviewed/human_protein_biomarkers.csv", "reviewed/glycan_biomarkers.csv"]
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet_as_dict(data_frame, in_file, ",", "biomarker_id")
        tmp_fl = data_frame["fields"]
        f_list =["xref_key","xref_id"]
        for main_id in data_frame["data"]:
            for tmp_row in data_frame["data"][main_id]:
                if main_id == "":
                    continue
                if main_id not in crossref_dict:
                    crossref_dict[main_id] = {}
                xref_key, xref_id = "protein_xref_biomarkerkb", main_id
                url = "https://biomarkerkb.org/biomarker/%s" % (xref_id)
                o = {"id":xref_id, "database":"BiomarkerKB","url":url,"categories":["Biomarkers"]}
                crossref_dict[main_id][xref_id] = o





    for main_id in load_obj:
        if main_id == "seen":
            continue
        obj = load_obj[main_id]
        #if main_id in biomarkerid2doid:
        #    obj["disease"] = []
        #    for do_id in list(set(biomarkerid2doid[main_id])):
        #        obj["disease"].append(get_disease_info(do_id, "", ""))
        obj["crossref"] = []
        if main_id in crossref_dict:
            for xref_id in crossref_dict[main_id]:
                obj["crossref"].append(crossref_dict[main_id][xref_id])

        obj["citation"] = []
        if main_id in pub_dict:
            for xref_id in pub_dict[main_id]:
                o = pub_dict[main_id][xref_id]
                oo = {}
                for k in o:
                    oo[k] = o[k]
                oo["evidence"] = []
                for bco_id in o["evidence"]:
                    ev_obj = o["evidence"][bco_id]
                    oo["evidence"].append(ev_obj)
                obj["citation"].append(oo)

        tmp_list = []
        for k in obj["biomarker_component"]:
            o = obj["biomarker_component"][k]
            m_id = o["assessed_biomarker_entity_id"]
            for oo in o["evidence_source"]:
                if "evidence_list" not in oo:
                    continue
                xref_id, xref_badge = oo["id"], oo["database"]
                cmp_evdn_comobo_id = "%s|%s|%s|%s" % (main_id,m_id,xref_id,xref_badge)
                if cmp_evdn_comobo_id in evdn_dict:
                    for note in evdn_dict[cmp_evdn_comobo_id]:
                        oo["evidence_list"].append({"evidence":note})
                if cmp_evdn_comobo_id in tag_dict:
                    for cmp_tags in tag_dict[cmp_evdn_comobo_id]:
                        for t in cmp_tags.split(";"):
                            if t not in oo["tags"]:
                                oo["tags"].append(t)
            tmp_list.append(o)
        obj["biomarker_component"] = tmp_list

        for oo in obj["evidence_source"]:
            if "evidence_list" not in oo:
                continue
            xref_id, xref_badge = oo["id"], oo["database"]
            top_evdn_comobo_id = "%s|%s|%s|%s" % (main_id,"top",xref_id,xref_badge)
            if top_evdn_comobo_id in evdn_dict:
                for note in evdn_dict[top_evdn_comobo_id]:
                    oo["evidence_list"].append({"evidence":note})
            if top_evdn_comobo_id in tag_dict:
                for top_tags in tag_dict[top_evdn_comobo_id]:
                    for t in top_tags.split(";"):
                        if t not in oo["tags"]:
                            oo["tags"].append(t)



        tmp_list = []
        for r in obj["best_biomarker_role"]:
            for role in r.split(";"):
                tmp_list.append({"role":role})
        obj["best_biomarker_role"] = tmp_list

        tmp_list = []
        seen_evdn = {}
        section_stats.get_sec_stats(obj, "biomarker")

        out_file = "jsondb/biomarkerdb/%s.json" % (main_id)
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(obj, indent=4)))
        continue
 
        new_obj_list = biomarker_util.convert(obj)
        for obj in new_obj_list:
            main_id = obj["biomarker_canonical_id"] 
            out_file = "jsondb/biomarkerdb/%s.json" % (main_id) 
            #with open(out_file, "w") as FW:
            #    FW.write("%s\n" % (json.dumps(obj, indent=4)))



            

if __name__ == '__main__':
	main()

