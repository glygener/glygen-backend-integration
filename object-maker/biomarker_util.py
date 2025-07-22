
def convert(glygen_biomarker):

    result_data = []

    # Grab common fields that will be shared across the records after being split.
    biomarker_id = glygen_biomarker["biomarker_id"]
    assessed_entity_type = glygen_biomarker["assessed_entity_type"]
    # The A0001 record only has one value for the components[protein] array. Not
    # sure how/what it means when there are multiple values here. For this example,
    # just hardcoding for the assumption of one value.
    assessed_biomarker_entity = ""
    assessed_biomarker_entity_id = ""
    if "components" in glygen_biomarker:
        for mol in ["protein", "glycan"]:
            if mol in glygen_biomarker["components"]:
                if len(glygen_biomarker["components"][mol]) > 0:
                    o = glygen_biomarker["components"][mol][0]
                    assessed_biomarker_entity = o["name"] if "name" in o else assessed_biomarker_entity
                    
                    # In our data model, the ID value is formatted as "{resource acronym}:{accession}", so
                    # this value should be "UPKB:P05231-1". For different assessed entity type biomarkers
                    # we have to support multiple databases so the resource acronym/namespace is required.
                    # Not sure if you will have to support this.
                    assessed_biomarker_entity_id = o["accession"] if "accession" in o else assessed_biomarker_entity_id 

    citation_data = []
    # Convert publication data to Biomarker-Partnership citation data.
    for publication_object in glygen_biomarker["publication"]:
        citation_data.append(build_citation_object(publication_object))

    idx = 1
    # Loop through instances.
    for instance in glygen_biomarker["instances"]:
        # Build biomarker component.
        biomarker_component = build_biomarker_component(
            assessed_biomarker_entity,
            assessed_biomarker_entity_id,
            assessed_entity_type,
            instance,
        )
        # Build best biomarker role object.
        best_biomarker_role_object = build_best_biomarker_role_obejct(
            instance["best_biomarker_type"]
        )
        # Build condition object.
        condition_object = build_condition_object(instance["disease"])

        # Build the full biomarker record.
        result_data.append(
            build_record(
                biomarker_id + "." + str(idx),
                biomarker_component,
                best_biomarker_role_object,
                condition_object,
                citation_data,
            )
        )
        idx += 1

    return result_data



def build_biomarker_component(
    assessed_biomarker_entity: str,
    assessed_biomarker_entity_id: str,
    assessed_entity_type: str,
    instance: dict,
) -> list:
    """Takes in the GlyGen biomarker instance and builds the biomarker component in accordance
    with the Biomarker-Partnership data model. This assumes a non-panel/multi-component biomarker
    meaning that it returns just one biomarker component entry.

    Parameters
    ----------
    assessed_biomarker_entity : str
        The assessed biomarker entity.
    assessed_biomarker_entity_id : str
        The assessed biomarker entity id.
    assessed_entity_type : str
        The assessed entity type.
    instance : dict
        The GlyGen biomarker instance.

    Returns
    -------
    list
        The biomarker component according to the Biomarker-Partnership data model.
    """

    biomarker_component = {
        "biomarker": instance["status"],
        "assessed_biomarker_entity": {
            "recommended_name": assessed_biomarker_entity,
            # Not sure if the GlyGen data model supports synonyms,
            # if not then this can be left empty. When the Biomarker-Partnership
            # data is converted it hits the corresponding data resource API to
            # automatically pull in synonyms. Don't have to populate this if
            # this won't be used for searching.
            "synonyms": [],
        },
        "assessed_biomarker_entity_id": assessed_biomarker_entity_id,
        "assessed_entity_type": assessed_entity_type,
        "specimen": [
            {
                "name": instance["tissue"]["name"],
                # Similar to the note about the assessed_biomarker_entity_id in the
                # conversion.py script, our ID field includes the namespace. So in
                # our data model this would be "UBERON:0000178".
                "id": instance["tissue"]["id"],
                "name_space": instance["tissue"]["namespace"],
                "url": instance["tissue"]["url"],
                "loinc_code": instance["loinc_code"],
            }
        ],
        "evidence": build_evidence_sources(instance["evidence"]),
    }

    return [biomarker_component]


def build_evidence_sources(evidences: list) -> list:
    """Builds the evidence source list in accordance to the Biomarker-Partnership data model.
    We talked about the difference in philosophy with how the evidence source object is
    structured. The frontend will likely have to be slightly different based on the difference
    in structure. Make edits to this as you see fit.

    Parameters
    ----------
    evidence : list
        The evidence list from the GlyGen data model.

    Returns
    -------
    list
        The Biomarker-Partnership formatted evidence list.
    """
    evidence_sources: list = []

    for evidence in evidences:
        evidence_source = {
            "id": evidence["id"],
            "database": evidence["database"],
            "url": evidence["url"]
            # Not sure what to map here. GlyGen has the literature evidence field
            # but it is unclear to me if all evidence values have that direct
            # quote in the literature evidence field. Leaving this blank.
            #"evidence_list": [],
            # Not sure if you will be able to retroactively add evidence tags.
            #"tags": [],
        }
        evidence_sources.append(evidence_source)

    return evidence_sources


def build_condition_object(disease_instance: dict) -> dict:
    """Builds the condition object in accordance to the Biomarker-Partnership
    data model.

    Parameters
    ----------
    disease_instance : dict
        The instance disease object from the GLyGen data model.

    Returns
    -------
    dict
        The condition object formatted according to the Biomarker-Partnership
        data model.
    """

    # Build top level condition object.
    condition_object = {
        "id": disease_instance["disease_id"],
        "recommended_name": {
            "id": disease_instance["recommended_name"]["id"],
            "name": disease_instance["recommended_name"]["name"],
            "description": disease_instance["recommended_name"]["description"],
            "resource": disease_instance["recommended_name"]["resource"],
            "url": disease_instance["recommended_name"]["url"],
        },
        "synonyms": [],
    }

    # Loop through condition synonyms and add to synonym array.
    for synonym_entry in disease_instance["synonyms"]:
        condition_object["synonyms"].append(
            {
                "id": synonym_entry["id"],
                "name": synonym_entry["name"],
                "resource": synonym_entry["resource"],
                "url": synonym_entry["url"],
            }
        )

    return condition_object


def build_best_biomarker_role_obejct(role: str) -> list:
    """Takes the best biomarker role (called "type" in GLyGen) and builds the
    best_biomarker_role object in accordance to the Biomarker-Partnership
    data model.

    Parameters
    ----------
    role : str
        The best biomarker role string.

    Returns
    -------
    list
        The formatted best_biomarker_role object.
    """
    return [{"role": role}]


def build_record(
    biomarker_id: str,
    biomarker_component: list,
    best_biomarker_role_object: list,
    condition_object: dict,
    citation_data: list,
) -> dict:
    """Takes in the individual formatted components and builds the full
    biomarker record in accordance to the Biomarker-Partnership data model.

    Parameters
    ----------
    biomarker_id : str
        The biomarker ID.
    biomarker_component : list
        The biomarker component array.
    best_biomarker_role_object : list
        The best biomarker role list.
    condition_object : dict
        The condition object.
    citation_data : list
        The citation list.

    Returns
    -------
    dict
        The formatted biomarker record.
    """
    return {
        "biomarker_canonical_id": biomarker_id,
        "biomarker_component": biomarker_component,
        "best_biomarker_role": best_biomarker_role_object,
        "condition": condition_object,
        "publication": citation_data,
    }


def build_citation_object(publication_object: dict) -> dict:
    """Creates an instance of the citation list in the Biomarker-Partnership
    data model.

    Parameters
    ----------
    publication_object : dict
        The publication object from the GLyGen record.

    Returns
    -------
    dict
        An entry in the reference list.
    """
    return {
        "title": publication_object["title"],
        "journal": publication_object["journal"],
        "authors": publication_object["authors"],
        "date": publication_object["date"],
        "reference": [
            {"id": x["id"], "type": x["type"], "url": x["url"]}
            for x in publication_object["reference"]
        ],
        "evidence": [
            {"id": x["id"], "database": x["database"], "url": x["url"]}
            for x in publication_object["evidence"]
        ],
    }
