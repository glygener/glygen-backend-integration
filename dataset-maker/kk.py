from pronto import Ontology
in_file = "downloads/mondo/current/mondo.owl"
out_file = "junk.obo"
edam = Ontology(in_file)
with open(out_file, "wb") as f:
    edam.dump(f, format="obo")


