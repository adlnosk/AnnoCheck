from Bio import Entrez
import sys

Entrez.email = "your.email@example.com"  # Replace with your email

def get_taxon_info(species_name):
    # Search for the species in NCBI Taxonomy
    search_handle = Entrez.esearch(db="taxonomy", term=species_name, retmode="xml")
    search_results = Entrez.read(search_handle)
    search_handle.close()

    if not search_results["IdList"]:
        print(f"Error: No Taxon ID found for {species_name}")
        return None, None, None

    taxon_id = search_results["IdList"][0]

    # Get the full lineage
    summary_handle = Entrez.esummary(db="taxonomy", id=taxon_id, retmode="xml")
    summary_results = Entrez.read(summary_handle)
    summary_handle.close()

    if "TaxId" not in summary_results[0] or "Lineage" not in summary_results[0]:
        print(f"Error: No lineage information found for {species_name}")
        return taxon_id, None, None

    lineage = summary_results[0]["Lineage"]

    # Extract family by splitting the lineage
    lineage_list = lineage.split("; ")
    family = lineage_list[-4] if len(lineage_list) >= 4 else "Unknown"

    return taxon_id, lineage, family

if __name__ == "__main__":
    species_name = sys.argv[1]
    taxon_id, lineage, family = get_taxon_info(species_name)

    if taxon_id:
        with open("taxon_id.txt", "w") as taxid_file:
            taxid_file.write(f"{taxon_id}\n")

    if lineage:
        with open("lineage.txt", "w") as lineage_file:
            lineage_file.write(f"{lineage}\n")

    if family:
        with open("family.txt", "w") as family_file:
            family_file.write(f"{family}\n")
