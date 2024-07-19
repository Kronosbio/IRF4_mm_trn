import os

from pipeline_tools.utils import utils
from pipeline_tools.utils.pipeline_utils import load_data_table


def make_gene_gffs(annot_file, gff_folder, species="HG18"):
    """Makes a tss gff with the given window size for all genes in the annotation.

    TSS +/-5kb,
    TSS +/-1kb,
    TSS +/-300,
    BODY +300/+3000.
    Can work on any genome build given the right annot file.
    """
    os.makedirs(gff_folder, exist_ok=True)

    start_dict = utils.make_start_dict(annot_file)
    gene_list = start_dict.keys()

    tss_loci = []
    for gene in gene_list:
        tss_locus = utils.Locus(
            start_dict[gene]["chr"],
            start_dict[gene]["start"][0] - 5000,
            start_dict[gene]["start"][0] + 5000,
            start_dict[gene]["sense"],
            gene,
        )
        tss_loci.append(tss_locus)
    tss_collection = utils.LocusCollection(tss_loci, 500)
    tss_gff_5kb = utils.locus_collection_to_gff(tss_collection)

    for gene in gene_list:
        tss_locus = utils.Locus(
            start_dict[gene]["chr"],
            start_dict[gene]["start"][0] - 1000,
            start_dict[gene]["start"][0] + 1000,
            start_dict[gene]["sense"],
            gene,
        )
        tss_loci.append(tss_locus)
    tss_collection = utils.LocusCollection(tss_loci, 500)
    tss_gff_1kb = utils.locus_collection_to_gff(tss_collection)

    tss_gff_300 = []
    txn_gff = []
    for line in tss_gff_5kb:
        gene = line[1]
        chrom = start_dict[gene]["chr"]
        start = start_dict[gene]["start"][0]
        end = start_dict[gene]["end"][0]
        sense = start_dict[gene]["sense"]
        tss_line = [chrom, gene, "", start - 300, start + 300, "", sense, "", gene]
        if sense == "+":
            txn_line = [chrom, gene, "", start + 300, end + 3000, "", sense, "", gene]
        else:
            txn_line = [chrom, gene, "", end - 3000, start - 300, "", sense, "", gene]

        tss_gff_300.append(tss_line)
        txn_gff.append(txn_line)

    utils.unparse_table(
        tss_gff_5kb,
        os.path.join(gff_folder, f"{species}_TSS_ALL_-5000_+5000.gff"),
        "\t",
    )
    utils.unparse_table(
        tss_gff_1kb,
        os.path.join(gff_folder, f"{species}_TSS_ALL_-1000_+1000.gff"),
        "\t",
    )
    utils.unparse_table(
        tss_gff_300, os.path.join(gff_folder, f"{species}_TSS_ALL_-300_+300.gff"), "\t"
    )
    utils.unparse_table(
        txn_gff, os.path.join(gff_folder, f"{species}_BODY_ALL_+300_+3000.gff"), "\t"
    )


def map_enriched_to_gff(
    data_file,
    set_name,
    gff_list,
    cell_type_list,
    enriched_folder,
    mapped_folder,
    macs=True,
    names_list=[],
    use_background=True,
):
    """
    Map enriched regions from a set of cell types to a set of gffs.

    Try to make a new folder for each gff.
    """
    data_dict = load_data_table(data_file)

    utils.format_folder(enriched_folder, True)
    utils.format_folder(mapped_folder, True)

    for gff_file in gff_list:
        gff_name = os.path.basename(gff_file).split(".")[0]
        print(f"making enriched regions to {gff_name}")

        outdir = os.path.join(mapped_folder, gff_name)
        os.makedirs(outdir, exist_ok=True)

        # first filter the name list
        cell_type_name_list = []
        if len(names_list) == 0:
            names_list = data_dict.keys()

        for name in names_list:
            # check to make sure in the right celltype
            # also make sure to not process WCEs
            if use_background and data_dict[name]["background"] == "NONE":
                continue
            cell_name = name.split("_")[0]
            if macs == True:
                if (
                    cell_type_list.count(cell_name) == 1
                    and data_dict[name]["enrichedMacs"] != "NONE"
                ):
                    cell_type_name_list.append(name)

            else:
                if (
                    cell_type_list.count(cell_name) == 1
                    and data_dict[name]["enriched"] != "NONE"
                ):
                    cell_type_name_list.append(name)

        cell_type_name_list.sort()

        mapped_gff = [["GFF_LINE", "ID"] + cell_type_name_list]
        # now we go through the gff and fill in stuff
        gff_table = utils.parse_table(gff_file, "\t")

        gff_loci = []
        # making the header line
        for line in gff_table:
            gff_locus = utils.Locus(line[0], line[3], line[4], line[6], line[8])
            gff_line = gff_locus.__str__()
            gff_id = line[1]

            gff_loci.append(gff_locus)
            mapped_gff.append([gff_line, gff_id])

        for name in cell_type_name_list:
            print(f"dataset {name}")
            if macs:
                enriched_collection = utils.import_bound_region(
                    os.path.join(enriched_folder, data_dict[name]["enrichedMacs"]), name
                )
            else:
                enriched_collection = utils.import_bound_region(
                    os.path.join(enriched_folder, data_dict[name]["enriched"]), name
                )
            for i in range(len(gff_loci)):
                if len(enriched_collection.get_overlap(gff_loci[i], "both")) > 0:
                    mapped_gff[i + 1].append(1)
                else:
                    mapped_gff[i + 1].append(0)

        utils.unparse_table(
            mapped_gff, os.path.join(outdir, f"{gff_name}_{set_name}.txt"), "\t"
        )


def make_gff_list_file(mapped_enriched_file, set_list, output, annot_file=""):
    """
    Create gene list file.

    set_list defines the dataset names to be used.
    AND operators within lists, OR operators outside of lists
    [[A,B],[C,D]] = (A AND B) OR (C AND D) for this row
    [[A],[B],[C],[D]] = A OR B OR C OR D for this row.
    """
    if len(annot_file) > 0:
        start_dict = utils.make_start_dict(annot_file)

    gene_list_file = []
    bound_gff_table = utils.parse_table(mapped_enriched_file, "\t")
    header = bound_gff_table[0]

    output_folder = os.path.dirname(output)
    os.makedirs(output_folder, exist_ok=True)

    # convert the set_list into column numbers
    column_set = []
    for binding_set in set_list:
        try:
            column_set.append([header.index(x) for x in binding_set])
        except ValueError:
            print("ERROR: not all datasets in binding table")
            exit()

    for i in range(1, len(bound_gff_table), 1):
        line = bound_gff_table[i]
        ref_id = line[1]

        # if any of these end up being true, the line gets added
        for and_columns in column_set:
            binding_vector = [int(line[x]) for x in and_columns]
            if ref_id == "NM_133941":
                print(binding_vector)
            if binding_vector.count(1) == len(binding_vector):
                if len(annot_file) > 0:
                    gene_list_file.append(
                        [
                            i,
                            bound_gff_table[i][1],
                            start_dict[bound_gff_table[i][1]]["name"],
                        ]
                    )
                else:
                    gene_list_file.append([i, bound_gff_table[i][1]])
                break
    print(len(gene_list_file))
    utils.unparse_table(gene_list_file, output, "\t")
