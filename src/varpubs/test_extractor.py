from varpubs.hgvs_extractor import extract_hgvsp_from_vcf


def main():
    vcf_path = "/Users/ilmaytas/Desktop/annotated.vcf"
    hgvsp_list = extract_hgvsp_from_vcf(vcf_path)

    for item in hgvsp_list:
        print(item)

    with open("output_pubmed_terms.txt", "w") as f:
        for item in hgvsp_list:
            f.write(item + "\n")


if __name__ == "__main__":
    main()
